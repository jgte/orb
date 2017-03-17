classdef num
  methods(Static)
    function out=cov(x,y)
      n=size(x,2);
      if n==0
        out=0;
        return
      end
      if any(size(x)~=size(y))
        out=Nan;
        return
      end
      out=size(1,n);
      for i=1:n
        out(i)=sum(...
          (x(:,i)-mean(x(:,i))).*...
          (y(:,i)-mean(y(:,i)))...
        )/size(x,1);
      end
    end
    function out=corrcoef(x,y)
      n=min([size(x,2),size(y,2)]);
      if n==0
        out=0;
        return
      end
      out=size(1,n);
      for i=1:n
        if any([std(x(:,i)),std(y(:,i))]==0)
          out(i)=0;
        else
          out(i)=num.cov(x(:,i),y(:,i))./std(x(:,i),1)./std(y(:,i),1);
        end
      end
    end
    function out=pardecomp(t,y,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't', @(i) isnumeric(i) && isvector(i) && size(i,2)==1          && ~any(isnan(i)) );
      p.addRequired( 'y', @(i) isnumeric(i) && isvector(i) && all(size(i)==size(t)) && ~any(isnan(i)) );
      p.addParameter('polynomial',[1 1],                                @(i) isnumeric(i) || isempty(i));
      p.addParameter('sinusoidal',[2*min(diff(t)),(t(end)-t(1))/2],     @(i) isnumeric(i) || isempty(i) || isduration(i));
      p.addParameter('phi',       [0,0],                                @(i) isnumeric(i) || isempty(i));
      p.addParameter('mode',      'struct',                             @(i) ischar(i));
      % parse it
      p.parse(t,y,varargin{:});
      % branch on mode
      %assign outputs
      switch p.Results.mode
        case 'test'
          n=1000;
          randn_scale=0.2;
          poly_scale=[0.3 0.8 1.6];
          sin_scale=[0.8 0.5];
          sin_T=n./[3 5];
          sin_phi=[pi/3 pi/4]; %randn(numel(sin_scale));
          %derived parameters
          legend_str=cell(1,numel(poly_scale)+numel(sin_scale));
          t=transpose(1:n);
          y=randn_scale*randn(n,1);
          for i=1:numel(poly_scale)
            legend_str{i}=['t^',num2str(i-1)];
            y=y+poly_scale(i)*(t/n).^(i-1);
          end
          for i=1:numel(sin_scale)
            legend_str{i+numel(poly_scale)}=['sin_',num2str(i)];
            y=y+sin_scale(i)*sin(2*pi/sin_T(i)*t+sin_phi(i));
          end
          %decompose
          a=num.pardecomp(t,y,'polynomial',ones(size(poly_scale)),'sinusoidal',sin_T,'mode','solve_phi');
          %show results
          figure
          plot(t,y,t,a.y_polynomial,t,a.y_sinusoidal,t,a.y_res)
          legend('original',legend_str{:},'res')
        case 'struct'
          %simpler names
          np=numel(p.Results.polynomial);
          ns=numel(p.Results.sinusoidal);
          ny=numel(y);
          %handle durations
          if isduration(p.Results.sinusoidal)
            T=seconds(p.Results.sinusoidal);
          else
            T=p.Results.sinusoidal;
          end
          %handle empty phi
          if isempty(p.Results.phi)
            phi=zeros(size(T));
          else
            phi=p.Results.phi;
          end
          % init design matrix
          A=zeros(ny,np+ns);
          S=zeros(np+ns);
          % build design matrix: polynomial coefficients
          for i=1:np
            S(i)=t(end)^(i-1);
            A(:,i)=t.^(i-1)/S(i);
          end
          % build design matrix: sinusoidal coefficients
          for i=np+(1:ns)
            S(i)=1;
            A(:,i)=sin(2*pi/T(i-np)*t+phi(i-np));
          end
          %solve the system of linear equations
          x=A\y.*S;
          %get modelled components
          y_mod=zeros(numel(y),numel(x));
          for j=1:numel(x)
            y_mod(:,j)=A(:,j)*x(j);
          end
          %get residuals
          y_res=y-sum(y_mod,2);
          %assign outputs
          out=struct(...
            'polynomial',x(1:np),...
            'sinusoidal',x(np+1:end),...
            'y_polynomial',y_mod(:,1:np),...
            'y_sinusoidal',y_mod(:,np+1:end),...
            'y_res',y_res,...
            'norm',norm(y_res)...
          );
        case 'solve_phi'
          %declare function to optimize
          fun=@(phi) (...
            num.pardecomp(t,y,...
              'polynomial',p.Results.polynomial,...
              'sinusoidal',p.Results.sinusoidal,...
              'phi',phi,...
              'mode','struct'...
            ).norm ...
          );
          %find the phase angles
          phi=fminsearch(fun,p.Results.phi);
          %get all info (the system is solved one more time, but with the correct phi)
          out=num.pardecomp(t,y,varargin{:},'phi',phi,'mode','struct');
          %append initial phase
          out.phi=phi(:);
      otherwise
        error([mfilename,': unknown output mode ''',p.Results.mode,'''.'])
      end
    end
    % shortcut for [s.field(:,column)]
    function out=struct_deal(s,field,line,column)
      if ~xor(isempty(line),isempty(column))
        error([mfilename,': (only) one of ''line'' or ''column'' must be empty.'])
      end
      if isempty(line)
        out=zeros(size(s(1).(field),1),numel(s));
        for i=1:numel(s)
          out(:,i)=s(i).(field)(:,column);
        end
      else
        out=zeros(numel(s),size(s(1).(field),2));
        for i=1:numel(s)
          out(i,:)=s(i).(field)(line,:);
        end
      end
    end
  end
end