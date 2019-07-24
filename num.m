classdef num
  methods(Static)
    %% utils
    % converts cell array to numeric array if all entries are numeric
    function inout=cell(inout)
      if iscell(inout) && all(cellfun(@isnumeric,inout))
        inout=cell2mat(inout);
      end
    end
    function out=odd(in)
      if ~isscalar(in)
        out=arrayfun(@num.odd,in);
        return
      end
      if mod(in,2)==0
        out=in-sign(in);
      else
        out=in;
      end
    end
    function out=even(in)
      if ~isscalar(in)
        out=arrayfun(@num.even,in);
        return
      end
      if mod(in,2)==0
        out=in;
      else
        out=in-sign(in);
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
    %removes NaNs from the extremeties of 'in'
    function out=trim_NaN(in)
      transpose_flag=false;
      %need column vector
      if isvector(in) && size(in,2)>1; transpose_flag=true; end
      %enforce transpose
      if transpose_flag; in=transpose(in); end
      %init loop
      out=true(size(in));
      in_isnan=isnan(in);
      for j=1:size(in,2)
        for i=1:size(in,1)
          if ~in_isnan(i,j); break; end
          out(i,j)=false;
        end
      end
      for j=size(in,2):-1:1
        for i=size(in,1):-1:1
          if ~in_isnan(i,j) || ~out(i,j); break; end
          out(i,j)=false;
        end
      end
    end
    %% pardecomp stuff (deprecated)
    function plot_pardecomp(pd_struct)
      s=200+[0,0,21,9]*50;
      figure('Position',s,'PaperPosition',s)
      legend_str=cell(1,numel(pd_struct.polynomial)+numel(pd_struct.sinusoidal)+numel(pd_struct.cosinusoidal)+2);
      plot(pd_struct.in.t,pd_struct.in.y,'b','LineWidth',2), hold on
      c=0;
      c=c+1;legend_str{c}='residual';
      plot(pd_struct.in.t,pd_struct.y_res,'k','LineWidth',2)
      c=c+1;legend_str{c}='original';
      for i=1:numel(pd_struct.polynomial)
        plot(pd_struct.in.t,pd_struct.y_polynomial(:,i),'r','LineWidth',2)
        c=c+1;legend_str{c}=['t^',num2str(i-1),':',num2str(pd_struct.polynomial(i))];
      end
      for i=1:numel(pd_struct.sinusoidal)
        plot(pd_struct.in.t,pd_struct.y_sinusoidal(:,i),'g','LineWidth',2)
        c=c+1;legend_str{c}=['sin_',num2str(i),':',num2str(pd_struct.sinusoidal(i))];
      end
      for i=1:numel(pd_struct.cosinusoidal)
        plot(pd_struct.in.t,pd_struct.y_cosinusoidal(:,i),'m','LineWidth',2)
        c=c+1;legend_str{c}=['cos_',num2str(i),':',num2str(pd_struct.cosinusoidal(i))];
      end
      legend(legend_str,'location','eastoutside')
      title(['norm(x+res-y)=',num2str(norm(sum([pd_struct.y_polynomial,pd_struct.y_sinusoidal,pd_struct.y_cosinusoidal,pd_struct.y_res,-pd_struct.in.y],2))),...
        newline,'T=',num2str(pd_struct.in.T(:)')])
      fs=16;
      set(    gca,          'FontSize',fs);
      set(get(gca,'Title' ),'FontSize',round(fs*1.3));
      set(get(gca,'XLabel'),'FontSize',round(fs*1.1));
      set(get(gca,'YLabel'),'FontSize',round(fs*1.2));
      grid on
    end
    function pd_struct=save_pardecomp(t,T,x,y,y_mod,y_res,np,ns)
      pd_struct=struct(...
        'in',struct(... %this is just record keeping
          't',t,...
          'y',y,...
          'T',T...
        ),...
        'polynomial',x(1:np),...          %constant, liner, quadratic, etc (reverse order than the poly* functions of matlab)
        'sinusoidal',x(np+1:np+ns),...    %in the same order as the periods defined in the 'sinusoidal' argument
        'cosinusoidal',x(np+ns+1:end),... %in the same order as the periods defined in the 'sinusoidal' argument
        'y_polynomial',y_mod(:,1:np),...
        'y_sinusoidal',y_mod(:,np+1:np+ns),...
        'y_cosinusoidal',y_mod(:,np+ns+1:end),...
        'y_res',y_res,...
        'rnorm',norm(y_res)./norm(y),...
        'norm',norm(y_res)...
      );
    end
    function out=pardecomp(t,y,varargin)
      %NOTICE: the timescale of t must be defined externally, it is implicit here
      if nargin==0
        num.pardecomp([0;1],[-1;-1],'mode','test');
        return
      end
% add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        {...
          't',         [], @(i) (isnumeric(i) && isvector(i) && all(size(i)==size(t)) && ~any(isnan(i))) || isempty(i) ;...
          'y',         [], @(i) (isnumeric(i) && isvector(i) && all(size(i)==size(t)) && ~any(isnan(i))) || isempty(i) ;...
          'polynomial', 2,                               @(i) (isnumeric(i) && isscalar(i)) || isempty(i);...
          'sinusoidal',[2*min(diff(t)),(t(end)-t(1))/2], @(i) isnumeric(i) || isempty(i);...
          't0',        [],                               @(i) isnumeric(i) && isscalar(i);...
          'mode',      'pd_struct',                      @(i) ischar(i);...
          'x',         [],                               @(i) isnumeric(i) || isempty(i);...
          'pd_struct', [],                               @(i) isstruct(i)  || isempty(i);...
        },...
      },varargin{:});
      %patch derived parameters
      if ~isempty(v.t) && isempty(v.t0)
        v.t0=t(1);
      end
      if ~isempty(v.pd_struct)
        if isempty(v.t); v.t=v.pd_struct.in.t; end
        if isempty(v.y); v.y=v.pd_struct.in.y; end
        if isempty(v.T); v.sinusoidal=v.pd_struct.in.T; end
        
%         continuar aqui
        
      end
      %some sanity
      if isempty(y)
        assert(strcmp(p.Results.mode,'model'),'if input ''y'' is empty, then mode must be ''model''')
        y=ones(size(t));
      end
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't', @(i)  isnumeric(i) && isvector(i) && size(i,2)==1          && ~any(isnan(i)) );
      p.addRequired( 'y', @(i) (isnumeric(i) && isvector(i) && all(size(i)==size(t)) && ~any(isnan(i))) || isempty(i) );
      %number of polynomial coefficients (not order): 1 constant, 2 constant+linear, 3 constant+linear+quadratic, etc
      p.addParameter('polynomial',2,                                @(i) (isnumeric(i) && isscalar(i)) || isempty(i)); 
      %sinusoidal periods (in implicit units):
      p.addParameter('sinusoidal',[2*min(diff(t)),(t(end)-t(1))/2], @(i) isnumeric(i) || isempty(i));
      p.addParameter('t0',        t(1),                             @(i) isnumeric(i) && isscalar(i));
      p.addParameter('mode',      'pd_struct',                      @(i) ischar(i));
      %these parameters are only valid for the "model" mode.
      p.addParameter('x',         [],                               @(i) isnumeric(i) || isempty(i));
      p.addParameter('pd_struct', [],                               @(i) isstruct(i)  || isempty(i));
      % parse it
      p.parse(t,y,varargin{:});

      %simpler names
      np=p.Results.polynomial;
      ns=numel(p.Results.sinusoidal);
      ny=numel(y);
      %start from t0
      t=t-p.Results.t0;
      % trivial call
      if ~any(y~=0)
        %assign outputs
        out=num.save_pardecomp(t,p.Results.sinusoidal,...
          zeros(np+2*ns,1),... x
          y,... y
          zeros(ny,np+2*ns),... y_mod
          zeros(ny,np+2*ns),... y_res
          np,ns...
        );
        %we're done
        return
      end
      % branch on mode
      switch p.Results.mode
        case 'test'
          %test parameters
          step=1;
          n=10000;
          poly_coeffs=[1 3 5]./[1 n n^2];
          sin_periods=n/step./[2 5];
          sin_periods_assumed=sin_periods+randn(size(sin_periods))*n/1e2;
           sin_coeffs=[0.5 3];
           cos_coeffs=[2 0.8];
          randn_scale=0.5;
          %inform
          disp(['sin_periods : ',num2str(sin_periods(:)')])
          disp(['poly_coeffs : ',num2str(poly_coeffs(:)')])
          disp(['sin_coeffs  : ',num2str(sin_coeffs(:)')])
          disp(['cos_coeffs  : ',num2str(cos_coeffs(:)')])
          disp(['randn_scale : ',num2str(randn_scale(:)')])
          %derived parameters
          t=transpose(1:step:(n*step));
          y=num.pardecomp(t,[],'mode','model',...
            'polynomial',numel(poly_coeffs),...
            'sinusoidal',sin_periods,...
            'x',[poly_coeffs,sin_coeffs,cos_coeffs]...
          );
          %sum all components
          y=sum(y,2);
          %add noise
          y=y+randn_scale*randn(size(y));
          %decompose
          pd_struct=num.pardecomp(t,y,'mode','pd_struct',...
            'polynomial',numel(poly_coeffs),...
            'sinusoidal',sin_periods_assumed...
          );
          %show results
          num.plot_pardecomp(pd_struct);
        % returns the design matrix
        case 'design'
          %get the period(s)
          if ~isempty(p.Results.pd_struct)
            %TODO
          end
          T=p.Results.sinusoidal;
          % init design matrix
          A=zeros(ny,np+2*ns);
          % build design matrix: polynomial coefficients
          for i=1:np
            A(:,i)=t.^(i-1);
          end
          % build design matrix: sinusoidal coefficients
          for i=np+(1:ns)
            A(:,i)=sin(2*pi/T(i-np)*t);
          end
          % build design matrix: co-sinusoidal coefficients
          for i=np+ns+(1:ns)
            A(:,i)=cos(2*pi/T(i-np-ns)*t);
          end
          %outputs
          out=A;
        % returns the y vector
        case 'model'
          %sanity
          assert(~isempty(p.Results.x),[mfilename,': need input argument ''x''.'])
          %get the design matrix
          A=num.pardecomp(t,y,varargin{:},'mode','design');
          %get modelled components
          x=p.Results.x;
          y_mod=zeros(numel(y),numel(x));
          for j=1:numel(x)
            y_mod(:,j)=A(:,j)*x(j);
          end
          %outputs
          out=y_mod;
        % returns struct with x coefficients and some more stuff
        case 'pd_struct'
          %get the design matrix
          A=num.pardecomp(t,y,varargin{:},'mode','design');
          %solve the system of linear equations
          x=A\y;
          %get modelled components
          y_mod=num.pardecomp(t,y,varargin{:},'mode','model','x',x);
          %get residuals
          y_res=y-sum(y_mod,2);
          %assign outputs
          out=num.save_pardecomp(t,p.Results.sinusoidal,x,y,y_mod,y_res,np,ns);
      otherwise
        error([mfilename,': unknown output mode ''',p.Results.mode,'''.'])
      end
    end
    function x_opt=param_search1(fun,x,varargin)
      % input parsing
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'fun',                   @(i) isa(i,'function_handle'));
      p.addRequired( 'xrange',                @(i) isnumeric(i));
      p.addParameter('searchspace', 'linear', @(i) ischar(i));
      p.addParameter('searchlen',numel(x)*10, @(i) isnumeric(i));
      p.addParameter('interpmethod','spline', @(i) ischar(i));
      p.addParameter('plot',           false, @(i) islogical(i));
      p.addParameter('enforce_scalar', false, @(i) islogical(i));
      p.addParameter('vectorize',     false, @(i) islogical(i));
      p.parse(fun,x,varargin{:})
      %define search space
      switch p.Results.searchspace
      case {'linear'}
        xd=linspace(x(1),x(end),p.Results.searchlen);
      case {'log'}
        xd=logspace(log10(x(1)),log10(x(end)),p.Results.searchlen);
      otherwise
        error(['Cannot handle argument ''searchspace'' with value ''',p.Results.searchspace,'''.'])
      end
      %use vectorize if allowed
      if p.Results.vectorize
        y=fun(x);
      else
        %init records
        y=nan(size(x));
        %loop over xrange
        for i=1:numel(x)
          y(i)=fun(x(i));
        end
      end
      %interpolate
      yd=interp1(x,y,xd,p.Results.interpmethod);
      %pick the minimum
      yi=yd==min(yd);
      x_opt=xd(yi);
      %warn user if minimum is at one extremety
      if yi(1) || yi(end)
        warning(['Could not find minimum of function ''',func2str(fun),''' in [',num2str(x(1)),',',num2str(x(end)),']'])
      end
      if p.Results.enforce_scalar
        x_opt=x_opt(1);
      end
      %debug plot
      if p.Results.plot
        figure
        plot(x_opt,yd(yi),'*','Markersize',10), hold on
        plot(x,y,'o')
        plot(xd,yd,'-')
        set(gca,'XScale',p.Results.searchspace)
        legend({'minimum','tries','interpolated'})
      end
    end
    function [x_opt,y,x]=param_brute(fun,x_lower,x_upper,varargin)
      % input parsing
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'fun',                   @(i) isa(i,'function_handle'));
      p.addRequired( 'x_upper',               @(i) isnumeric(i));
      p.addRequired( 'x_lower',               @(i) isnumeric(i));
      p.addParameter('searchspace', 'linear', @(i) ischar(i));
      p.addParameter('searchlen',        100, @(i) isnumeric(i));
      p.addParameter('searchlenmaxfactor',2e10,@(i) isnumeric(i));
      p.addParameter('searchiter',         0, @(i) isnumeric(i));
      p.addParameter('searchiter_counter', 0, @(i) isnumeric(i));
      p.addParameter('searchpdf',     'rand', @(i) ischar(i));
      p.addParameter('searchshrinkfactor', 2, @(i) isnumeric(i));
      p.addParameter('vectorize',      false, @(i) islogical(i));
      p.parse(fun,x_upper,x_lower,varargin{:})
      %sanity
      assert(numel(x_upper)==numel(x_lower),...
        'Arguments ''x0'', ''x_upper'' and ''x_lower'' must have the same number of elements')
      %simpler names
      l=p.Results.searchlen;
      n=numel(x_upper);
      searchpdf=str2func(p.Results.searchpdf);
      %define search space
      switch p.Results.searchspace
      case {'linear'}
        d0=x_lower;
        dx=x_upper-d0;
        x=searchpdf(l,n).*(ones(l,1)*dx)+ones(l,1)*d0;
      case {'log'}
        d0=log10(x_lower);
        dx=log10(x_upper)-d0;
        x=10.^(searchpdf(l,n).*(ones(l,1)*dx)+ones(l,1)*d0);
      otherwise
        error(['Cannot handle argument ''searchspace'' with value ''',p.Results.searchspace,'''.'])
      end
      %limit the search length: get an estimate of the memory used by fun
      fun_info=functions(fun);
      fun_mem=sum(structfun(@(i) numel(i),fun_info.workspace{1}));
      if l*fun_mem^2>p.Results.searchlenmaxfactor
        l_old=l;
        l=floor(p.Results.searchlenmaxfactor/fun_mem^2);
        str.say('stack_delta',1,['searchlen reset to ',num2str(l),...
          ' because searchlen*fun_mem>searchlenmaxfactor (',...
          num2str(l_old),'*',num2str(fun_mem),'>',num2str(p.Results.searchlenmaxfactor),')'])
      end
      %use vectorize if allowed
      if p.Results.vectorize
        y=fun(x);
      else
        %init records
        y=nan(size(x,1),1);
        %loop over xrange
        for i=1:size(x,1)
          y(i)=fun(x(i,:));
        end
      end
      assert(~any(isnan(y)),'Found NaNs in the output of fun')
      %retrieve lowest value from all y and get corresponding x
      x_opt=x(find(min(y)==y,1,'first'),:);
      %iterate if requested
      if p.Results.searchiter>0
        sd=1+p.Results.searchiter_counter;
        str.say('stack_delta',sd,str.tablify(16,'----- iter ----- ',p.Results.searchiter_counter+1))
        str.say('stack_delta',sd,str.tablify(16,'res',min(y)))
        str.say('stack_delta',sd,str.tablify(16,'res',min(y)))
        str.say('stack_delta',sd,str.tablify(16,'x_lower',x_lower))
        str.say('stack_delta',sd,str.tablify(16,'x_opt',x_opt))
        str.say('stack_delta',sd,str.tablify(16,'x_upper',x_upper))
        if p.Results.searchiter_counter<=p.Results.searchiter
          %define next iter search space amplitude
          dx=dx/p.Results.searchshrinkfactor;
          %define next iter search space
          switch p.Results.searchspace
          case {'linear'}
            x_lower=x_opt-dx/2;
            x_upper=x_opt+dx/2;
          case {'log'}
            x_lower=x_opt./10.^(dx/2);
            x_upper=x_opt.*10.^(dx/2);
          otherwise
            error(['Cannot handle argument ''searchspace'' with value ''',p.Results.searchspace,'''.'])
          end
          %iteration number
          iter=p.Results.searchiter_counter+1;
          %recursive call
          [~,y_iter,x_iter]=num.param_brute(fun,...
            x_lower,...
            x_upper,...
            varargin{:},...
            'searchlen',l,...
            'searchiter_counter',iter ...
          );
          %append to existing set of y,x
          y=[y(:);y_iter(:)];x=[x;x_iter];
          %update lowest value from all y and get corresponding x (may be the same)
          x_opt=x(find(min(y)==y,1,'first'),:);
        end
      end
    end
    %% memory utils (everything is in Mb)
    function out=memory_system
      if ismac
        com='sysctl hw.memsize';
        [status, cmdout]=system(com);
        if status; error(['could not issue the command ''',com,'''; debug needed.']);end
        cmdout=cells.rm_empty(strsplit(cmdout));
        out=str2double(cmdout{2})/1024^2;
      elseif isunis
        error('implementation needed')
      elseif ispc
        error('implementation needed')
      else
        error('cannot determine the operating system type')
      end
      if isnan(out);keyboard;end
    end
    function out=memory_matlab_fraction
      if ismac
        com='ps -c -o ''pmem comm'' | grep maci64';
        [status, cmdout]=system(com);
        if status; error(['could not issue the command ''',com,'''; debug needed.']);end
        cmdout=cells.rm_empty(strsplit(cmdout));
        out=str2double(cmdout{1})/100;
      elseif isunis
        error('implementation needed')
      elseif ispc
        error('implementation needed')
      else
        error('cannot determine the operating system type')
      end
      if isnan(out);keyboard;end
    end
    %% statistics
    function out = rms(x,w,dim)
      if ~exist('dim','var') || isempty(dim)
        dim=0;
      else
        assert(dim <= numel(size(w)),'the value of input <dim> must not be larger than the number of dimensions in <x>')
      end
      if ~exist('w','var') || isempty(w)
        w=ones(size(x));
      else
        assert(all(size(w) == size(x)),'inputs <x> and <w> must have the same size')
      end
      switch dim
      case 0
        %filtering NaNs
        i = ~isnan(x) & ~isnan(w);
        %compute
        out= sqrt(sum(w(i).*(x(i).^2))/sum(w(i)));
      case 1 %computes the RMS along columns, so a row is returned
        out=arrayfun(@(i) num.rms(x(:,i),w(:,i),0),1:size(x,2));
      case 2 %computes the RMS along rows, so a column is returned
        out=transpose(arrayfun(@(i) num.rms(x(i,:),w(i,:),0),1:size(x,1)));
      otherwise
        error(['dim=',num2str(dim),' not implemented.'])
      end
    end
    function out = mean(x,w,dim)
      if ~exist('dim','var') || isempty(dim)
        dim=0;
      else
        assert(dim <= numel(size(w)),'the value of input <dim> must not be larger than the number of dimensions in <x>')
      end
      if ~exist('w','var') || isempty(w)
        w=ones(size(x));
      else
        assert(all(size(w) == size(x)),'inputs <x> and <w> must have the same size')
      end
      switch dim
      case 0
        %filtering NaNs
        i = ~isnan(x) & ~isnan(w);
        %compute
        out= sum(w(i).*(x(i)))/sum(w(i));
      case 1 %computes the RMS along columns, so a row is returned
        out=arrayfun(@(i) num.mean(x(:,i),w(:,i),0),1:size(x,2));
      case 2 %computes the RMS along rows, so a column is returned
        out=transpose(arrayfun(@(i) num.mean(x(i,:),w(i,:),0),1:size(x,1)));
      otherwise
        error(['dim=',num2str(dim),' not implemented.'])
      end
    end
    function out = std(x,w,dim)
      out=sqrt(num.rms(x,w,dim).^2-num.mean(x,w,dim).^2);
    end
    function out=cov(x,y)
      n=size(x,2);
      if n==0
        out=0;
        return
      end
      if any(size(x)~=size(y))
        out=NaN;
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
  end
end
