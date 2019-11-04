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
    %NOTICE: this may look redundant but Matlab's linspace is not very accurate and introduces small errors
    function out=linspace(l,u,n)
      out=l:(u-l)/(n-1):u;
      assert(num.ishomogeneous(out),...
        ['num.linspace failed to produce homoegeneous domain, min/max diff is ',num2str(diff(num.minmax(diff(out))))])
    end
    %NOTICE: unlike matlab's minmax, this returns the min/max of all entries
    function out=minmax(in)
      out=minmax(transpose(in(:)));
    end
    function out=ishomogeneous(x,tol)
      if ~exist('tol','var') || isempty(tol)
        tol=1e-15;
      end
      out=diff(num.minmax(diff(x)))<tol;
    end
    function io=makehomogeneous(io)
      if ~num.ishomogeneous(io)
        io=num.linspace(min(io),max(io),numel(io));
      end
    end
    function out=deltax(x)
      x=x(:);
      out=interp1(mean([x(1:end-1),x(2:end)],2),diff(x),x,'linear','extrap');
    end
    function M=diff_coeffs(dt,xpoints,nderiv)
      n=length(xpoints);
      i=0:n-1;
      w=warning('error','MATLAB:nearlySingularMatrix'); %#ok<NASGU,CTPCT>
      try 
        M=inv((xpoints(:)*ones(1,n)).^(ones(n,1)*i).*(dt*ones(n)).^(ones(n,1)*i)).*(factorial(0:n-1)'*ones(1,n));
      catch ME
        switch ME.identifier
        case 'MATLAB:nearlySingularMatrix'
          %reduce number of points
          if abs(xpoints(1))>=abs(xpoints(end))
            xpoints_new=xpoints(2:end);
          else
            xpoints_new=xpoints(1:end-1);
          end
          assert(numel(xpoints_new)>nderiv,'When trying to fix nearly singular matrix issue, ended up removing too many xpoints')
%           warning(['Using xpoints = [',strjoin(cellfun(@num2str,cells.m2c(xpoints),'Uniform',false),' '),'] for the ',...
%             str.th(nderiv),' derivative leads to nearly singular matrix; trying the following xpoints = [',...
%             strjoin(cellfun(@num2str,cells.m2c(xpoints_new),'Uniform',false),' '),']'])
           M=num.diff_coeffs(dt,xpoints_new,nderiv);
        otherwise
          error(['Cannot handle error ', ME.identifier,'; DEBUG NEEDED!'])
        end
        return
      end
      %negative sign needed here.
      M=M(nderiv+1,:)*(-1)^nderiv;
    end
    %applies convolution for each rows of 'in'
    function out=conv(in,coeffs)
      out=zeros(size(in));
      for i=1:size(in,1)
        out(i,:)=conv(in(i,:),coeffs,'same');
      end
    end
    %computes the derivative with order 'nderiv' of each row of 'f'
    %assuming homogeneous x-domain step equal to 'dt'
    %using a central stencial scheme with npoints*2+1
    function [df,out]=diff(f,dt,npoints,nderiv,varargin)
      %special calls
      if ischar(f)
        switch lower(f)
        case 'dt';       df=0.1;   return
        case 'npoints';  df=3;     return
        case 'nderiv';   df=1;     return
        case {'sin','p'};
          %need dt,npoints,nderiv
          out.dt=dt;
          out.npoints=npoints;
          out.nderiv=nderiv;
          switch lower(f)
          case 'sin'
            out.t=0:dt:dt*100*pi/2/dt;
            out.f=[sin(out.dt*out.t);cos(out.dt*out.t)];
            syms dts ts
            out.fa=[sin(dts*ts);cos(dts*ts)];
            out.dfa=double(subs(diff(out.fa,nderiv),{'dts','ts'},{dt,out.t}));
            out.dfn1=num.diff(out.f,dt,npoints,nderiv,varargin{:});
          case 'p'
            out.degree=2;
            %NOTICE: as of now, crossing 45 latitude creates problems
            out.t=sin(linspace(-pi*0.5,pi*0.5,round(2*pi/dt)));
            out.f=Ppq(out.t,out.degree,0);
            out.fa='legendre';
            out.dfa=gradient(out.f,1)./(ones(size(out.f,1),1)*transpose(num.deltax(out.t)));
            out.dfn1=Ppq(out.t,out.degree,1);out.dfn1=out.dfn1{2};
          end
          out.dfn2=num.diff(out.f,dt,npoints,nderiv,'fix_extremities');
          for i=1:size(out.f,1)
            delta=out.dfn1(i,:)-out.dfa(i,:);
            delta=delta(~isnan(delta));
            out.std1(i)=std(delta);
            out.max1(i)=max(abs(delta));
            delta=out.dfn2(i,:)-out.dfa(i,:);
            delta=delta(~isnan(delta));
            out.std2(i)=std(delta);
            out.max2(i)=max(abs(delta));
          end
          df=[];
          return
        case {'test-sin','test-p'}
          if ~exist('dt','var') || isempty(dt)
            dt=num.diff('dt');
          end
          if ~exist('nderiv','var') || isempty(nderiv)
            nderiv=num.diff('nderiv');
          end
          if ~exist('npoints','var') || isempty(npoints)
            npoints=num.diff('npoints');
          end
          [~,out]=num.diff(strrep(f,'test-',''),dt,npoints,nderiv,varargin{:});
          if ~cells.isincluded(varargin,'noplot')
            for i=1:size(out.dfn1,1)
              figure
              subplot(2,1,1)
              plot(out.t,out.dfn1(i,:)),hold on
              plot(out.t,out.dfa(i,:),'k')
              legend('numeric','analytic')
              ylabel(char(out.fa(i))),xlabel('ts')
              title([num2str(npoints),'-point stencil central ',num2str(nderiv),'-th derivative'])
              subplot(2,1,2)
              semilogy(abs(out.dfn1(i,:)-out.dfa(i,:)))
              title(['numeric - analytic, std=',num2str(out.std1(i)),', maxabs=',num2str(out.max1(i))])
            end
          end
          df=[];
          return
        case {'error-sin','error-p'}
          if ~exist('dt','var') || isempty(dt)
            dt=num.diff('dt');
          end
          if ~exist('nderiv','var') || isempty(nderiv)
            nderiv=num.diff('nderiv');
          end
          if ~exist('npoints','var') || isempty(npoints)
            min_npoints=max(1,nderiv);
            max_npoints=max([nderiv,num.diff('npoints')])*3;
            npoints=min_npoints:2:max_npoints;
          end
          stds1=zeros(size(npoints));
          maxs1=zeros(size(npoints));
          stds2=zeros(size(npoints));
          maxs2=zeros(size(npoints));
          for i=1:length(npoints)
            if i==1 || i==length(npoints)
              vararg={};
            else
              vararg={'noplot'};
            end
            [~,out]=num.diff(strrep(f,'error','test'),dt,npoints(i),nderiv,vararg{:},varargin{:});
            stds1(i)=out.std1(1)/diff(minmax(out.dfa(1,:)))*100;
            maxs1(i)=out.max1(1)/diff(minmax(out.dfa(1,:)))*100;
            stds2(i)=out.std2(1)/diff(minmax(out.dfa(1,:)))*100;
            maxs2(i)=out.max2(1)/diff(minmax(out.dfa(1,:)))*100;
          end
          plotting.figure;
          subplot(2,2,1)
          semilogy(npoints,stds1)
          xlabel('n')
          ylabel('std(error) [% of df]')
          plotting.enforce('plot_title',[num2str(nderiv),'-th order central difference'],'plot_legend_location','none');
          subplot(2,2,2)
          semilogy(npoints,maxs1)
          xlabel('n')
          ylabel('max(error) [% of df]')
          plotting.enforce('plot_title',[num2str(nderiv),'-th order central difference'],'plot_legend_location','none');
          subplot(2,2,3)
          semilogy(npoints,stds2)
          title('edges fixed')
          xlabel('n')
          ylabel('std(error) [% of df]')
          plotting.enforce('plot_title','edges fixed','plot_legend_location','none');
          subplot(2,2,4)
          semilogy(npoints,maxs2)
          xlabel('n')
          ylabel('max(error) [% of df]')
          plotting.enforce('plot_title','edges fixed','plot_legend_location','none');
          return
        end
      end
      %handle default inputs
      if ~exist('npoints','var') || isempty(npoints)
        npoints=num.diff('npoints');
      end
      if ~exist('nderiv','var') || isempty(nderiv)
        nderiv=num.diff('nderiv');
      end
      %handle shortcuts calls
      if npoints==0
        df=f;
        for i=1:nderiv
          df=gradient(df,dt);
        end
        return
      end
      %some sanity
      if isscalar(npoints)
        dpoints=transpose(-npoints:npoints);
      else
        assert(isvector(npoints),'input <npoints>, if not scalar, must be a 1D vector.')
        dpoints=npoints(:);
      end
      assert(nderiv <= npoints,['cannot calculate the ',num2str(nderiv),'-th derivative with only ',num2str(npoints),' points.'])
      nf=size(f,2);
      assert(nf>=npoints,['cannot use ',num2str(npoints),' points stencil with data with only ',num2str(nf),' points'])
      %compute derivative
      df=num.conv(f,num.diff_coeffs(dt,dpoints,nderiv));
      %find edges
      idx=[floor(npoints),ceil(npoints)];
      if cells.isincluded(varargin,'verbose')
        tab=16;
        disp(['start=',num2str(idx(1)),' end=',num2str(idx(2)),...
          ' dpoints=',num2str(dpoints(1)),':',num2str(dpoints(end))])
        for i=1:size(df,1)
          disp([' df(',num2str(i),',[1:',num2str(idx(1)),',',num2str(nf-idx(2)+1),':',num2str(nf),'])=',str.tablify(tab,df(i,[1:idx(1),end-idx(2)+1:end]))])
        end
      end
      %fixing warm-up/cool-down points, if requested
      if cells.isincluded(varargin,'crop_extremities')
        %crop it
        df(:,[1:idx(1),end-idx(2)+1:end])=[];
      elseif cells.isincluded(varargin,'nan_extremities')
        %nan it
        df(:,[1:idx(1),end-idx(2)+1:end])=NaN;
      elseif cells.isincluded(varargin,'fix_extremities')
        for i=0:npoints-1
          if cells.isincluded(varargin,'verbose')
            disp(['--- idx=',num2str(i+1)])
            disp(['dpoints=',num2str(-i),':',num2str(2*npoints-i)])
            disp(str.tablify(tab,['old df(:,',num2str(i+1),')'],transpose(df(:,i+1))))
          end
          tmp=num.conv(f(:,1:2*npoints+1),num.diff_coeffs(dt,-i:2*npoints-i,nderiv));
          df(:,i+1)=tmp(:,npoints+1);
          if cells.isincluded(varargin,'verbose')
            disp(str.tablify(tab,['fixed df(:,',num2str(i+1),')'],transpose(df(:,i+1))))
            disp(['dpoints=',num2str(-2*npoints+i),':',num2str(i)])
            disp(str.tablify(tab,['old df(:,',num2str(nf-i),')'],transpose(df(:,end-i))))
          end
          tmp=num.conv(f(:,end-2*npoints:end),num.diff_coeffs(dt,-2*npoints+i:i,nderiv));
          df(:,end-i)=tmp(:,npoints+1);
          if cells.isincluded(varargin,'verbose')
            disp(str.tablify(tab,['fixed df(:,',num2str(nf-i),')'],transpose(df(:,end-i))))
          end
        end
      end
    end
    function out=zeros(in)
      %also handes symbolic zeros
      switch class(in)
      case 'sym';   out=sym(zeros(size(in)));
      otherwise;    out=    zeros(size(in));
      end
    end
    function out=triu(in)
      switch ndims(in)
      case 1
        %this shouldn't happen, because matlab:
        %ndims([])=2
        %ndims(0)=2
        %ndims(ones(1,1))=2
        %ndims(ones(1,2))=2
        %ndims(ones(2,1))=2
        %ndims(ones(1,1,1))=2
        %ndims(ones(1,1,2))=3
      case 2
        if isvector(in)
          out=in;
        else
          out=triu(in);
        end
      case 3
        out=num.zeros(in);
        nk=size(in,3);
        for k=1:nk
          out(k:end,k:end,k)=num.triu(in(k:end,k:end,k));
        end
      otherwise
        error(['Cannot handle matrices with ',num2str(ndims(in)),' dimensions: implementation needed'])
      end
    end
    function out=symmetricidx(in)
      switch ndims(in)
      case 2
        if isvector(in)
          %trivial call, no symmetry
          out=reshape(1:numel(in),size(in));
          return
        else
          n1=size(in,1);
          n2=size(in,2);
          ci=cell(n1,n2);
          s1=str.characters(1:n1);
          s2=str.characters(1:n2);
          for i=1:n1
            for j=1:n2
              ci{i,j}=[s1(i),s2(j)];
            end
          end
        end
      case 3  
        n1=size(in,1);
        n2=size(in,2);
        n3=size(in,3);
        ci=cell(n1,n2,3);
        s1=str.characters(1:n1);
        s2=str.characters(1:n2);
        s3=str.characters(1:n3);
        for i=1:n1
          for j=1:n2
            for k=1:n3
              ci{i,j,k}=[s1(i),s2(j),s3(k)];
            end
          end
        end
      end
      out=num.zeros(ci);
      ci=cellfun(@sort,ci,'UniformOutput',false);
      done=false(size(ci));
      for i=1:numel(ci)
        if done(i), continue, end
        same_idx=strfind(ci,ci{i});
        same_idx=cellfun(@(i) ~isempty(i),same_idx);
        out(same_idx)=i;
        done(same_idx)=true;
      end
    end
    function out=makesymmetric(in)
      idx=num.symmetricidx(in);
      out=num.zeros(in);
      for i=1:numel(in)
        out(i)=in(idx(i));
      end
    end
    function out=issymmetric(in)
      idx=num.symmetricidx(in);
      idx_flat=unique(idx(:));
      for i=1:numel(idx_flat)
        f=(idx==idx_flat(i));
        if sum(f)==1, continue, end
        fv=in(f);
        if ~all(fv==fv(1))
          out=false;
          return
        end
      end
      out=true;
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
    %% symbolic stuff
    function [common,map,pow]=factorize(expr)
      %NOTE: expr = (map*comm).^pow
      %gather factors
      fact=cell(size(expr));
      for i=1:numel(expr)
        fact{i}=factor(expr(i));
      end
      %get common factors
      common=unique([fact{:}]);
      %get power
      pow=zeros(numel(expr),numel(common));
      for i=1:numel(expr)
        for j=1:numel(common)
          for l=1:numel(fact{i})
            if isequal(fact{i}(l),common(j))
              pow(i,j)=pow(i,j)+1;
            end
          end
        end
      end
      %get mapping of common factors
      map=(pow~=0);
    end
    %NOTE: exact inverse operation to factorize
    function expr=defactorize(common,map,pow)
      expr=(map*common).^pow;
    end
    function out=matvecmul(expr,var)
      %NOTE: expr=out*var (matrix multiplication)
      out=sym(zeros(size(expr)));
      %loop through the expressions
      for i=1:numel(expr)
        for j=1:numel(var)
          out(i,j)=diff(expr(i),var(j));
        end
      end
    end
    % Christoffel Symbol of the Second Kind: Lambda^m_{i,j}
    % http://mathworld.wolfram.com/ChristoffelSymboloftheSecondKind.html
    function out=cs2k(g,u,m,invg)
      n=numel(u);
      %this is only used internally to avoid re-computing the inverse when the full tensor is requested
      if ~exist('invg','var') || isempty(invg)
        invg=inv(g);
      end
      if ~exist('m','var') || isempty(m)
        %full tensor requested, loop over m
        out=sym(zeros(n,n,n));
        for m=1:n
          out(:,:,m)=num.cs2k(g,u,m,invg);
        end
      else
        out=sym(zeros(n));
        for i=1:n
          for j=1:n
            for k=1:n
              out(i,j)=out(i,j)+invg(k,m)*(diff(g(i,k),u(j))+diff(g(j,k),u(i))-diff(g(i,j),u(k)));
            end
          end
        end
        out=simplify(0.5*out);
      end
    end
    function out=latex_def(name,value,repstr,transpose_flag)
      %NOTE: defined a latex macro called <name> with the latex string of the symbolic expression <value>
      %NOTE: <repstr> is the replacement string cell array that defines how matlab's latex command is fine-tuned
      if exist('transpose_flag','var') && ~isempty(transpose_flag) && transpose_flag
        value=[latex(transpose(value)),'^T'];
      else
        value=latex(value);
      end
      if exist('repstr','var') && ~isempty(repstr)
        value=str.rep(['$',value,'$'],repstr{:});
      end
      out=['\def\',name,'{',value,'}'];
    end
    function out=dydx(y,x,Y)
      %NOTE: derivatives on the symbol in 'Y' are taken implicitly, i.e. dY/dx=Y_x
      %INFO: a 3D tensor of second order can be reprenented as:
      %  syms x y z V; reshape(num.dydx(num.dydx(V,[x;y;z],V),[x;y;z],V),3,3)
      assert(strcmp(class(y),'sym') && strcmp(class(x),'sym'),'Need symbolics') %#ok<STISA>
      if ~exist('Y','var')
        %patch something that will fail the contains test
        Y=sym(str.rand(numel(char(y))));
      end
      if isscalar(x) && isscalar(y)
        if contains(char(y),char(Y)) 
          out=sym([char(y),char(x)]);
        else
          out=diff(y,x);
        end
      elseif isscalar(y)
        out=cells.c2m(arrayfun(@(i) num.dydx(y,i,Y),x(:),'UniformOutput',false));
      elseif isscalar(x)
        out=cells.c2m(arrayfun(@(i) num.dydx(i,x,Y),y(:),'UniformOutput',false));
      else
        [XM,YM]=meshgrid(transpose(y(:)),x(:));
        out=cells.c2m(arrayfun(@(i,j) num.dydx(i,j,Y),XM(:),YM(:),'UniformOutput',false));
      end
    end
  end 
end
