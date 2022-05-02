classdef pardecomp
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list={...
      'debug',     false,               @(i) islogical(i) && isscalar(i);...
      'peekwidth', 10,       @num.isscalar;...
    };
  end
  properties(GetAccess=public,SetAccess=public)
    t  %time domain
    y  %data vector
    T  %periods
    np %polynomidal order
    t0 %zero-epoch
    xi %regression coefficients
    peekwidth
    Qy %observation/data covariance matrix (a.k.a. Cd)
  end
  properties(Dependent)
    ns %number of sine and cosine coefficients
    ny %data length
    nx %total number of coefficients
    x %regression coefficients as set/get method
    p %polynomial coefficients
    s %sinusoidal coefficients
    c %co-sinusoidal coefficients
    y_comp %all contributions in separate columns
    y_sum  %sum of all contributions
    yp %polynomial contribution
    ys %sinusoidal contribution
    yc %co-sinusoidal contribution
    yr %unmodelled residual
    rrn %relative norm of the residuals
    rn %norm of the residuals
    A  %design matrix
    AtQy
    N  %normal matrix
    h  %right-hand vector
    Qx %variance-covariance
    e  %postfit residuals
    Pa %projection matrix
  end
  methods(Static)
    %% object static methods
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(pardecomp.parameter_list); end
      out=v.picker(varargin{:});
    end
    function out=valid_t(x)
      out=isnumeric(x) && ~isempty(x) && isvector(x);
    end
    function out=valid_y(y)
      out=isnumeric(y) || isempty(y);
    end
    function out=xlength(np,T)
      assert(isnumeric(T)  || isempty(T),  'Ilegal input T')
      assert(isnumeric(np) && isscalar(np),'Ilegal input np')
      out=np+2*numel(T);
    end
    function out=xnames(np,T)
      out=cell(1,pardecomp.xlength(np,T));
      for i=1:np
        out{i}=['p',num2str(i-1)];
      end
      for i=1:numel(T)
        out{i+np}=['s',num2str(i)];
      end
      for i=1:numel(T)
        out{i+np+numel(T)}=['c',num2str(i)];
      end
    end
    function out=xstr(xname)
      switch xname(1)
        case 'p'
          switch xname(2)
            case '0'; out='bias';
            case '1'; out='trend';
            case '2'; out='quadratic';
            case '3'; out='cubic';
            otherwise; out=[xname(2),'-th order polynomial'];
          end
        case 's'; out='sine';
        case 'c'; out='cosine';
        otherwise
          error(['Cannot handle xname with value ''',xname,'''.'])
      end
    end
    function out=idx(np,T,name)
      assert(isnumeric(T)  || isempty(T),  'Ilegal input T')
      assert(isnumeric(np) && isscalar(np),'Ilegal input np')
      if iscellstr(name)
        %NOTICE: if the error 'Non-scalar in Uniform output, at index 3, output 1.' is triggered,
        %        it means one of the entries in 'name' is not a valid x-name
        out=cellfun(@(i) pardecomp.idx(np,T,i),name,'UniformOutput',true);
      else
        out=find(cells.isstrequal(pardecomp.xnames(np,T),name));
      end
    end
    function out=to_timescaled(t,timescale)
      out=time.num2duration(simpletimeseries.timescale(t),timescale);
    end
    function out=from_timescaled(t,timescale)
      out=simpletimeseries.timescale(time.duration2num(t,timescale));
    end
    %% decomposition into pd_set
    function pd_set=split(obj,varargin)
      % NOTICE: the following parameters are pretty much mandatory to define the parametric regression (see below)
      % 'T',[2*min(diff(t)),(t(end)-t(1))/2], @(i) isnumeric(i) || isempty(i);...
      % 'np',0,                               @num.isscalar;...
      % NOTICE: this method expect obj to be simpledata-esque
      v=varargs.wrap('sources',{....
        {...
          'timescale','seconds', @ischar;...
          'quiet',        false, @islogical;...
          'parallel',     false, @islogical;...
        },...
      },varargin{:});
      %need time series
      assert(ismethod(obj,'t_masked'),['Need ''simpletimeseries'' or derived object, not ''',class(obj),'''.'])
      %init loop vars
      t_pd=pardecomp.to_timescaled(obj.x_masked,v.timescale);
      y_pd=obj.y_masked;

      % TODO: pardecomp().lsq may well handle y_pd as an array

      % parfor wide objects
      % TODO: can't seem to get this to work: Error using pardecomp (line 313) Not enough input arguments.)
      if v.parallel
        parfor i=1:obj.width
          %create temporary var
          varargin_now=varargin;
          % call mother routine, solve and implement higher x-domain resolution
          % for the model if needed (trivial check done inside)
          % TODO: the latter needs checking
          d(i)=pardecomp(t_pd,y_pd(:,i),varargin_now{:});
          d(i)=d(i).lsq;
        end
      else
        s.msg=['Parametric decomposition of ',obj.descriptor]; s.n=obj.width;
        for i=1:obj.width
          % call mother routine, solve and implement higher x-domain resolution
          % for the model if needed (trivial check done inside)
          % TODO: the latter needs checking
          d(i)=pardecomp(...
            t_pd,y_pd(:,i),varargin{:}...
          ).lsq; %#ok<AGROW>
          if ~v.quiet; s=time.progress(s,i); end
        end
      end
      %init containers
      init=str2func(class(obj)); %use correct constructor
      pd_args=cell(1,4*pardecomp.xlength(d(1).np,d(1).T)); c=0;
      coeffnames=pardecomp.xnames(d(1).np,d(1).T);
      for j=1:numel(coeffnames)
        %retrieve coefficient index within its type
        switch coeffnames{j}(1)
        case 'p'
          i=str2double(coeffnames{j}(2:end))+1;
          labels=[str.th(i-1),'-order polynomial term'];
          units=[obj.units{1},'/',v.timescale,'^',num2str(i-1)];
        case {'c','s'}
          i=str2double(coeffnames{j}(2:end));
          labels=['sine term for period with ',num2str(d(1).T(i)),' ',v.timescale];
          if coeffnames{j}(1)=='c'; labels=['co',labels]; end %#ok<AGROW>
          units=obj.units{1};
        otherwise
          error(['Cannot understand the coefficient name ''',coeffnames{j},'''.'])
        end
        %save coefficients (indexed to zero-date, since they are time-invariant)
        o=init(time.zero_date,transpose(num.struct_deal(d,coeffnames{j}(1),i,[])),obj.varargin{:});
        o.descriptor=[coeffnames{j},' of ',str.clean(obj.descriptor,'file')];
        o.units(:)={units}; o.labels(:)={labels};
        %append to pd_args
        pd_args{c+1}=coeffnames{j};
        pd_args{c+2}=o;
        c=c+2;
        %reconstruct the time domain in the default timescale
        x_now=pardecomp.from_timescaled(d(1).t,v.timescale);
        %make sure time domain is not borked
        assert(all(simpledata.isx('==',obj.x_masked,x_now)),['time domain discrepancy in ',o.descriptor])
        %save timeseries represented by each coefficient
        o=init(obj.t_masked,num.struct_deal(d,['y',coeffnames{j}(1)],[],i),obj.varargin{:});
        o.descriptor=['p',num2str(i-1),' of ',str.clean(obj.descriptor,'file')];
        %restore gaps
        o=o.t_merge(obj.t);
        %append to pd_args
        pd_args{c+1}=['ts_',coeffnames{j}];
        pd_args{c+2}=o;
        c=c+2;
      end
      %save everything into pd_set record
      pd_set=pardecomp.assemble(...
        'T',d(1).T,...,
        'np',d(1).np,...
        't0',d(1).t0,...
        't',d(1).t,...   %this is the abcissae used in the regression (actually: t-t0)
        'time',obj.t,... %this is the abcissae defined in obj, including gaps
        'epoch',obj.epoch,...
        'descriptor',obj.descriptor,...
        'timescale',v.timescale,...
        pd_args{:});
      %this is the abcissae defined in obj, excluding gaps
      pd_set.t_masked=obj.t_masked;
      %save residuals
      o=init(obj.t_masked,num.struct_deal(d,'yr',[],1),obj.varargin{:});
      o.descriptor=['residual of ',str.clean(obj.descriptor,'file')];
      %restore gaps
      o=o.t_merge(obj.t);
      pd_set.res=o;
      %save norms
      o=init(time.zero_date,num.struct_deal(d,'rn',[],1),obj.varargin{:});
      o.descriptor=['norm of the residuals of ',str.clean(obj.descriptor,'file')];
      pd_set.norm=o;
      %save norm ratio
      o=init(time.zero_date,num.struct_deal(d,'rrn',[],1),obj.varargin{:});
      o.descriptor=['signal and residual norms ratio for ',str.clean(obj.descriptor,'file')];
      pd_set.rnorm=o;
    end
    %% reconstruction from pd_set
    function obj=join(pd_set,varargin)
      %sanity
      assert(pardecomp.ispd_set(pd_set),'Input ''pd_set'' must validate the pardecomp.ispd_set method')
      %easiner names
      coeffnamall=pardecomp.xnames(pd_set.np,pd_set.T);
      %parse inputs
      v=varargs.wrap('sources',{....
        {...
          'time',      pd_set.time, @isdatetime;...
          'coeffnames',coeffnamall,    @iscellstr;...
          'add_res',   false,       @islogical;...
        },...
      },varargin{:});
      %init output and copy metadata
      pd_set.metadata=varargs(pd_set.metadata).rm_empty;
      pd_set.varargin=pd_set.metadata.varargin;
      switch func2str(pd_set.init)
        case 'gravity';obj=eval([func2str(pd_set.init),'.unit(gravity.width2lmax(pd_set.width),''t'',v.time,pd_set.varargin{:},''scale'',0);']);
        otherwise;     obj=eval([func2str(pd_set.init),'.zero(v.time,pd_set.width,pd_set.varargin{:});']);
      end
      %match epochs (very important and not done with copying the metadata)
      obj.epoch=pd_set.epoch;
      %update descriptor (not done with copying the metadata1)
      obj.descriptor=pd_set.descriptor;
      %initialize coefficient containers
      coeffval=zeros(pardecomp.xlength(pd_set.np,pd_set.T),pd_set.width);
      coeffidx=zeros(pardecomp.xlength(pd_set.np,pd_set.T),1);
      % go over all coefficients and collect that component (if possible)
      for i=1:numel(v.coeffnames)
        %get time series name
        tname=['ts_',v.coeffnames{i}];
        %check if this time series is available and is in the same time domain
        if isfield(pd_set,tname) && pd_set.(tname).isxequal(v.time)
          %if so, just increment it
          obj=obj+pd_set.(tname);
        %check if this parameter is available
        elseif isfield(pd_set,v.coeffnames{i})
          y_now=pd_set.(v.coeffnames{i}).y;
          if size(y_now,1)>1
            assert(all(mean(y_now(2:end,:))==y_now(1,:)),[v.coeffnames{i},' has time-varying values, this is ilegal'])
          end
          %get parameter index
          idx=pardecomp.idx(pd_set.np,pd_set.T,v.coeffnames{i});
          %save the parameter values
          coeffval(idx,:)=y_now(1,:);
          coeffidx(idx)=1;
        end
      end
      %rebuild component timeseries if needed
      if any(coeffidx==1)
        % call mother routine
        s.msg=['Parametric reconstruction of ',obj.descriptor,' with components ',strjoin(coeffnamall(coeffidx==1),', ')]; s.n=obj.width;
        for i=1:obj.width
          obj=obj.set_cols(i,pardecomp(...
            pardecomp.to_timescaled(obj.x,pd_set.timescale),[],'T',pd_set.T,'np',pd_set.np,'t0',pd_set.t0,'x',coeffval(:,i)...
          ).y_sum);
          s=time.progress(s,i);
        end
      end
      %add residual, if available and in the same time domain
      if isfield(pd_set,'res') && v.add_res
        obj=obj+pd_set.res.interp(v.time);
      end
    end
    function out=default
      out={...
        't',                     [], @(i)  isnumeric(i) || isempty(i);...
        'T',                     [], @(i)  isnumeric(i) || isempty(i);...
        'np',                     0, @num.isscalar;...
        't0',                     0, @num.isscalar;...
        'epoch',     time.zero_date, @(i)  isdatetime(i) && isscalar(i);...
        'timescale',      'seconds', @ischar;...
        'descriptor',            '', @ischar;...
      };
    end
    function out=str(varargin)
      switch nargin
        case 1
          out=varargs(varargin{1}).str;
        otherwise
          out=varargs.wrap('sources',{pardecomp.default},varargin{:}).str;
      end
    end
    function pd_set=assemble(varargin)
      %NOTICE; this is needed so that the coeffnames get into v
      v=varargs(varargin);
      v=varargs.wrap('sources',{v,pardecomp.default},varargin{:});
      v=varargs.wrap('sources',{v,....
        {...
          'time',v.epoch+pardecomp.from_timescaled(v.t,v.timescale), @isdatetime;...
        },...
      },varargin{:});
      %initialize output
      pd_set=struct('t',v.t,'T',v.T,'np',v.np,'t0',v.t0,'time',v.time);
      %init looping
      records=struct([]);
      rnames={'width','class','metadata'};
      coeffnames=pardecomp.xnames(pd_set.np,pd_set.T);
      for i=1:numel(coeffnames)
        %check if this field exists
        if v.isparameter(coeffnames{i})
          d=v.(coeffnames{i});
          %check that all objects have the same character
          for f=1:numel(rnames)
            try
              rvalue=d.(rnames{f});
            catch
              rvalue=eval([rnames{f},'(d)']);
            end
            %need to delete the some metadata entries (it changes for every coefficient or type of coefficient)
            switch rnames{f}
            case 'metadata'
              for m={'descriptor','units','x_units','labels'}
                rvalue=rvalue.delete(m{1});
              end
            end
            if ~isfield(records,rnames{f})
              records(1).(rnames{f})=rvalue;
            else
              assert(cells.isequal(records.(rnames{f}),rvalue),...
                ['Conflict in ',rnames{f},' between the coefficient objects.'])
            end
          end
          pd_set.(coeffnames{i})=d;
          %these coefficients are static, so set the time domain to zero
          pd_set.(coeffnames{i}).t=time.zero_date;
        end
        %check if the timeseries of this field exists
        if v.isparameter(['ts_',coeffnames{i}])
          pd_set.(['ts_',coeffnames{i}])=v.(['ts_',coeffnames{i}]);
          assert(pd_set.(['ts_',coeffnames{i}]).istequal(v.time),'Time domain discrepancy: debug needed!')
        end
      end
      pd_set.init=str2func(records.class);
      pd_set.width=records.width;
      pd_set.metadata=records.metadata.rm_empty;
      pd_set.descriptor=v.descriptor;
      pd_set.start=v.time(1);
      pd_set.stop=v.time(end);
      pd_set.timescale=v.timescale;
      pd_set.epoch=v.epoch;
    end
    %% checking pd_set
    function out=ispd_set(in)
      out=false;
      for i={'t','T','np','t0','init','width','metadata','descriptor','start','stop','time'}
        if ~isfield(in,i{1}); disp(['ERROR: field ',i{1},' missing from pd_set']);return; end
      end
      out=true;
    end
    %% data visualization
    function out=table(pd_set,varargin)
      v=varargs.wrap('sources',{...
        {...
          'tabw',          20, @num.isscalar;...
          'time_unit', @years, @(i) isa(i,'function_handle');...
          'tablify',     true, @islogical;...
          'latex_table',false, @islogical;...
          'cols',           1, @isnumeric;...
        }},...
      varargin{:});
      ps=0;pc=0;c=0;
      if v.tablify
        out=cell(pd_set.np+numel(pd_set.T),1);
      else
        out=cell(pd_set.np+numel(pd_set.T),2+numel(v.cols));
      end
      if v.latex_table
        assert(~v.tablify,'options ''tablify'' and ''latex_table'' cannot both be true')
        c=c+1;
        out(c,1:2)={'Component',['Period [',func2str(v.time_unit),']']};
        for colsi=1:numel(v.cols)
%           out{c,colsi+2}=['Value ',num2str(v.cols(colsi))];
          out{c,colsi+2}='Value';
        end
      end
      for i=pardecomp.xnames(pd_set.np,pd_set.T)
        c=c+1;
        s=pardecomp.xstr(i{1});
        name=s;
        switch s
          case 'sine'
            ps=ps+1;
            sec=seconds(pardecomp.from_timescaled( pd_set.T(ps),pd_set.timescale ));
            period=num2str(v.time_unit(sec));
          case 'cosine'
            pc=pc+1;
            sec=seconds(pardecomp.from_timescaled( pd_set.T(pc),pd_set.timescale ));
            period=num2str(v.time_unit(sec));
          otherwise
            period='-';
        end
        if v.tablify
          out{c}=str.tablify(v.tabw,name,period,pd_set.(i{1}).y(v.cols));
        else
          out(c,1:2)={name,period};
          for colsi=1:numel(v.cols)
            out{c,colsi+2}=pd_set.(i{1}).y(v.cols(colsi));
          end
        end
      end
      if v.latex_table
        out=str.latex_table(out);
      end
    end
    %% general test for the current object
    function obj=test(varargin)
        v=varargs.wrap('sources',{...
        {...
          'plot',   false, @islogical;...
          'print',   true, @islogical;...
        }},...
      varargin{:});
      %test parameters
      step=1;
      n=10000;
      poly_coeffs=[1 3 5]./[1 n n^2];
      sin_periods=n/step./[2 5];
      sin_periods_assumed=sin_periods;
       sin_coeffs=[0.5 3];
       cos_coeffs=[2 0.8];
      %TODO: implement multiple columns in Y
%       randn_scale=[0.1,1,10];
      randn_scale=0.1;
      %derived parameters
      t=transpose(1:step:(n*step));
      %forward modelling
      ref=pardecomp(t,[],...
        'np',numel(poly_coeffs),...
        'T',sin_periods,...
        'p',poly_coeffs,...
        's',sin_coeffs,...
        'c',cos_coeffs...
      );
      y=ref.y_sum;
      %add noise
      y=y*ones(size(randn_scale))+randn(size(y))*randn_scale;
      %inversion
      obj=pardecomp(t,y,...
        'np',numel(poly_coeffs),...
        'T',sin_periods_assumed...
      ).lsq;
      %show results
      if v.plot
        for i=1:numel(randn_scale)
          obj.plot(varargin{:},'columns',{i});
        end
      end
      %inform
      if v.print
        obj.print([],ref)
      end
    end
  end
  methods
    function obj=pardecomp(t,y,varargin)
      p=machinery.inputParser;
      p.addRequired( 't', @(i) pardecomp.valid_t(i));
      p.addRequired( 'y', @(i) pardecomp.valid_y(i));
      p.addParameter('Qy',[],@(i) @isnumeric);
      [~,~,obj]=varargs.wrap('sinks',{obj},'parser',p,'sources',{....
        pardecomp.parameters('obj'),...
        {...
          'T',[2*min(diff(t)),(t(end)-t(1))/2], @(i) isnumeric(i) || isempty(i);...
          'np',0,   @num.isscalar;...
          't0',t(1),@num.isscalar;...
          'p', [],  @(i) (isnumeric(i) && isvector(i)) || isempty(i);...
          's', [],  @(i) (isnumeric(i) && isvector(i)) || isempty(i);...
          'c', [],  @(i) (isnumeric(i) && isvector(i)) || isempty(i);...
          'x', [],  @(i) (isnumeric(i) && isvector(i)) || isempty(i);...
        },...
      },'mandatory',{t,y},varargin{:});
      %need to assign the mandatory arguments
      obj.t=t;
      obj.y=y;
      %patch missing information
      if isempty(obj.Qy)
        %WARNING: this is a ny-by-ny matrix (i.e. potentially huge)
        obj.Qy=eye(obj.ny);
      end
    end
    %% info functions
    function disp_field(obj,field,tab,value,label,fmt)
      if ~exist('value','var') || isempty(value)
        value=obj.(field);
        if isnumeric(value) || iscell(value)
          value=value(1:min([numel(value),obj.peekwidth]));
        end
      end
      if ~exist('fmt','var')
        fmt='';
      end
      if ~exist('label','var') || isempty(label)
        disp([str.tabbed(field,tab),' : ',str.show(transpose(value(:)),fmt)])
      else
        disp([str.tabbed(label,tab),' : ',str.show(transpose(value(:)),fmt)])
      end
    end
    function print(obj,tab,ref)
      if ~exist('tab','var') || isempty(tab)
        tab=10;
      end
      if exist('ref','var')
        flag=true;
      else
        flag=false;
      end
      fmtf=['%',num2str(tab),'.',num2str(tab-2),'f '];
      fmte=['%',num2str(tab),'.',num2str(tab-7),'e '];
               obj.disp_field('T'      ,tab,[]               ,''     ,fmtf);
      if flag; ref.disp_field('T'      ,tab,[]               ,'T ref',fmtf);
               obj.disp_field('T delta',tab,obj.T(:)-ref.T(:),''     ,fmte); end
               obj.disp_field('p'      ,tab,[]               ,''     ,fmtf);
      if flag; ref.disp_field('p'      ,tab,[]               ,'p ref',fmtf)
               obj.disp_field('p delta',tab,obj.p(:)-ref.p(:),''     ,fmte); end
               obj.disp_field('s'      ,tab,[]               ,''     ,fmtf);
      if flag; ref.disp_field('s'      ,tab,[]               ,'s ref',fmtf)
               obj.disp_field('s delta',tab,obj.s(:)-ref.s(:),''     ,fmte); end
               obj.disp_field('c'      ,tab,[]               ,''     ,fmtf);
      if flag; ref.disp_field('c'      ,tab,[]               ,'c ref',fmtf)
               obj.disp_field('c delta',tab,obj.c(:)-ref.c(:),''     ,fmte); end
      obj.disp_field('nx',tab)
      obj.disp_field('ny',tab)
    end
    %% length functions
    function out=get.ns(obj); out=numel(obj.T);  end
    function out=get.nx(obj); out=obj.np+2*obj.ns; end
    function out=get.ny(obj)
      if isempty(obj.y)
        assert(~isempty(obj.t),'Both y and t are empty, this is ilegal')
        out=numel(obj.t);
      else
        out=size(obj.y,1);
      end
    end
    %% coeffcients idx functions
    function idx=p_idx(obj);  idx=1:obj.np;  end
    function idx=s_idx(obj);  idx=obj.np+1:obj.np+obj.ns;  end
    function idx=c_idx(obj);  idx=obj.np+obj.ns+1:obj.np+2*obj.ns;  end
    %% coeffcients get functions
    function out=get.x(obj)
      if isempty(obj.xi)
        assert(~isempty(obj.p) || ~isempty(obj.s),'Cannot determine the x-vector because no coefficient record is defined')
        obj.xi=[obj.p(:);obj.s(:);obj.c(:)];
      else
        assert(numel(obj.xi)==obj.nx,'Ilegal sizes of the x-vector')
      end
      out=obj.xi(:);
    end
    function out=get.p(obj)
      if obj.np==0; out=[]; return; end
      assert(~isempty(obj.xi),'Cannot determine the polynomial coefficients because the x-vector is undefined')
      out=obj.x(obj.p_idx,:);
    end
    function out=get.s(obj)
      if obj.ns==0; out=[]; return; end
      assert(~isempty(obj.xi),'Cannot determine the sine coefficients because the x-vector is undefined')
      out=obj.x(obj.s_idx,:);
    end
    function out=get.c(obj)
      if obj.ns==0; out=[]; return; end
      assert(~isempty(obj.xi),'Cannot determine the cosine coefficients because the x-vector is undefined')
      out=obj.x(obj.c_idx,:);
    end
    %% coefficients set functions

    %TODO: these functions break when there are multiple columns in y

    function obj=set.x(obj,x_in)
      if isempty(x_in); return; end
      assert(size(x_in,1)==obj.nx,['input x_in must have ',num2str(obj.nx),' rows, not ',num2str(size(x_in,1))])
      obj.xi=x_in;
    end
    function obj=set.p(obj,p_in)
      if isempty(p_in); return; end
      assert(numel(p_in)==obj.np,['Input p_in must have length ',num2str(obj.np),', not ',num2str(numel(p_in)),'.'])
      obj.xi(obj.p_idx)=p_in;
    end
    function obj=set.s(obj,s_in)
      if isempty(s_in); return; end
      assert(numel(s_in)==obj.ns,['Input s_in must have length ',num2str(obj.ns),', not ',num2str(numel(s_in)),'.'])
      obj.xi(obj.s_idx)=s_in;
    end
    function obj=set.c(obj,c_in)
      if isempty(c_in); return; end
      assert(numel(c_in)==obj.ns,['Input c_in must have length ',num2str(obj.ns),', not ',num2str(numel(c_in)),'.'])
      obj.xi(obj.c_idx)=c_in;
    end
    %% build the design matrix
    function out=get.A(obj)
      % get time domain for inversion
      t_now=obj.t-obj.t0;
      % init design matrix
      out=zeros(obj.ny,obj.nx);
      % build design matrix: polynomial coefficients
      for i=obj.p_idx
        out(:,i)=t_now.^(i-1);
      end
      % build design matrix: sinusoidal coefficients
      for i=obj.s_idx
        out(:,i)=sin(2*pi/obj.T(i-obj.np)*t_now);
      end
      % build design matrix: co-sinusoidal coefficients
      for i=obj.c_idx
        out(:,i)=cos(2*pi/obj.T(i-obj.np-obj.ns)*t_now);
      end
    end
    function out=get.AtQy(obj)
      out=transpose(obj.A)/obj.Qy;
    end
    function out=get.N(obj)
      out=obj.AtQy*obj.A;
    end
    function out=get.h(obj)
      out=obj.AtQy*obj.y;
    end
    function out=get.Qx(obj)
      out=inv(obj.N);
    end
    %notice, this is a ny-by-ny matrix
    function out=get.Pa(obj)
      out=obj.A/obj.N*obj.AtQy;
    end
    function out=get.e(obj)
      out=(eye(obj.ny)-obj.Pa)*obj.y;
    end
    %% forward modelling
    function out=get.y_sum(obj)
      out=obj.A*obj.x;
    end
    function out=get.y_comp(obj)
      out=zeros(obj.ny,obj.nx);
      for j=1:obj.nx
        out(:,j)=obj.A(:,j)*obj.x(j,:);
      end
    end
    function out=get.yp(obj); out=obj.y_comp(:,obj.p_idx); end
    function out=get.ys(obj); out=obj.y_comp(:,obj.s_idx); end
    function out=get.yc(obj); out=obj.y_comp(:,obj.c_idx); end
    function out=get.yr(obj); out=obj.y-obj.y_sum; end
    function out=get.rn(obj); out=norm(obj.yr); end
    function out=get.rrn(obj);out=obj.rn./norm(obj.y); end
    function obj=resample(obj,nt)
      if obj.ny==nt
        return
      end
      t_new=linspace(obj.t(1),obj.t(end),round(nt));
      y_new=interp1(obj.t,obj.y,t_new);
      obj=pardecomp(t_new,y_new,'T',obj.T,'np',obj.np,'t0',obj.t0,'x',obj.x);
    end
    %% inverse modelling
    function out=issolved(obj)
      out=~isempty(obj.xi);
    end
    function obj=lsq(obj)
      if ~obj.issolved
        %solve the system of linear equations
        obj.x=obj.N\obj.h;
      end
    end
    %% general plotting
    function plot(obj,varargin)
      % parse input arguments
      v=varargs.wrap('sources',{...
        {...
          'columns',        {1} , @iscell; ...
        },...
      },varargin{:});
      screen_position=200+[0,0,21,9]*50;
      for j=1:numel(v.columns)
        figure('Position',screen_position,'PaperPosition',screen_position)
        legend_str=cell(1,obj.nx+2);
        counter=0;
        plot(obj.t,obj.y(:,j),'b','LineWidth',2), hold on
        counter=counter+1;legend_str{counter}='original';
        plot(obj.t,obj.yr(:,j),'k','LineWidth',2)
        counter=counter+1;legend_str{counter}='residual';
        for i=1:numel(obj.p)
          plot(obj.t,obj.yp(:,i),'r','LineWidth',2)
          counter=counter+1;legend_str{counter}=['t^',num2str(i-1),':',num2str(obj.p(i))];
        end
        for i=1:numel(obj.s)
          plot(obj.t,obj.ys(:,i),'g','LineWidth',2)
          counter=counter+1;legend_str{counter}=['sin_',num2str(i),':',num2str(obj.s(i))];
        end
        for i=1:numel(obj.c)
          plot(obj.t,obj.yc(:,i),'m','LineWidth',2)
          counter=counter+1;legend_str{counter}=['cos_',num2str(i),':',num2str(obj.c(i))];
        end
        legend(legend_str,'location','eastoutside')
        title(['norm(x+res-y)=',num2str(norm(sum([obj.yp,obj.ys,obj.yc,obj.yr,-obj.y],2))),...
          newline,'T=',num2str(obj.T(:)')])
        fs=16;
        set(    gca,          'FontSize',fs);
        set(get(gca,'Title' ),'FontSize',round(fs*1.3));
        set(get(gca,'XLabel'),'FontSize',round(fs*1.1));
        set(get(gca,'YLabel'),'FontSize',round(fs*1.2));
        grid on
      end
    end
  end
end