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
    %% decomposition into pd_set
    function pd_set=split(obj,varargin)
      % NOTICE: the following input arguments are pretty much mandatory to define the 
      %         parametric regression (see below):
      % - 'T',[2*min(diff(t)),(t(end)-t(1))/2], @(i) isnumeric(i) || isempty(i);...
      % - 'np',0,                               @num.isscalar;...
      % NOTICE: the following input argumentsare good to define externally (so you can 
      %         interpret the estimated parameters)
      % - 'timescale',obj.x_units, @ischar;...
      % - 'epoch',      obj.epoch, @isdatetime;... 
      % NOTICE: this method expect obj to be simpledata-esque 
      % NOTICE: this is a wrapper for a obj.y with several columns (TODO: fix 
      %         pardecomp.init to handle multiple y columns) and re-assembled the estimated
      %         parameters in timeseries/gravity objects (grouped in the pd_set)
      v=varargs.wrap('sources',{....
        {...
          'timescale',obj.x_units, @ischar;...
          'epoch',      obj.epoch, @isdatetime;... %NOTICE: need simpletimeseries or derived class
          'quiet',          false, @islogical;...
          'parallel',       false, @islogical;...
        },...
      },varargin{:});
      %init loop vars
      t_pd=simpletimeseries.time2num(obj.t_masked,obj.epoch,v.timescale);
      y_pd=obj.y_masked;
      %convert t0 to numeric, honour v.timescale (obj.t2x only know about the obj.x_units timescale and obj.epoch)
      v=v.join({'t0',time.swap_units(obj.t2x(v.epoch),obj.x_units,v.timescale)});

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
            t_pd,y_pd(:,i),v.varargin{:}...
          ).lsq; %#ok<AGROW>
          if ~v.quiet; s=time.progress(s,i); end
          %sanity
          assert(rms(y_pd(:,i)-d(i).y_sum-d(i).yr)<1e-12,...
            ['Parameter decomposition failed for ',obj.labels{i},', RMS=',num2str(rms(y_pd(:,i)-d(i).y_sum-d(i).yr))])
        end
      end
      %init containers
      init=str2func(class(obj)); %use correct constructor
      pd_args=cell(1,4*pardecomp.xlength(d(1).np,d(1).T)); 
      coeffnames=pardecomp.xnames(d(1).np,d(1).T);
      common_args=[obj.varargin,{'silent',true,'x_units',v.timescale,'epoch',obj.epoch}];
      c=0;
      clearvars s; s.msg=['ts constituents pardecomp   ',obj.descriptor]; s.n=numel(coeffnames);
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
        o=init(time.zero_date,transpose(num.struct_deal(d,coeffnames{j}(1),i,[])),common_args{:});
        o.descriptor=[coeffnames{j},' of ',str.clean(obj.descriptor,'file')];
        o.units(:)={units}; o.labels(:)={labels};
        %append to pd_args
        pd_args{c+1}=coeffnames{j};
        pd_args{c+2}=o;
        c=c+2;
        %save timeseries represented by each coefficient
        clearvars o
        o=init(obj.t_masked,num.struct_deal(d,['y',coeffnames{j}(1)],[],i),common_args{:});
        o.descriptor=['p',num2str(i-1),' of ',str.clean(obj.descriptor,'file')];
        %restore gaps
        o=o.t_merge(obj.t);
        %append to pd_args
        pd_args{c+1}=['ts_',coeffnames{j}];
        pd_args{c+2}=o;
        c=c+2;
        if ~v.quiet; s=time.progress(s,j); end
      end
      %some sanity
      pdm.metadata=varargs(pd_args);
      pdm.fn=fieldnames(pdm.metadata);
      pdm.fg=pdm.fn(cells.strfind(pdm.fn,'ts_'));
      pdm.sum=num.struct_deal(d,'yr',[],1);
      for i=1:numel(pdm.fg)
        pdm.sum=pdm.sum+pdm.metadata.(pdm.fg{i}).y_masked;
      end
      pdm.check=pdm.sum-y_pd;
      if rms(pdm.check(:))>1e-12
        plotting.figure
        semilogy(rms(pdm.check))
        title(['RMS=',num2str(rms(pdm.check(:)))])
        error('Buildind time series constituents failed.')
      end
      %save everything into pd_set record (including metadata, done internally)
      pd_set=pardecomp.assemble(...
        'T',d(1).T,...,
        'np',d(1).np,...
        't0',d(1).t0,...
        't',d(1).t,...   %this is the abcissae used in the regression (actually: t-t0)
        'time',obj.t,... %this is the abcissae defined in obj, including gaps
        'epoch',obj.epoch,...
        'descriptor',obj.descriptor,...
        'timescale',v.timescale,...
        'quiet',v.quiet,...
        pd_args{:});
      %this is the abcissae defined in obj, excluding gaps
      pd_set.t_masked=obj.t_masked;
      clearvars s; s.msg=['statistics of  pardecomp of ',obj.descriptor]; s.n=3;
      %save residuals
      clearvars o
      o=init(obj.t_masked,num.struct_deal(d,'yr',[],1),common_args{:});
      o.descriptor=['residual of ',str.clean(obj.descriptor,'file')];
      %restore gaps
      o=o.t_merge(obj.t);
      pd_set.res=o;
      if ~v.quiet; s=time.progress(s,1); end
      %save norms
      clearvars o
      o=init(time.zero_date,num.struct_deal(d,'rn',[],1),common_args{:});
      o.descriptor=['norm of the residuals of ',str.clean(obj.descriptor,'file')];
      pd_set.norm=o;
      if ~v.quiet; s=time.progress(s,2); end
      %save norm ratio
      clearvars o
      o=init(time.zero_date,num.struct_deal(d,'rrn',[],1),common_args{:});
      o.descriptor=['signal and residual norms ratio for ',str.clean(obj.descriptor,'file')];
      pd_set.rnorm=o;
      if ~v.quiet; s=time.progress(s,3); end
    end
    %% reconstruction from pd_set
    function obj=join(pd_set,varargin)
      %sanity
      assert(pardecomp.ispd_set(pd_set),'Input ''pd_set'' must validate the pardecomp.ispd_set method')
      %easiner names
      coeffnameall=pardecomp.xnames(pd_set.np,pd_set.T);
      %parse inputs
      v=varargs.wrap('sources',{....
        {...
          'time',      pd_set.time,  @isdatetime;...
          'coeffnames',coeffnameall, @iscellstr;...
          'add_res',   false,        @islogical;...
        },...
      },varargin{:});
      %init output and copy metadata
      pd_set.metadata=varargs(pd_set.metadata).rm_empty;
      pd_set.varargin=pd_set.metadata.varargin;
      switch func2str(pd_set.init)
        case 'gravity';obj=eval([func2str(pd_set.init),'.unit(gravity.width2lmax(pd_set.width),''t'',v.time,pd_set.varargin{:},''scale'',0);']);
        otherwise;     obj=eval([func2str(pd_set.init),'.zero(v.time,pd_set.width,pd_set.varargin{:});']);
      end
      %NOTICE: do not match epochs explicitly, it should be done at initialization (above), otherwise you get artificial shifts
      %obj.epoch=pd_set.epoch;
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
        if isfield(pd_set,tname) && pd_set.(tname).istxequal(v.time)
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
        s.msg=['Parametric reconstruction of ',obj.descriptor,' with components ',strjoin(coeffnameall(coeffidx==1),', ')]; s.n=obj.width;
        for i=1:obj.width
          obj=obj.set_cols(i,pardecomp(...
            simpletimeseries.time2num(obj.t,pd_set.epoch,pd_set.timescale),[],...
            'T' ,pd_set.T,...
            'np',pd_set.np,...
            't0',pd_set.t0,...
            'x',coeffval(:,i)...
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
        'timescale',      simpletimeseries.parameters('x_units'), @ischar;...
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
          'time',simpletimeseries.num2time(v.t,v.epoch,v.timescale), @isdatetime;...
          'quiet',        false, @islogical;...
        },...
      },varargin{:});
      %initialize output
      pd_set=struct('t',v.t,'T',v.T,'np',v.np,'t0',v.t0,'time',v.time);
      %init looping
      records=struct([]);
      rnames={'width','class','metadata'};
      coeffnames=pardecomp.xnames(pd_set.np,pd_set.T);
      s.msg=['Assemble pd_set  pardecomp  ',v.descriptor]; s.n=numel(coeffnames);
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
              for m={'descriptor','units','labels'}
                rvalue=rmfield(rvalue,m{1});
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
        if ~v.quiet; s=time.progress(s,i); end
      end
      %NOTICE: this is important to join the pd_set correctly
      pd_set.metadata=structs.rm_empty(records.metadata);
      %NOTICE: this is largely for information only (at least I think so)
      pd_set.init=str2func(records.class);
      pd_set.width=records.width;
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
          'time_unit', pd_set.timescale, @ischar;...
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
            period=num2str(time.duration2num(simpletimeseries.timescale(...
              pd_set.T(ps),pd_set.timescale...
            ),v.time_unit));
          case 'cosine'
            pc=pc+1;
            period=num2str(time.duration2num(simpletimeseries.timescale(...
              pd_set.T(pc),pd_set.timescale...
            ),v.time_unit));
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
    function out=test(varargin)
      v=varargs.wrap('sources',{...
        {...
          'mode',   'all', @ischar;...
          'plot',   false, @islogical;...
          'print',   true, @islogical;...
        }},...
      varargin{:});
      switch v.mode
      case 'all'
        mode_list={'for','back','split','join'};
        for i=1:numel(mode_list)
          pardecomp.test('mode',mode_list{i},'plot',true,'print',true);
        end
      case 'parameters'
        %test parameters
        out.step=1;
        out.n=10000;
        out.poly_coeffs=[1 3 5]./[1 out.n out.n^2];
        out.sin_periods=out.n/out.step./[2 5];
        out.sin_coeffs=[0.5 3];
        out.cos_coeffs=[2 0.8];
        out.sin_periods_assumed=out.sin_periods;
        %derived parameters
        out.t=transpose(1:out.step:(out.n*out.step));
%TODO: implement multiple columns in Y
  %       obj.randn_scale=[0.1,1,10];
        out.randn_scale=0.1;
      case {'for','forward','forwards'}
        %retrieve parameterss
        p=pardecomp.test('mode','parameters');
        %forward modelling
        out=pardecomp(p.t,[],...
          'np',numel(p.poly_coeffs),...
          'T',p.sin_periods,...
          'p',p.poly_coeffs,...
          's',p.sin_coeffs,...
          'c',p.cos_coeffs...
        ).forward;
        %show results
        if v.plot
          for i=1:numel(p.randn_scale)
            out.plot(varargin{:},'columns',{i});
          end
          plotting.enforce;
        end
      case {'back','backward','backwards'}
        %retrieve parameters
        p=pardecomp.test('mode','parameters');
        %retrieve forward model
        ref=pardecomp.test('mode','forward');
        %add noise
        y=ref.y_sum;
        y=y*ones(size(p.randn_scale))+randn(size(y))*p.randn_scale;
        %backwards modelling
        out=pardecomp(p.t,y,...
          'np',numel(p.poly_coeffs),...
          'T',p.sin_periods_assumed...
        ).lsq;
        %show results
        if v.plot
          for i=1:numel(p.randn_scale)
            out.plot(varargin{:},'columns',{i});
          end
          plotting.enforce;
        end
        %inform
        if v.print
          out.print([],ref);
        end
      case {'split','split-unstable'}
        t_start=datetime('now');
        out_timescale='hours';
        tst_timescale='seconds';
        tsx_timescale='days';
        switch v.mode
          case 'split-unstable'
            epoch_pd=t_start-years(1);
            epoch_ts=t_start-days(0.5);
          case 'split'
            delta=years(randn(1)*10);
            epoch_pd=t_start-delta;
            epoch_ts=t_start-delta;
        end
        %get timeseries objects
        pd=pardecomp.test('mode','forwards');
        ts=pd.ts(tst_timescale,epoch_ts,'x_units',tsx_timescale);
        %split it
        out=pardecomp.split(ts,...
          'np',pd.np,...
          'T',time.swap_units(pd.T,tst_timescale,out_timescale),...
          'timescale',out_timescale,...
          'epoch',epoch_pd...
        );
        %show results
        if v.plot
          plotting.figure;
          fn=fieldnames(out);
          fn=fn(contains(fn,'ts_'));
          legend_str=cell(size(fn));
          for i=1:numel(fn)
            out.(fn{i}).plot
            legend_str{i}=strrep(fn{i},'ts_','');
          end
          ts.plot
          legend_str{i+1}='ts';
          plot(ts.t,pd.y_sum)
          legend_str{i+2}='pd.y_sum';
          plotting.enforce(...
            'plot_title',['RMS(pd.y_sum-ts.y)=',num2str(rms(pd.y_sum-ts.y))],...
            'plot_legend',legend_str...
          );
        end
        %inform
        if v.print
          ts.print
          disp(pardecomp.table(out));
          disp(['RMS(pd.y_sum-ts.y)=',num2str(rms(pd.y_sum-ts.y))])
        end
      case {'join','join-unstable'}
        pd_set=pardecomp.test('mode',strrep(v.mode,'join','split'));
        pd=pardecomp.test('mode','forwards');
        out=pardecomp.join(pd_set);
        %show results
        if v.plot
          plotting.figure;
          subplot(2,1,1)
          out.plot
          plot(out.t,pd.y_sum)
          plotting.enforce(...
            'plot_title',['RMS(split/joined-forwards)=',num2str(rms(out.y-pd.y_sum))],...
            'plot_legend',{'split/joined','forwards'}...
          );
          subplot(2,1,2)
          plot(out.t,out.y-pd.y_sum)
          plotting.enforce(...
            'plot_title','split/joined-forwards',...
            'plot_legend_location','none'...
          );
        end
        %inform
        if v.print
          out.print
          disp(pardecomp.table(pd_set));
          disp(['RMS(split/joined-forwards)=',num2str(rms(out.y-pd.y_sum))])
        end
      end
      if nargout==0
        clearvars obj
      end
    end
  end
  methods
    function obj=pardecomp(t,y,varargin)
      %NOTICE: t0 is zero-value of the (numeric) time domain that is used to define the
      %        regression parameters; therefore, it is arbitrary.
      %NOTICE: t0 has nothing to do with the epoch from simpletimeseries and its children 
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
    function out=varargin(obj)
      out={...
        'T' ,obj.T;...
        'np',obj.np;...
        't0',obj.t0;...
        'p' ,obj.p;...
        's' ,obj.s;...
        'c' ,obj.c;...
        'x' ,obj.x;...
        't' ,obj.t;...
        'y' ,obj.y;...
      };
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
        disp([str.tabbed(field,tab),' : ',str.show(transpose(value(:)),'fmt',fmt)])
      else
        disp([str.tabbed(label,tab),' : ',str.show(transpose(value(:)),'fmt',fmt)])
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
               obj.disp_field('t0'     ,tab,[]               ,''     ,fmtf);
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
    function obj=forward(obj)
      obj.y=obj.y_sum;
    end
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
        warning off
        obj.x=obj.N\obj.h;
        warning on
      end
    end
    %% convert to time series
    function out=ts(obj,t_units,epoch,varargin)
      w=size(obj.y,2);
      out=simpletimeseries(...
        epoch+time.num2duration(obj.t,t_units),...
        obj.y,...
        'labels',strcat(cellstr(repmat('label-',w,1)),cellstr(num2str((1:w)'))),...
        'units', strcat(cellstr(repmat('unit-', w,1)),cellstr(num2str((1:w)'))),...
        'descriptor','pardecomp obj',...
        varargin{:},...
        'format','datetime'...
      );
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
        legend_str=cell(1,obj.nx);
        counter=0; hold on
        if ~isempty(obj.y)
          plot(obj.t,obj.y(:,j),'b','LineWidth',2), hold on
          counter=counter+1;legend_str{counter}='original';
          plot(obj.t,obj.yr(:,j),'k','LineWidth',2)
          counter=counter+1;legend_str{counter}='residual';
        end
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
        if ~isempty(obj.y)
          title(['norm(x+res-y)=',num2str(norm(sum([obj.yp,obj.ys,obj.yc,obj.yr,-obj.y],2))),...
            newline,'T=',num2str(obj.T(:)')])
        else
          title(['T=',num2str(obj.T(:)')])
        end
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