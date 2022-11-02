classdef simpletimeseries < simpledata
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list={...
      'format',    'modifiedjuliandate',@(i)ischar(i)||isdatetime(i);...
      't_tol',     2e-6,                @num.isscalar;...
      'timesystem','utc',               @ischar;...
      'data_dir'   file.orbdir('data'), @ischar;...
      'x_units',   'seconds',           @ischar;...
      'epoch',     time.zero_date,      @isdatetime;...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'timesystem','epoch'};
    %define periods when CSR calmod is upside down
    csr_acc_mod_invert_periods=datetime({...
      '2016-01-28','2016-03-02';...
    });
  end
  properties(Constant)
    valid_timesystems={'utc','gps'};
  end
  %NOTE: edit this if you add a new parameter (if read only)
  properties(SetAccess=private)
    step
  end
 %These parameters should not modify the data in any way; they should
  %only describe the data or the input/output format of it.
  %NOTE: edit this if you add a new parameter (if read/write)
  properties(GetAccess=public,SetAccess=public)
    format
    t_tol
    timesystem
    data_dir
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    epochi %absolute epoch (datetime class), from which x in simpledata is relative to
  end
  %calculated only when asked for
  properties(Dependent)
    t
    t_formatted   %this handles the numeric/char version of t
    epoch
    start
    stop
    tsys
    first
    last
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(simpletimeseries.parameter_list); end
      out=v.picker(varargin{:});
    end
    function out=timescale(in,time_units)
      if ~exist('time_units','var') || isempty(time_units)
        time_units=simpletimeseries.parameters('x_units');
      end
      %NOTICE: this method handles both duration and numeric inputs, returning the other data type.
      out=time.num2duration(in,time_units);
    end
    function out=valid_t(in)
      out=isdatetime(in);
    end
    function out=valid_epoch(in)
      out=isdatetime(in) && isscalar(in);
    end
    function out=valid_timesystem(in)
      switch lower(in)
      case simpletimeseries.valid_timesystems
        out=true;
      otherwise
        out=false;
      end
    end
    function out=time2num(in,epoch,time_units)
      %NOTICE: when calling obj.x_units=..., there can be the case that 'in' is empty
      if isempty(in)
        out=[];
        return
      end
      if ~exist('epoch','var') || isempty(epoch)
        epoch=in(1);
      end
      if ~exist('time_units','var') || isempty(time_units)
        time_units=simpletimeseries.parameters('x_units');
      end
      if isfinite(epoch)
        out=simpletimeseries.timescale(in-epoch,time_units);
      elseif any(isfinite(in))
        out=simpletimeseries.timescale(in-min(in),time_units);
      else
        out=Inf(size(in));
      end
    end
    function out=num2time(in,epoch,time_units)
      if ~exist('epoch','var') || isempty(epoch)
        error('need input ''epoch''.')
      end
      if ~exist('time_units','var') || isempty(time_units)
        time_units=simpletimeseries.parameters('x_units');
      end
      out=epoch+simpletimeseries.timescale(in,time_units);
    end
    function out=ist(mode,t1,t2,tol)
      %expect vectors as well
      if numel(t1)~=numel(t2)
        out=false;
        return
      end
			%this handles infinites
      if t1(:)~=t2(:)
        out=simpledata.isx(mode,seconds(t1(:)-t2(:)),0,tol);
      else
        out=true;
      end
    end
    function presence=ispresent(parser,fields)
      % defaults
      if ~exist('fields','var') || isempty(fields)
        fields={'t','x'};
        check_for_concurrence=true;
      else
        check_for_concurrence=false;
      end
      %sanity
      if ~iscell(fields)
        error('input argument ''fields'' must be a cell array.')
      end
      % look for existence
      for i=1:numel(fields)
        if any(strcmp(parser.Parameters,fields{i}))
          presence.(fields{i})=~any(strcmp(parser.UsingDefaults,fields{i}));
        else
          presence.(fields{i})=isfield(parser.Unmatched,fields{i});
        end
      end
      %this is often how this routine is called
      if check_for_concurrence
        %cannot have both 't' and 'x'
        if presence.x && presence.t
          error('cannot handle both inputs ''x'' and ''t''.')
        end
      end
    end
    function out=transmute(in)
      if isa(in,'simpletimeseries')
        %trivial call
        out=in;
      else
        %transmute into this object
        if obj.is_timeseries
          out=simpletimeseries(in.t,in.y,in.varargin{:});
        else
          out=simpletimeseries(in.x,in.y,in.varargin{:});
        end
      end
    end
    function out=timestep(in,varargin)
      p=machinery.inputParser;
      p.addRequired( 'in',                @isdatetime);
      p.addParameter('nsigma',    4,      @num.isscalar);
      p.addParameter('max_iter',  10,     @num.isscalar);
      p.addParameter('sigma_iter',2,      @num.isscalar);
      p.addParameter('sigma_crit',1e-9,   @num.isscalar);
      p.addParameter('max_mean_ratio',1e3,@num.isscalar);
      p.addParameter('curr_iter', 0,      @num.isscalar);
      p.addParameter('disp_flag', false,  @islogical);
      p.addParameter('time_units',simpletimeseries.parameters('x_units'), @ischar);
      % parse it
      p.parse(in,varargin{:});
      %handle singularities
      switch numel(in)
        case 0
          error('cannot handle empty time stamps')
        case 1
          out=0;
          return
      end
      %get numeric diff of time
      tdiff=simpletimeseries.timescale(diff(in),p.Results.time_units);
      %large jumps produce erroneous results, so get rid of those first
      while std(tdiff)~=0 && max(tdiff)/mean(tdiff)>p.Results.max_mean_ratio
        %save stats
        stdtdiff=std(tdiff);
        ratiotdiff=max(tdiff)/mean(tdiff);
        %remove large gaps
        tdiff=simpledata.rm_outliers(tdiff,varargin{:});
        %send feedback
        if p.Results.disp_flag
          disp([' removed ',num2str(sum(isnan(tdiff))),' large gaps, since ',...
            'std(delta t) is ',num2str(stdtdiff),' and ',...
            'max(delta t) is ',num2str(ratiotdiff),' times larger than mean(delta).'])
        end
        %remove nans
        tdiff=tdiff(~isnan(tdiff));
      end
      %get diff of time domain without jumps
      outdiff=simpledata.rm_outliers(tdiff,varargin{:});
      %get rid of nans
      outdiff=outdiff(~isnan(outdiff));
      %check if there are still lots of gaps in the data
      if std(outdiff)>p.Results.sigma_crit*mean(outdiff) && p.Results.curr_iter < p.Results.max_iter
        %reduce sigma
        nsigma_new=p.Results.nsigma/p.Results.sigma_iter;
        %send feedback
        if p.Results.disp_flag
          disp([' failed to determine the timestep, since std(delta t) is ',num2str(std(outdiff)),...
            '. Reducing NSIGMA from ',num2str(p.Results.nsigma),' to ',num2str(nsigma_new),'.'])
        end
        %recursive call
        vararginnow=cells.vararginclean(varargin,{'nsigma','curr_iter','disp_flag'});
        out=simpletimeseries.timestep(in,...
          'nsigma',nsigma_new,...
          'curr_iter',p.Results.curr_iter+1,...
          'disp_flag',false,...
          vararginnow{:});
      elseif isempty(outdiff)
        %dead end, sigma was reduced too much and all data is flagged as
        %outliers: nothing to do but to give some estimated of the previous
        %sigma (rounded to micro-seconds to avoid round off errors)
        vararginnow=cells.vararginclean(varargin,{'nsigma'});
        outdiff=simpledata.rm_outliers(tdiff,...
          'nsigma',p.Results.nsigma*p.Results.sigma_iter,...
          vararginnow{:});
        out=simpletimeseries.timescale(...
          round(...
            mean(...
              outdiff(~isnan(outdiff))...
            )*1e6...
          )*1e-6...
        );
      else
        out=simpletimeseries.timescale(outdiff(1));
      end
      %send feedback if needed
      if p.Results.disp_flag
        disp([' final timestep is ',char(out),'.'])
      end
    end
    function v=fix_interp_over_gaps_narrower_than(v)
      if ~iscell(v)
        error(['expecting input ''v'' to be a cell array, not a ',class(v),'.'])
      end
      for i=1:numel(v)
        if strcmp(v{i},'interp_over_gaps_narrower_than')
          if isduration(v{i+1})
            v{i+1}=simpletimeseries.timescale(v{i+1});
          end
          break
        end
      end
    end
    %% internally-consistent naming of satellites
    function out=translatesat(in)
      %search for satellite name
      switch lower(in)
        case {'champ','ch'}
          out='ch';
        case {'swarm-a','swarma','swarm a','swma','sa','l47'}
          out='sa';
        case {'swarm-b','swarmb','swarm b','swmb','sb','l48'}
          out='sb';
        case {'swarm-c','swarmc','swarm c','swmc','sc','l49'}
          out='sc';
        case {'goce','go'}
          out='go';
        case {'unknown','test'}
          out=in;
        otherwise
          success=false;
          %try grace.m
          [out,success]=machinery.trycatch(success,'grace:BadSat',@grace.translatesat,{in});
          %add additional class calls to translatesat (using the same structure above)

          %check something worked
          assert(success,['Cannot translate satellite ''',in,'''.'])
      end
    end
    function out=translatesatname(in)
      %search for satellite name
      switch simpletimeseries.translatesat(in)
        case 'ch'; out='CHAMP';
        case 'sa'; out='Swarm-A';
        case 'sb'; out='Swarm-B';
        case 'sc'; out='Swarm-C';
        case 'go'; out='GOCE';
        case {'unknown','test'}; out=in;
        otherwise
          success=false;
          %try grace.m
          [out,success]=machinery.trycatch(success,'grace:BadSatName',@grace.translatesatname,{in});
          %add additional class calls to translatesat (using the same structure above)

          %check something worked
          assert(success,['Cannot translate satellite ''',in,'''.'])
      end
    end
    %% consistent reference frame names
    function out=translateframe(in)
      if ischar(in)
        switch lower(in)
        case {'crs','crf','eci','icrf','gcrf','j2000','eme2000','celestial','inertial'}
          out='crf';
        case {'m50'}
          out='m50';
        case {'teme'}
          out='teme';
        case {'trs','trf','ecf','ecef','itrf','terrestrial','rotating','co-rotating',}
          out='trf';
        case {'body','satellite','srf'}
          out='srf';
        otherwise
          out='';
        end
      else
        out='';
      end
    end
    function out=isframe(in)
      out=~isempty(simpletimeseries.translateframe(in));
    end
    %% import methods
    %NOTICE: the mat-file handling in this method is so that there are mat files duplicating the raw data: one raw file, one mat file. This is not a datastorage-type of structuring the data.
    %TODO: consider retiring cut24hrs
    function obj=batch_import(filename,varargin)
      %additional options, from dependencies:
      % - simpletimeseries.import:
      % p.addParameter('format','',@ischar);
      % - file.load_mat:
      % p.addParameter('data_var' ,'out', @ischar);
      % - file.save_mat:
      % p.addParameter('save_mat', true, @(i) isscalar(i) && islogical(i))
      % p.addParameter('data_var','out', @ischar);
      % - file.delete_compressed:
      % p.addParameter('del_arch', true, @(i) isscalar(i) && islogical(i))
      p=machinery.inputParser;
      p.addRequired( 'filename',       @(i) ischar(i) || iscellstr(i));
      p.addParameter('cut24hrs', true, @(i) isscalar(i) && islogical(i))
      p.parse(filename,varargin{:})
      %unwrap wildcards and place holders (output is always a cellstr)
      filename=file.unwrap(filename,varargin{:});
      %loop over all files (maybe only one, resolved by cells.scalar)
      for i=1:numel(filename)
        disp([' reading data from file ',filename{i}])
        %read the data from a single file
        obj_now=simpletimeseries.import(filename{i},varargin{:});
        %skip if empty
        if isempty(obj_now)
          continue
        end
        %handle cutting data to requested periods
        if p.Results.cut24hrs
          %determine current day
          day_now=datetime(yyyymmdd(obj_now.t(round(obj_now.length/2))),'ConvertFrom','yyyymmdd');
          %get rid of overlaps
          obj_now=obj_now.trim(day_now,day_now+hours(24)-obj_now.step);
        end
        %append or initialize
        if ~exist('obj','var')
          obj=obj_now;
        else
          try
            obj=obj.append(obj_now);
          catch
            obj=obj.augment(obj_now);
          end
        end
      end
      %in case there are no files, 'filename' will be empty and the loop will be skipped
      if ~exist('obj','var')
        obj=[];
      end
    end
    function obj=import(filename,varargin)
      %additional options, from dependencies
      % - file.load_mat:
      % p.addParameter('data_var' ,'out', @ischar);
      % - file.save_mat:
      % p.addParameter('save_mat', true, @(i) isscalar(i) && islogical(i))
      % p.addParameter('data_var','out', @ischar);
      % - file.delete_compressed:
      % p.addParameter('del_arch', true, @(i) isscalar(i) && islogical(i))
      p=machinery.inputParser;
      p.addParameter( 'format','',@ischar);
      p.parse(varargin{:})
      %resolve scalar cell
      filename=cells.scalar(filename,'get');
      %make sure this is a filename
      assert(ischar(filename),['Input ''filename'' must be of class ''char'', not ''',...
        class(filename),'''; consider using the simpletimeseries.batch_import method.'])
      %try to load mat file
      [obj,loaded_flag]=file.load_mat(filename,...
        'default_dir',simpletimeseries.parameters('value','data_dir'),...
        'data_var','obj',...
        varargin{:}...
      );
      %if loaded something, we're done
      if loaded_flag, return, end
      %enforce format given as argument
      if ~isempty(p.Results.format)
        format=p.Results.format;
      else
        %split into parts and propagate the extension as the format
        [~,~,format]=fileparts(filename);
      end
      %call format interface aggregator
      obj=simpletimeseries.import_format(filename,format);
      %possibly save this object as a mat file
      file.save_mat(obj,filename,'data_var','obj',varargin{:});
      %possibly delete uncompressed file
      file.delete_uncompressed(filename,varargin{:});
    end
    function obj=import_format(filename,format)
      %branch on extension/format ID
      switch format
      case 'seconds'
        %define data cols
        data_cols=2:4;
        %load the data
        raw=file.textscan(filename,'%f %f %f %f');
        %get the labels (from the header)
        labels={'RL06_grav','RL05_grav','res_grav'};
        %build units
        units=cell(size(data_cols));
        units(:)={'m/s'};
        t=datetime([2000 01 01 0 0 0])+seconds(raw(:,1));
        %building object
        obj=simpletimeseries(t,raw(:,data_cols),...
          'format','datetime',...
          'units',units,...
          'labels', labels,...
          'timesystem','gps',...
          'descriptor',['KBR postfits ',filename]...
        );
      case {'TAS-2','TAS-12'}
        %define data cols
        data_cols=2:1+str2double(strrep(format,'TAS-',''));
        %load the data
        raw=file.textscan(filename,'%f %f %f %f');
        switch numel(data_cols)
        case 2
          labels={'SS_DistanceMeasurementError','NonGravDiffAccelMeasErr'};
           units={'m','m/s^2'};
        case 12
          labels={...
            'TestMass2_LinAcc_Noise_X','TestMass2_LinAcc_Noise_Y','TestMass2_LinAcc_Noise_X',...
            'TestMass3_LinAcc_Noise_X','TestMass3_LinAcc_Noise_Y','TestMass3_LinAcc_Noise_X',...
            'TestMass5_LinAcc_Noise_X','TestMass5_LinAcc_Noise_Y','TestMass5_LinAcc_Noise_X',...
            'TestMass6_LinAcc_Noise_X','TestMass6_LinAcc_Noise_Y','TestMass6_LinAcc_Noise_X',...
            };
          units=cell(size(labels));
          units(:)={'m/s^2'};
        end
        t=datetime([2000 01 01 0 0 0])+seconds(raw(:,1));
        %building object
        obj=simpletimeseries(t,raw(:,data_cols),...
          'format','datetime',...
          'units',units,...
          'labels', labels,...
          'timesystem','gps',...
          'descriptor',[format,' ',filename]...
        );
      otherwise
        success=false;
        %grace.m
        [obj,success]=machinery.trycatch(success,'grace:BadImportFormat',@grace.import_format,{filename,format});
        %add additional class calls to translatesat (using the same structure above)

        %check something worked
        assert(success,['Cannot handle files of type ''',format,'''.'])
      end
    end
    function obj=GRACEaltitude(varargin)
      p=machinery.inputParser;
      p.addParameter('datafile',file.resolve_home(fullfile('~','data','grace','altitude','GRACE.altitude.dat')));
      p.parse(varargin{:});
      obj=simpletimeseries.import(p.Results.datafile,...
        'format','mjd',...
        'cut24hrs',false...
      );
		end
    %% utilities
    function out=list(start,stop,period)
      p=machinery.inputParser;
      p.addRequired( 'start',   @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'stop',    @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'period',  @(i) isscalar(i) && isduration(i));
      p.parse(start,stop,period)
      out=datetime([],[],[]);
      for i=1:ceil((stop-start)/period)+1
        out(i)=start+(i-1)*period;
      end
      %trim end if after stop date
      if out(end)>stop
        out=out(1:end-1);
      end
    end
    function out=t_mergev(obj_list)
      for i=2:numel(obj_list)
        obj_list{1}=obj_list{1}.t_merge(obj_list{i}.t);
      end
      out=obj_list{1}.t;
    end
    %% constructors
    function out=one(t,width,varargin)
      out=simpletimeseries(t(:),ones(numel(t),width),'descriptor','unit',varargin{:});
    end
    function out=zero(t,width,varargin)
      out=simpletimeseries(t(:),zeros(numel(t),width),'descriptor','zero',varargin{:});
    end
    function out=randn(t,width,varargin)
      out=simpletimeseries(t(:),randn(numel(t),width),'descriptor','randn',varargin{:});
    end
    function out=sinusoidal(t,w,varargin)
      y=cell2mat(arrayfun(@(i) sin((t(:)-t(1))*pi/i),w,'UniformOutput',false));
      out=simpletimeseries(t(:),y,'descriptor','sinusoidal',varargin{:});
    end
    %% general test for the current object
    function test(method,l,w)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      if ~exist('l','var') || isempty(l)
        l=1000;
      end
      if ~exist('w','var') || isempty(w)
        w=3;
      end

      %get common parameters
      args=varargs(simpledata.test_parameters('args',l,w));
      %need to replace x_units
      args.x_units='days';
      now=juliandate(datetime('now'),'modifiedjuliandate');
      t=datetime(now,           'convertfrom','modifiedjuliandate'):...
        datetime(now+round(l)-1,'convertfrom','modifiedjuliandate');
      %init object
      a=simpletimeseries(...
          t,...
          simpledata.test_parameters('y_all',l,w),...
          'mask',simpledata.test_parameters('mask',l,w),...
          args.varargin{:}...
        );

      switch method
        case 'all'
          %TODO: add 'component_split'
          for i={'calibrate_poly','median','fill-resample','append','trim','resample','extend','slice','pardecomp'}
            simpletimeseries.test(i{1},l);
          end
        case 'component_split'
          a=simpletimeseries.randn(t,w,args.varargin{:});
          a=a.scale(0.05);
          idx=2:4;
          for i=1:numel(idx)
            as=simpletimeseries.sinusoidal(t,days(l/3)/idx(i)*ones(1,w),args.varargin{:});
            as=as.scale(rand(1,w));
            a=a+as;
          end
          a.descriptor='original';
          i=0;c=1;
          i=i+1;h{i}=plotting.figure('visible','on');
          subplot(2,2,1)
          a.plot('column',c)
          legend off
          subplot(2,2,2)
          [out,res,J,segs]=a.component_ampl(days(2*l/5));
          for si=1:numel(segs)
            segs{si}.plot('column',c)
          end
          title('segmented'); legend off
          subplot(2,2,3)
          plot(out(:,1))
          title('average and repeated'); legend off
          subplot(2,2,4)
          plot(res(:,1))
          title(['residuals (relative norm=',num2str(J),')']); legend off
          %TODO: finish this test
        case 'calibrate_poly'
          a=simpletimeseries.sinusoidal(t,days(l./(1:w)),args.varargin{:});
          bn=simpletimeseries.randn(t,w,args.varargin{:});
          a=a+bn.scale(0.05);
          b=a.scale(rand(1,w))+bn.scale(0.1)+ones(a.length,1)*randn(1,w);
          c=a.calibrate_poly(b);
          figure
          a.plot('columns',1)
          b.plot('columns',1)
          c.plot('columns',1)
          legend('uncal','target','cal')
          title(method)
        case 'median'
          figure
          a.plot(            'columns',1,'line',{'o-'});
          a.medfilt(10).plot('columns',1,'line',{'+-'});
          a.median( 10).plot('columns',1,'line',{'x-'});
          legend('origina','median','medfilt')
          title(method)
        case 'fill-resample'
          %TODO: make this test more illustrative
          b=a.resample;
          c=a.fill;
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          c.plot('columns',1,'line',{'x-'})
          legend('original','fill','resample')
          title(method)
        case 'append'
          b=a.append(...
            simpletimeseries(...
              a.stop+(round(l/3):round(4*l/3)-1),...
              simpledata.test_parameters('y_all',l,w),...
              'mask',simpledata.test_parameters('mask',l,w),...
              args.varargin{:}...
            )...
          );
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          legend('original','appended')
          title(method)
        case 'trim'
          b=a.trim(...
            datetime(now+round(l*0.3),'convertfrom','modifiedjuliandate'),...
            datetime(now+round(l*0.7),'convertfrom','modifiedjuliandate')...
          );
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          legend('original','trimmed')
          title(method)
        case 'resample'
          b=a.resample(...
            days(0.8) ...
          );
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          legend('original','resampled')
          title(method)
        case 'extend'
          b=a.extend(...
            round(l/4) ...
          ).extend(...
            -round(l/2) ...
          );
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          legend('original','extended')
          title(method)
        case 'slice'
          b=a.slice(...
            datetime(now+round( l/5),'convertfrom','modifiedjuliandate'),...
            datetime(now+round( l/2),'convertfrom','modifiedjuliandate')...
          );
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          legend('original','sliced')
          title(method)
        case 'pardecomp'
          a=simpletimeseries(...
              t,...
              simpledata.test_parameters('y_all_T',l,w),...
              'mask',simpledata.test_parameters('mask',l,w),...
              args.varargin{:}...
            );
          %test parameters
          poly_coeffs=simpledata.test_parameters('y_poly_scale');
          sin_periods=simpledata.test_parameters('T',l);
           sin_coeffs=simpledata.test_parameters('y_sin_scale');
           cos_coeffs=simpledata.test_parameters('y_cos_scale');
          %inform
          disp(['poly_coeffs : ',num2str(poly_coeffs(:)')])
          disp(['sin_periods : ',num2str(sin_periods(:)')])
          disp(['sin_coeffs  : ',num2str(sin_coeffs(:)')])
          disp(['cos_coeffs  : ',num2str(cos_coeffs(:)')])
          %derived parameters
          [b,pd_set]=a.pardecomp(...
            'np',numel(poly_coeffs),...
            'T',sin_periods,...
            'timescale','days'...
          );
          figure;
          a.plot('columns',1,'line',{'o-'})
          b.plot('columns',1,'line',{'+-'})
          legend('original','pardecomp')
          title(method)

          disp(pardecomp.table(pd_set,'tablify',true))
        case 'x_units'
          t=datetime('now')-years(10:-1:0);
          a=simpletimeseries.zero(t,2,'x_units','seconds');
          b=simpletimeseries.zero(t,2,'x_units','years');
          disp('First column was created with x_units=seconds and second column with x_units=years')
          disp('t-domains:')
          disp([a.t,b.t])
          disp('x-domains:')
          disp([a.x,b.x])
          assert(all(a.t==b.t),'test failed')
        case {'minus','minus-scalar'}
          switch method
            case 'minus',        t=datetime('now')-years(10:-1:0);
            case 'minus-scalar', t=datetime('now');
          end
          a=simpletimeseries.one(t,2,'x_units','seconds','epoch',datetime('now'));
          b=simpletimeseries.one(t,2,'x_units','years',  'epoch',datetime('now')+days(1));
          c=a-b;
          disp('First  column: x_units=seconds,epoch=now')
          disp('Second column: x_units=year   ,epoch=now+1 day')
          disp('Third  column: first - second objects')
          disp('epochs:')
          disp(datestr([a.epoch,b.epoch,c.epoch]))
          disp('x_units:')
          disp([a.x_units,' ',b.x_units,' ',c.x_units])
          disp('t-domains:')
          disp([a.t,b.t,c.t])
          disp('x-domains:')
          disp([a.x,b.x,c.x])
          disp('time2num:')
          disp([...
            simpletimeseries.time2num(a.t,a.epoch,a.x_units),...
            simpletimeseries.time2num(b.t,b.epoch,b.x_units),...
            simpletimeseries.time2num(c.t,c.epoch,c.x_units)...
          ])
          disp('num2time:')
          disp([...
            simpletimeseries.num2time(a.x,a.epoch,a.x_units),...
            simpletimeseries.num2time(b.x,b.epoch,b.x_units),...
            simpletimeseries.num2time(c.x,c.epoch,c.x_units)...
          ])
          assert(all(c.y(:)==0),'test failed')
      end
    end
  end
  methods
    %% constructor
    function obj=simpletimeseries(t,y,varargin)
      % input parsing
      p=machinery.inputParser;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y', @(i) simpledata.valid_y(i));
      %create argument object, declare and parse parameters, save them to obj
      [v,p]=varargs.wrap('parser',p,'sources',{simpletimeseries.parameters('obj')},'mandatory',{t,y},varargin{:});
      % get datetime
      [t,f]=time.ToDateTime(t,p.Results.format);
      %check if epoch and/or x_units are given
      if time.iszero(v.epoch); v.epoch=t(1); end
      %make sure the x_units are relevant to time
      v.x_units=time.translate_units(v.x_units);
      %call superclass (create empty object, assignment comes later)
      obj=obj@simpledata(simpletimeseries.time2num(t,v.epoch,v.x_units),y,...
        v.varargin{:}...
      );
      % save the arguments v into this object
      obj=v.save(obj,{'t','y'});
      %save input format (can be different from p.Results.format)
      obj.format=f;
    end
    function obj=assign(obj,y,varargin)
      p=machinery.inputParser;
      p.addRequired( 'y'         ,          @(i) simpledata.valid_y(i));
      p.addParameter('t'         ,obj.t,    @(i) simpletimeseries.valid_t(i));
      p.addParameter('epoch'     ,obj.epoch,@(i) simpletimeseries.valid_epoch(i));
      % parse it
      p.parse(y,varargin{:});
      % simpler names
      presence=simpletimeseries.ispresent(p);
      %if 't' is not present, then pass it on to simple data
      if ~presence.t
        obj=assign@simpledata(obj,y,varargin{:});
        %if there is no 'x', then this is a simple assignment of y
        if ~presence.x; return; end
      end
      %if 't' is present, assign it to 'x'
      if presence.t
        obj=assign@simpledata(obj,y,'x',obj.t2x(p.Results.t),varargin{:});
      end
      %update epoch (needed to derive obj.t from obj.x)
      %NOTICE: don't use obj.epoch= here, because at init that is not possible
      if ~isempty(p.Results.epoch)
        obj.epochi=p.Results.epoch;
      elseif presence.t
        obj.epochi=p.Results.t(1);
      else
        error('cannot derive epoch without either input ''epoch'' or ''t''.')
      end
      %update local records
      obj.step=simpletimeseries.timestep(obj.t);
      %sanitize (don't pass t, since it can be deliberatly changed)
      obj.check_st
    end
    function obj=copy_metadata(obj,obj_in,more_parameters,less_parameters)
      if ~exist('less_parameters','var')
        less_parameters={};
      end
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      obj=copy_metadata@simpledata(obj,obj_in,[simpletimeseries.parameters('list');more_parameters(:)],less_parameters);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      out=metadata@simpledata(obj,[simpletimeseries.parameters('list');more_parameters(:)]);
    end
    %the varargin method can be called directly
    %% info methods
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      relevant_parameters={'step','format','epoch','start','stop','timesystem'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print superclass
      print@simpledata(obj,tab)
    end
    function out=stats(obj,varargin)
      p=machinery.inputParser;
      p.addParameter('period', seconds(inf), @isduration);
      p.addParameter('overlap',seconds(0),   @isduration);
      p.addParameter('mode',  'struct',      @ischar);
      % parse it
      p.parse(varargin{:});
      % call upstream method if period is infinite
      if ~isfinite(p.Results.period)
        out=stats@simpledata(obj,varargin{:});
        return
      end
      % separate time series into segments
      ts=segmentedfreqseries.time_segmented(obj.t,p.Results.period,p.Results.overlap);
      % derive statistics for each segment
      s.msg=['deriving segment-wise statistics for ',str.clean(obj.descriptor,'file')]; s.n=numel(ts);
      for i=1:numel(ts)
        %call upstream procedure
        dat(i)=stats@simpledata(obj.trim(ts{i}(1),ts{i}(end)),varargin{:},'mode','struct');  %#ok<AGROW>
        % inform about progress
        s=time.progress(s,i);
      end
      % add time stamps
      for i=1:numel(ts)
        dat(i).t=mean(ts{i});
      end
      % unwrap data
      fn=fieldnames(dat);
      for i=1:numel(fn)
        %skip time
        if strcmp(fn{i},'t')
          continue
        end
        %resolving units
        switch lower(fn{i})
          case {'min','max','mean','std','rms','meanabs','stdabs','rmsabs'}
            units=obj.units;
          case {'length','gaps'}
            units=repmat({' '},1,obj.width);
        end
        out.(fn{i})=simpletimeseries(...
          transpose([dat.t]),...
          transpose(reshape([dat.(fn{i})],size(dat(1).(fn{i}),2),numel(dat))),...
          'format','datetime',...
          'labels',obj.labels,...
          'timesystem',obj.timesystem,...
          'units',units,...
          'descriptor',[fn{i},' ',str.clean(obj.descriptor,'file')]...
        );
      end
    end
    function out=stats2(obj1,obj2,varargin)
      p=machinery.inputParser;
      p.addParameter('period', seconds(inf), @isduration); %30*max([obj1.step,obj2.step])
      p.addParameter('overlap',seconds(0),   @isduration);
      % parse it
      p.parse(varargin{:});
      % call upstream method if period is infinite
      if ~isfinite(p.Results.period)
        out=stats2@simpledata(obj1,obj2,varargin{:});
        return
      end
      % separate time series into segments
      ts=segmentedfreqseries.time_segmented(...
        simpledata.union(obj1.t,obj2.t),...
        p.Results.period,...
        p.Results.overlap...
      );
      % derive statistics for each segment
      s.msg=['deriving segment-wise statistics for ',...
        str.clean(obj1.descriptor,'file'),' and ',...
        str.clean(obj2.descriptor,'file')...
      ]; s.n=numel(ts);
      for i=1:numel(ts)
        %call upstream procedure
        dat(i)=stats2@simpledata(...
          obj1.trim(ts{i}(1),ts{i}(end)),...
          obj2.trim(ts{i}(1),ts{i}(end)),...
          'mode','struct',varargin{:}...
        ); %#ok<AGROW>
        % inform about progress
        s=time.progress(s,i);
      end
      % add time stamps
      for i=1:numel(ts)
        dat(i).t=mean(ts{i});
      end
      % unwrap data and build timeseries obj
      fn=fieldnames(dat);
      for i=1:numel(fn)
        %skip time
        if strcmp(fn{i},'t')
          continue
        end
        %resolving units
        units=cell(1,obj1.width);
        for j=1:numel(units)
          switch lower(fn{i})
          case 'cov'
            units{j}=[obj1.units{j},'.',obj2.units{j}];
          case {'corrcoef','length'}
            units{j}=' ';
          end
        end
        out.(fn{i})=simpletimeseries(...
          [dat.t],...
          transpose(reshape([dat.(fn{i})],size(dat(1).(fn{i}),2),numel(dat))),...
          'format','datetime',...
          'labels',obj1.labels,...
          'timesystem',obj1.timesystem,...
          'units',units,...
          'descriptor',[fn{i},' ',str.clean(obj1.descriptor,'file'),'x',str.clean(obj2.descriptor,'file')]...
        );
      end
    end
    function out=str(obj)
      out=[datestr(obj.start),' -> ',datestr(obj.stop),' (',num2str(obj.nr_gaps),' gaps)'];
    end
    %% t methods
    function x_out=t2x(obj,t_now)
      if ~exist('t_now','var')
        t_now=obj.t;
      end
      if simpletimeseries.valid_t(t_now)
        x_out=simpletimeseries.time2num(t_now,obj.epoch,obj.x_units);
      else
        x_out=t_now;
      end
    end
    function t_out=x2t(obj,x_now)
      if ~exist('x_now','var')
        x_now=obj.x;
      end
      switch class(x_now)
      case 'datetime'; t_out=x_now;
      otherwise
        if simpledata.valid_x(x_now)
          t_out=simpletimeseries.num2time(x_now,obj.epoch,obj.x_units);
        else
          t_out=x_now;
        end
      end
    end
    function obj=set.t(obj,t_now)
      %NOTICE: this blindly changes the time domain!
      obj=obj.assign(obj.y,'t',t_now);
      obj.epoch=t_now(1);
    end
    function obj=set_t(obj,t_now)
      obj.t=t_now;
    end
    function out=get.t(obj)
      if isempty(obj.x)
        out=[];
      else
        out=obj.x2t(obj.x);
      end
    end
    function out=isx1zero(obj)
      %handle empty object
      if isempty(obj.x)
        out=true;
        return
      end
      %this function checks that:
      %if obj.x(1) is zero, then obj.epoch and obj.t(1) are equal
      test=[obj.x(1)==0,obj.start==obj.epoch];
      %sanity
      if test(1)~=test(2)
        error([...
          'obj.x(1)=',num2str(obj.x(1)),10,...
          'obj.start=',datestr(obj.start),10,...
          'obj.epoch=',datestr(obj.epoch),10,...
          'This combination is ilegal.'...
        ])
      end
      %outputs
      out=test(1);
    end
    function obj=t_reset(obj)
      %if needed, this function:
      % - resets obj.x, given the current obj.t, so that obj.x(1)=0
      % - recomputes obj.step
      obj=obj.step_update.epoch_update;
      %sanity
      obj.check_st
    end
    function out=span(obj)
      out=obj.stop-obj.start;
    end
    function out=t_domain(obj,step_now)
      if ~exist('step_now','var') || isempty(step_now)
        step_now=obj.step;
      end
      out=time.domain(obj.start,obj.stop,step_now);
    end
    function out=ishomogeneous(obj)
      htd=obj.t_domain;
      out=(numel(htd)==numel(obj.t)) && all(obj.t(:)==htd(:));
    end
    function out=istavail(obj,t)
      if isscalar(t)
        out=any(simpletimeseries.ist('==',t,obj.t,obj.t_tol));
      else
        out=arrayfun(@(i) obj.istavail(i),t);
      end
    end
    function out=isxavail(obj,x)
      %need to overload isxavail so that calls from simpledata come through
      %here and not through the function defined there
      out=obj.istavail(obj.x2t(x));
    end
    %NOTICE: to check for identical time domain, use obj.istequal (defined below)
    function obj=set.t_formatted(obj,t_now)
      [obj.t,format_now]=time.ToDateTime(t_now,obj.format);
      if ~strcmp(format_now,format_in)
        obj.format=format_now;
      end
      %sanitize
      obj.check_st(t_now)
    end
    function out=get.t_formatted(obj)
      out=time.FromDateTime(obj.t,obj.format);
    end
    function out=t_masked(obj,mask,mode)
      %NOTICE: do not use obj.t_masked(end), it is bugged. Use obj.t_masked([],'stop') instead
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      end
      if ~exist('mode','var') || isempty(mode)
        mode='all';
      end
      out=obj.t(mask);
      if isempty(out)
        out=time.zero_date;
      end
      switch mode
      case 'start'; out=out(1);
      case 'stop';  out=out(end);
      case 'all'    %do nothing
      otherwise; error(['unknown mode ''',mode,'''.'])
      end
    end
    function out=idx(obj,t_now,varargin)
      %need to handle doubles, to make it compatible with simpledata
      if isdatetime(t_now)
        out=idx@simpledata(obj,obj.t2x(t_now),varargin{:});
      else
        out=idx@simpledata(obj,t_now,varargin{:});
      end
    end
    function obj=at(obj,t_now,varargin)
      i=unique(obj.idx(t_now,varargin{:}));
      obj=obj.assign(...
        obj.y(i,:),...
        't',obj.t(i,:),...
        'mask',obj.mask(i,:)...
      );
    end
    function [obj,idx_add,idx_old,t_old]=t_merge(obj,t_add)
      %update epoch if needed (this is not really necessary, it just keeps x starting at zero)
      if t_add(1)<obj.start
        obj.epoch=t_add(1);
      end
      %call upstream method
      [obj,idx_add,idx_old,x_old]=obj.x_merge(obj.t2x(t_add));
      %convert outputs
      t_old=obj.x2t(x_old);
    end
    function out=mjd(obj,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=true(size(obj.t));
      end
      out=time.FromDateTime(obj.t(mask),'modifiedjuliandate');
    end
    function obj=addgaps(obj,min_gap_len)
      %find time domain indexes where explicit gaps need to be added
      gap_idx=find(diff(obj.t)>min_gap_len);
      %trivial call
      if isempty(gap_idx); return; end
      %determine mid-epoch of explicit gaps
      t_gaps=mean([obj.t(gap_idx),obj.t(gap_idx+1)],2);
      %augment object with newly found explicit gaps
      obj=obj.t_merge(t_gaps);
      % alternatively:
      % %compute mid-epoch for each gap
      % t_new=obj.t(gap_idx)+(obj.t(gap_idx+1)-obj.t(gap_idx))/2;
      % %add NaNs
      % obj=obj.set_at(t_new,nan(numel(t_new),obj.width));
    end
    function sanity_t(obj)
      assert(obj.isxequal(simpletimeseries.time2num(obj.t,obj.epoch,obj.x_units)),...
        'Discrepancy between t-domain and x-domain, debug needed!')
    end
    function assert_t_domain(obj1,obj2)
      obj1.sanity_t;
      obj2.sanity_t;
      if ~obj1.istequal(obj2)
        disp([...
          'Cannot operate on scalar objects with different t-domains: ',newline,...
          'obj1.t obj2.t',...
          ])
        idx=obj1.peek_idx;
        for i=1:numel(idx)
          disp(datestr([obj1.t(i),obj2.t(i)]))
        end
        error('Cannot continue')
      end
    end
    %% step methods
    function out=step_num(obj)
      out=simpletimeseries.time2num(obj.step,0,obj.x_units);
    end
    function out=step_get(obj)
      out=simpletimeseries.timestep(obj.t);
    end
    function obj=step_update(obj)
      obj.step=simpletimeseries.timestep(obj.t);
    end
    function out=step_gcd(obj1,obj2)
      out=simpletimeseries.timescale(gcd(...
        simpletimeseries.timescale(obj1.step),...
        simpletimeseries.timescale(obj2.step)...
      ));
    end
    function out=step_lcm(obj1,obj2)
      if obj1.step==0 || obj2.step==0
        out=1;
      else
        out=simpletimeseries.timescale(lcm(...
          simpletimeseries.timescale(obj1.step),...
          simpletimeseries.timescale(obj2.step)...
        ));
      end
    end
    function [obj1,obj2]=match_step(obj1,obj2)
      %sanity
      if ~obj1.ishomogeneous || ~obj2.ishomogeneous
        error('can only handle homogeneous time domains.')
      end
      %trivial call
      if simpletimeseries.ist('==',obj1.step,obj2.step,min([obj1.t_tol,obj2.t_tol]))
        return
      end
      %new timestep is the greatest common divisor
      step_now=step_gcd(obj1,obj2);
      if obj1.debug || obj2.debug
        disp(['WARNING: Reset step in ',...
          'obj1 (',obj1.descriptor,') from ',char(obj1.epoch),' and ',...
          'obj2 (',obj2.descriptor,') from ',char(obj2.epoch),...
          'to ',char(step_now)])
      end
      %resample to the common step size
      obj1=obj1.resample(step_now);
      obj2=obj2.resample(step_now);
    end
    %% epoch methods
    function obj=set.epoch(obj,epoch)
      assert(simpletimeseries.valid_epoch(epoch),'invalid input ''epoch''.')
      %get current time domain
      t_old=obj.t;
      %set epoch
      obj.epochi=epoch;
      %shift x
      obj=obj.assign_x(simpletimeseries.time2num(t_old,epoch,obj.x_units));
      %sanity
      assert(obj.istequal(t_old),'changing epoch caused the time domain to also change.')
    end
    function out=get.epoch(obj)
      out=obj.epochi;
    end
    function obj=epoch_update(obj)
      obj.epoch=obj.t(1);
    end
    function [obj1,obj2]=match_epoch(obj1,obj2)
      %trivial call
      if simpletimeseries.ist('==',obj1.epoch,obj2.epoch,min([obj1.t_tol,obj2.t_tol]))
        return
      end
      %match epochs, avoid zero and inf epochs
      if time.isfinite(obj1.epoch)
        if obj1.debug || obj2.debug
          disp(['WARNING: Reset epoch in obj2 (',obj2.descriptor,') to ',datestr(obj1.epoch),' from ',datestr(obj2.epoch)])
        end
        obj2.epoch=obj1.epoch;
      else
        if obj1.debug || obj2.debug
          disp(['WARNING: Reset epoch in obj1 (',obj1.descriptor,') to ',datestr(obj2.epoch),' from ',datestr(obj1.epoch)])
        end
        obj1.epoch=obj2.epoch;
      end
    end
    %% start/stop methods
    function out=get.start(obj)
      out=obj.t(1);
    end
    function out=get.first(obj)
      out=obj.t_masked([],'start');
    end
    function out=get.stop(obj)
      out=obj.t(obj.length);
    end
    function out=get.last(obj)
      out=obj.t_masked([],'stop');
    end
    function obj=set.start(obj,start)
      if isempty(start) || ...                                        %ignore empty inputs
          simpletimeseries.ist('==',start,obj.start,obj.t_tol) || ... %shortcut
      all(simpletimeseries.ist('==',obj.t,time.zero_date,obj.t_tol))  %do not operate on static-esque objects
        %trivial call
        return
      %check if required start is before the start of the current time series
      elseif simpletimeseries.ist('<',start,obj.start,obj.t_tol)
        %preppend a single epoch
        obj=obj.append_epochs(start,nan);
      %check if required start is after the end of the current time series
      elseif simpletimeseries.ist('>',start,obj.stop,obj.t_tol)
        %build pseudo-empty data (all previous data is discarded)
        obj=obj.assign(nan(1,obj.width),'t',start);
      else
        %trim object
        obj=obj.trim(start,obj.stop);
        %recursive call in case start is at the middle of an epoch
        if ~simpletimeseries.ist('==',start,obj.start,obj.t_tol)
          %preppend a single epoch
          obj=obj.append_epochs(start,nan);
        end
      end
    end
    function obj=set.stop(obj,stop)
      if isempty(stop) || ...                                        %ignore empty inputs
          simpletimeseries.ist('==',stop,obj.stop,obj.t_tol) || ...  %shortcut
      all(simpletimeseries.ist('==',obj.t,time.zero_date,obj.t_tol)) %do not operate on static-esque objects
        %trivial call
        return
      %check if required stop is after the end of the current time series
      elseif simpletimeseries.ist('>',stop,obj.stop,obj.t_tol)
        %append a single epoch
        obj=obj.append_epochs(stop,nan);
      %check if required stop is before the start of the current time series
      elseif simpletimeseries.ist('<',stop,obj.start,obj.t_tol)
        %build pseudo-empty data (all previous data is discarded)
        obj=obj.assign(nan(1,obj.width),'t',stop);
      else
        %trim object
        obj=obj.trim(obj.start,stop);
        %recursive call in case start is at the middle of an epoch
        if ~simpletimeseries.ist('==',stop,obj.stop,obj.t_tol)
         %append a single epoch
          obj=obj.append_epochs(stop,nan);
        end
      end
    end
    %% tsys methods
    function out=get.tsys(obj)
      out=obj.timesystem;
    end
    function obj=set.tsys(obj,in)
      if ~simpletimeseries.valid_timesystem(in)
        error(['need a valid time system, i.e. one of ',strjoin(simpletimeseries.valid_timesystems,', '),'.'])
      end
      if ~strcmpi(obj.timesystem,in)
        obj.t=time.([obj.timesystem,'2',lower(in)])(obj.t);
        obj.timesystem=in;
      end
    end
    %% management methods
    function check_st(obj,t_now)
      %check consistency in the values of obj.start and obj.epoch
      obj.isx1zero;
      %check for monotonously increasing time domain
      if any(diff(obj.x)<=0)
        error('the time domain is not monotonously increasing.')
      end
      if exist('t_now','var') && ~isempty(t_now)
        %check for consistency in the time domain
        assert(obj.istequal(t_now),'the time domain is not consistent with input ''t_now''.')
      end
    end
    %% edit methods (overloaded with simpledata)
    %the remove method can be called directly
    function obj=trim(obj,start,stop)
      %do not operated on static-esque objects
      if obj.length==1 && simpletimeseries.ist('==',obj.t,time.zero_date,obj.tol)
        return
      end
      obj=trim@simpledata(obj,obj.t2x(start),obj.t2x(stop));
    end
    function obj=slice(obj,start,stop)
      obj=slice@simpledata(obj,obj.t2x(start),obj.t2x(stop));
    end
    function obj=interp(obj,t_now,varargin)
      %convert duration to numeric
      varargin=simpletimeseries.fix_interp_over_gaps_narrower_than(varargin);
      %call superclass
      obj=interp@simpledata(obj,obj.t2x(t_now),varargin{:});
      %update step
      obj=obj.t_reset;
    end
    function [obj,S]=polyfit(obj,order,t_now)
      %call superclass
      [obj,S]=polyfit@simpledata(obj,order,obj.t2x(t_now));
      %update step
      obj=obj.t_reset;
    end
    %the detrend method can be called directly
    %the outlier method can be called directly
    function obj=mean(obj,n,keep_time_domain)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      if ~exist('keet_time_domain','var')
        obj=obj.segstat(@mean,n);
      else
        obj=obj.segstat(@mean,n,keep_time_domain);
      end
    end
    function obj=median(obj,n,keep_time_domain)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      if ~exist('keet_time_domain','var')
        obj=obj.segstat(@median,n);
      else
        obj=obj.segstat(@median,n,keep_time_domain);
      end
    end
    function obj=rms(obj,n,keep_time_domain)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      if ~exist('keet_time_domain','var')
        obj=obj.segstat(@rms,n);
      else
        obj=obj.segstat(@rms,n,keep_time_domain);
      end
    end
    function obj=srd(obj,n,keep_time_domain)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      if ~exist('keet_time_domain','var')
        obj=obj.segstat(@srd,n);
      else
        obj=obj.segstat(@srd,n,keep_time_domain);
      end
    end
    function obj=segstat(obj,op,n,keep_time_domain)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      if ~exist('keet_time_domain','var') || isempty(keep_time_domain)
        keep_time_domain=false;
      end
      if keep_time_domain
        %save current time domain
        t_now=obj.t;
      end
      %handle periods
      if isduration(n)
        %compute (average) number of epochs within the requested t_span
        n=round(n/obj.step);
      end
      %trivial call: ignore irrelevant spans
      if n <= 1
        return
      end
      %call superclass
      obj=segstat@simpledata(obj,op,n);
      if keep_time_domain
        %resample (if needed, which is checked inside resample)
        obj=obj.interp(t_now,...
          'interp_over_gaps_narrower_than',0,...
          'interp1_args',{'linear'}...
        );
      end
    end
    %% edit methods (specific to this class)
    function obj=extend(obj,nr_epochs)
%       %sanity
%       if ~obj.ishomogeneous
%         error(['cannot handle non-homogeneous time domains.'])
%       end
      switch class(nr_epochs)
      case 'double'
        if nr_epochs==0
          return
        end
        if (nr_epochs~=round(nr_epochs))
          error(['input ''nr_epochs'' must be an integer, not ',num2str(nr_epochs),'.'])
        end
        %define
        if nr_epochs>0
          %extend
          t_new=[obj.t;transpose(obj.stop+obj.step:obj.step:obj.stop+nr_epochs*obj.step)];
          y_new=[obj.y;nan(nr_epochs,obj.width)];
        else
          nr_epochs=-nr_epochs;
          %prepend
          t_new=[transpose(obj.start-nr_epochs*obj.step:obj.step:obj.start-obj.step);obj.t];
          y_new=[nan(nr_epochs,obj.width);obj.y];
        end
        %propagate
        obj=obj.assign(y_new,'t',t_new);
      case 'datetime'
        t_now=nr_epochs;
        if t_now <obj.start
          t_ref=obj.start;
        elseif t_now> obj.stop
          t_ref=obj.stop;
        elseif t_now==obj.start || t_now==obj.stop
          %do nothing
          return
        else
          error(['input ''t'' (',datestr(t_now),') ',...
            'must be larger than obj.stop (',datestr(obj.stop),') ',...
            'or smaller than than obj.start (',datestr(obj.start),').'...
          ]);
        end
        obj=extend(obj,floor((t_now-t_ref)/obj.step));
      otherwise
        error(['cannot handle input ''nr_epochs'' of class ',class(nr_epochs),'.'])
      end
    end
    function obj=extend_or_trim_end(obj,t_end)
      assert(isdatetime(t_end),['Input ''t_cut'' must of of class datetime, not ',class(t_end)])
      %check if t_ref is outside obj timespan
      if t_end < obj.start
        error(['Expecting input t_ref (',datestr(t_end),') to be larger than obj.start (',datestr(obj.start),').'])
      elseif obj.stop <= t_end
        %extend
        obj=obj.extend(t_end);
      else
        %trim end section
        obj=obj.trim(obj.start,t_end);
      end
    end
    function obj=extend_or_trim_start(obj,t_start)
      assert(isdatetime(t_start),['Input ''t_cut'' must of of class datetime, not ',class(t_start)])
      %check if t_ref is outside obj timespan
      if obj.stop < t_end
        error(['Expecting input t_ref (',datestr(t_end),') to be smaller than obj.stop(',datestr(obj.stop),').'])
      elseif t_start <= obj.start
        %extend
        obj=obj.extend(t_start);
      else
        %trim start section
        obj=obj.trim(t_start,obj.stop);
      end
    end
    function obj=append_epochs(obj,t_now,y_now)
      %trivial call
      if isempty(t_now)
        return
      end
      %shortcuts
      if isscalar(y_now)
        y_now=y_now*ones(numel(t_now),obj.width);
      end
      %sanity
      assert(numel(t_now)==size(y_now,1) || size(y_now,2)~=obj.width,...
        'inputs ''t_now'' and/or ''y_now'' have sizes inconsistent with this obj.')
      %try to avoid sorting
      if all(t_now<obj.start)
        %preppend the data (t_now should be sorted)
        y_now=[y_now;obj.y];
        t_now=[t_now(:);obj.t];
      elseif all(t_now>obj.stop)
        %append the data (t_now should be sorted)
        y_now=[obj.y;y_now];
        t_now=[obj.t;t_now(:)];
      else
        %concatenate data
        y_now=[obj.y;y_now];
        t_now=[obj.t;t_now(:)];
        %sort along time
        [t_now,sort_idx]=sort(t_now);
        y_now=y_now(sort_idx,:);
      end
      %append the epoch
      obj=obj.assign(y_now,'t',t_now);
    end
    function [obj,idx]=fill(obj)
      %NOTICE: this method is similar to resample in the sense it creates a complete time domain
      %        but it differs because the added time entries are set as explicit gaps.
      %TODO: handle non-homogeneous time domains
      %trivial call
      if obj.ishomogeneous
        if nargout > 1, idx=true(obj.length,1);end
        return
      end
      %build complete time domain
      t_new=obj.t_domain;
      t_old=obj.t;
      % sanity
      if numel(t_new) < numel(t_old)
        error('BUG TRAP: complete time domain has less entries than current time domain. Debug needed!')
      end
      %find out where there are gaps larger than the step size
      gap_idx=find(diff(obj.t)>obj.step);
      %if there are no gaps and the time series is not homogeneous, we have a problem that needs fixing
      if isempty(gap_idx)
        error('implementation needed!')
      end
      disp(['Need to fill in missing epochs: ',num2str(numel(t_new)-obj.length),' ('...
        num2str((numel(t_new)-obj.length)/numel(t_new)*1e2),'%).'])
      %loop over all implicit gaps (i.e. missing epochs)
      s.msg=[' populating missing epochs (',datestr(obj.start),' to ',datestr(obj.stop),')',...
        ' of ',obj.descriptor];s.n=numel(gap_idx);
      while ~isempty(gap_idx)
        %create patch
        t_patch=transpose((obj.t(gap_idx(1))+obj.step):obj.step:(obj.t(gap_idx(1)+1)-obj.step));
        %if t_patch is empty, then this loop goes forever
        if isempty(t_patch)
          t_patch=obj.t(gap_idx(1))+obj.step;
        end
        %save data with patch (it gets deleted when assigning to x)
        y_patched=[obj.y(1:gap_idx(1),:);...
                   nan(numel(t_patch),obj.width);...
                   obj.y(gap_idx(1)+1:end,:)];
        %create patched t
        t_patched=[obj.t(1:gap_idx(1));...
                  t_patch;...
                  obj.t(gap_idx(1)+1:end)];
        %propagate y
        obj=obj.assign(y_patched,'t',t_patched);
        %re-discover gaps
        gap_idx=find(diff(obj.t)>obj.step);
        %user feedback
        s=time.progress(s);
      end
      %trim end (there might be a single epoch dangling at the end)
      obj.stop=t_new(end);
      %sanitize
      obj.check_st(t_new);
      %additional output arguments
      if nargout > 1
        [~,idx]=simpledata.union(t_old,t_new);
      end
    end
    function obj=resample(obj,step_now)
      % this function is a special case of interpolation
      if ~exist('step_now','var') || isempty(step_now)
        step_now=obj.step_get;
      end
      if ~isduration(step_now)
        error(['expecting input ''step_now'' to be duration, not ',class(step_now),'.'])
      end
      % build/retrieve relevant time domain
      t_now=obj.t_domain(step_now);
      % trivial call
      if numel(obj.t)==numel(t_now) && all(obj.t==t_now)
        return
      end
      % interpolate over new time domain
      obj=obj.interp(t_now,...
        'interp_over_gaps_narrower_than',3*step_now,...
        'interp1_args',{'linear'}...
      );
    end
    function obj=resample_full(obj)
      % this function is a special case of resample
      obj=obj.interp(obj.t_domain,...
        'interp_over_gaps_narrower_than',0,...
        'interp1_args',{'linear'}...
      );
    end
    function obj=fstep(obj,step_prev)
      %adds data entries that are equal to the preceeding value, but one
      %step_prev before the following epoch (also for explicit gaps)
      obj_new=simpletimeseries(...
               obj.t(   2:end    )-step_prev,... %time domain is the time domain of obj shifted by step_prev
               obj.y(   1:end-1,:),...
        'mask',obj.mask(1:end-1),...
        'format','datetime',...
        'timesystem',obj.timesystem...
      );
      %merge the two objects
      obj=obj.augment(obj_new);
    end
    function [obj_clean,obj_outlier]=despike(obj,n,varargin)
      if ~exist('n','var') || isempty(n)
        n=ceil(obj.length*0.05);
      end
      %get medianed timeseries
      obj_median=obj.median(n);
      %compute residual to median
      obj_res=obj-obj_median;
      %remove outliers from residual
      [obj_res_clean,obj_res_outlier]=obj_res.outlier(varargin{:});
      %restore median
      obj_clean=obj_median+obj_res_clean;
      obj_outlier=obj_median+obj_res_outlier;
    end
    function obj=repeat(obj,tf)
      while obj.stop < tf
        %append data translated in time by obj.span+obj.step
        obj=obj.append(...
          obj.set_t(obj.t+obj.span+obj.step)...
        );
      end
      %crop excess
      if obj.stop>tf
        obj=obj.trim(obj.start,tf);
      end
    end
    %% multiple object manipulation
    function out=istequal(obj1,obj2)
      %NOTICE: this also handles the single-object operation
      if isdatetime(obj2)
        out=obj1.length==numel(obj2) & ~any(~simpletimeseries.ist('==',obj1.t,obj2,obj1.t_tol));
      else
        out=obj1.length==obj2.length & ~any(~simpletimeseries.ist('==',obj1.t,obj2.t,min([obj1.t_tol,obj2.t_tol])));
      end
    end
    function compatible(obj1,obj2,varargin)
      p=machinery.inputParser;
      p.addParameter('skip_par_check',{''},@iscellstr)
      p.parse(varargin{:});
      %call mother routine
      compatible@simpledata(obj1,obj2,varargin{:});
      %shorter names
      par=simpletimeseries.compatible_parameter_list;
      for i=1:numel(par)
        % if a parameter is empty, no need to check it
        if all(cells.isempty(obj1.(par{i})))|| all(cells.isempty(obj2.(par{i})))
          continue
        end
        if ~cells.isincluded(p.Results.skip_par_check,par{i}) && ~isequal(obj1.(par{i}),obj2.(par{i}))
          error(['discrepancy in parameter ',par{i},'.'])
        end
      end
    end
    %the merge method can be called directly
    function [obj1,obj2]=interp2(obj1,obj2,varargin)
      %trivial call
      if obj1.istequal(obj2)
        return
      end
      %extends the t-domain of both objects to be in agreement
      %with the each other. The resulting t-domains possibly have
      %numerous gaps, which are interpolated over (interpolation
      %scheme and other options can be set in varargin).
      %handle default optional arguments
      if ~exist('varargin','var') || isempty(varargin)
        varargin={...
          'interp_over_gaps_narrower_than',3*min([obj1.step,obj2.step]),...
          'interp1_args',{'linear'}...
        };
      end
      %call upstream method
      [obj1,obj2]=interp2@simpledata(obj1,obj2,varargin{:});
      %sanity
      obj1.assert_t_domain(obj2);
    end
    %the append method can be called directly
    %the augment method can be called directly
    %the glue method can be called directly
    function [obj1,obj2]=interp2_lcm(obj1,obj2,varargin)
      %NOTICE: this function used to be called consolidate
      %extends the time domain of both objects to be in agreement
      %with the each other
      [obj1,obj2]=obj1.match_tx_domain(obj2);
      obj1.compatible(obj2,varargin{:})
      %trivial call
      if istequal(obj1,obj2)
        return
      end
      %build extended time domain, with lcm timestep, rounded to the nearest second
      t_now=dateshift(min([obj1.start,obj2.start]),'start','second'):...
           step_lcm(obj1,obj2):...
           dateshift(max([obj1.stop, obj2.stop]),  'end',  'second');
%       h=figure;
%       obj1.plot('column',1,'line',{'o-'}), hold on
%       obj2.plot('column',1,'line',{'x-'})
      %interpolate to new time domain
      obj1=obj1.interp(t_now,'interp_over_gaps_narrower_than',3*obj1.step,'interp1_args',{'spline'});
      obj2=obj2.interp(t_now,'interp_over_gaps_narrower_than',3*obj2.step,'interp1_args',{'spline'});
%       figure(h)
%       obj1.plot('column',1,'line',{'*-'}), hold on
%       obj2.plot('column',1,'line',{'+-'})
%       legend('o1 original','o2 original','o1 interp','o2 interp')
    end
    %% calibration
    function obj1=calibrate_poly(obj1,obj2,order)
      %need to match the epoch
      if isa(obj1,'simpletimeseries') && isa(obj2,'simpletimeseries')
        [obj1,obj2]=obj1.match_epoch(obj2);
      end
      if ~exist('order','var') || isempty(order)
        order=1;
      end
      %call mother routine
      obj1=calibrate_poly@simpledata(obj1,obj2,order);
    end
    %% wrappers
    function obj=smooth(obj,span,varargin)
      %handle periods
      if isduration(span)
        %compute (average) number of epochs within the requested t_span
        span=round(span/obj.step);
      end
      %trivial call: ignore irrelevant spans
      if span <= 1
        return
      end
      %call mother routine
      obj=smooth@simpledata(obj,span,varargin{:});
    end
    %% frequency analysis
    function [obj,filter_response]=bandpass(obj,T,varargin)
      p=machinery.inputParser;
      % add stuff as needed
      p.addRequired('T',                 @(i) isduration(i) && numel(i)==2 && i(1)>i(2));
      p.addParameter('gaps',  'zeroed',  @ischar);
      p.addParameter('debug_plot',false, @(i) islogical(i) && isscalar(i));
      p.addParameter('soft_transition_radius', 0.1, @num.isscalar);
      % parse it
      p.parse(T,varargin{:});
      %handle gaps
      switch p.Results.gaps
      case 'trunc'
        %truncating bad data (not a good idea)
        data_in=obj.y_masked;
      case 'zeroed'
        %zeroing bad data
        data_in=obj.y;
        data_in(~obj.mask,:)=0;
      otherwise
        error(['unknown gap handling mode ''',p.Results.gaps,'''.'])
      end
      %sanity
      if any(isnan(data_in(:)))
          error('found NaNs in the input data.')
      end
      %computational length
      n = 2^nextpow2(size(data_in,1));
      %build long filter domain
      ff=1/obj.step_num/2*linspace(0,1,n/2);
      fP=zeros(size(ff));
      %convert to Herts
      Wn=1./seconds(T);
      disp(str.show({'bandpass periods    : [',T,']'}))
      disp(str.show({'bandpass frequencies: [',Wn,'] Hz'}))
      %assign filter factors
      fP(ff>=Wn(1) & ff<=Wn(2))=1;
      %trivial call: cut-off frequencies are outside nyquist and lowest frequency (e.g. [inf,0])
      if sum(fP)==0; return; end
      %parameters (min is needed in case the pass-band is very wide)
      smooth_radius=min([...
        ceil(sum(~fP)*p.Results.soft_transition_radius),...
        ceil(sum( fP)*p.Results.soft_transition_radius)...
      ]); %data points
      %smooth transitions
      idx={...
        find(ff<Wn(1),1,'last'),...
        find(ff>Wn(2),1,'first')...
      };
      % figure
      % semilogx(ff,fP), hold on
      for i=1:2
        if ~isempty(idx{i})
          idx_out=(idx{i}-smooth_radius):(idx{i}+smooth_radius+1);
          idx_in =[idx_out(1),idx_out(end)];
          fP(idx_out)=spline(ff(idx_in),[0 fP(idx_in) 0],ff(idx_out));
        end
      end
      % semilogx(ff,fP), hold on
      % keyboard
      %mirror the filter
      fP=[fP,fliplr(fP)];
      %apply the filter
      fX=fft(data_in,n).*(fP(:)*ones(1,size(data_in,2)));
      fx=ifft(fX,'symmetric');
      %trim excess
      fx=fx(1:size(data_in,1),:);
      if p.Results.debug_plot
        m=numel(ff);
        X=fft(data_in(:,1),n);
        PX=X(1:m).*conj(X(1:m));
        PfX=fX(1:m,1).*conj(fX(1:m,1));
        figure
        subplot(2,1,1)
        title('frequency domain')
        loglog(ff,PX), hold on
        loglog(ff,PfX)
        loglog(ff,fP(1:m)*max([max(PX),max(PfX)]))
        legend('original','filtered','filter')

        subplot(2,1,2)
        title('time domain')
        plot(fx(:,1)), hold on
        plot(obj.y(:,1))
        legend('filtered','original')
        keyboard
      end
      %propagate
      obj=obj.assign(fx,'t',obj.t,'mask',obj.mask);
      %additional outputs
      if nargout>1
        filter_response.f=ff;
        filter_response.a=fP(1:numel(ff));
      end
    end
    %% derivative
    function obj=deriv(obj,varargin)
      %TODO: test this

      %get data
      y=obj.y;
      %call primitive
      %NOTICE: these input arguments are here to set a reasonable default for the important parameters of num.diff
      dy=num.diff(y,obj.step,'npoints',5,'nderiv',1,'extremeties','nan',varargin{:});
      %back propagate
      obj=obj.assign(dy);
    end
    %% segment analysis
    function out=segmentate(obj,seg_length,seg_overlap,varargin)
      if ~isduration(seg_length)
        error(['input ''seg_length'' must be of class ''duration'', not ''',class(seg_length),'''.'])
      end
      if ~isduration(seg_overlap)
        error(['input ''seg_overlap'' must be of class ''duration'', not ''',class(seg_overlap),'''.'])
      end
      if seg_overlap>=seg_length
        error(['input ''seg_overlap'' (',num2str(seg_overlap),') must be smaller than input ''seg_length'' (',num2str(seg_length),').'])
      end
      %handle infinite segment length
      if ~isfinite(seg_length)
        seg_length=obj.span;
        seg_overlap=seconds(0);
      end
      %guess the number of segments
      n=ceil((obj.span)/seg_length*2);
      %init outputs
      out=cell(1,n);
      %init loop vars
      ti=obj.start;
      tf=ti;
      c=0;
      %loop it
      while tf < obj.stop
        %increment counter
        c=c+1;
        %determine stop t
        tf=ti+seg_length;
        %build segment
        out{c}=obj.interp(ti:obj.step:tf);
        %update descriptor
        out{c}.descriptor=['segment ',num2str(c),' of ',obj.descriptor];
        %set start t for next iter
        ti=tf-seg_overlap;
      end
      %remove empty entries
      out=out(~cells.isempty(out));
    end
    %1. splits the time series into segments
    %2. computes the mean of all those segements
    %3. repeats that mean segment over the complete time domain (discontinuities likely)
    function [out,res,J,segs]=component_ampl(obj,seg_length,tf,varargin)
      %defaul value for some inputs
      if ~exist('tf','var') || isempty(tf)
        tf=obj.stop;
      end
      %segmentate
      segs=obj.segmentate(seg_length,obj.step,varargin{:});
      %build output time domain (done after getting the segments of the unextended object)
      obj=obj.extend_or_trim_end(tf).resample;
      %init output
      m=zeros(segs{1}.length,obj.width);
      %traverse data width
      for ci=1:obj.width
        %make room for segment data
        ynow=nan(segs{1}.length,numel(segs));
        %agregate segments
        for si=1:numel(segs)
          if segs{si}.length==size(ynow,1)
            ynow(:,si)=segs{si}.y(:,ci);
          else
            ynow(1:segs{si}.length,si)=segs{si}.y(:,ci);
          end
        end
        %compute mean and std
        m(:,ci)=mean(ynow, 2,'omitnan');
      end
      %build timeseries objects
      switch class(segs{1})
      case 'simpletimeseries'
        out=simpletimeseries(segs{1}.t,m,obj.varargin{:});
      case 'gravity'
        out=gravity(segs{1}.t,m,obj.varargin{:},varargin{:});
      otherwise
        error(['Class ',class(segs{1}),' not yet implemented (easy fix!)'])
      end
      out.descriptor=['mean amplitude of ',char(seg_length),' period of ',obj.descriptor];
      %repeat onto input (and extended) time domain
      out=out.repeat(tf).interp(obj.t);
      if nargout > 1
        %compute residuals
        res=out-obj;
      end
      if nargout > 2
        %compute relative norm of residuals (scalar)
        J=res.stats('mode','norm')./obj.stats('mode','norm');
      end
    end
    function component_split_plot(stats)
      plotting.figure('plot_visible',false);
      legend_str=cell(0); i=0;
      for j=1:numel(stats.periodic)
        stats.periodic{j}.plot('column',1,'line',{'--'});
        i=i+1;legend_str{i}=['period ',str.show(stats.parameters.seg_lengths(j))];
        stats.aperiodic{j}.plot('column',1,'line',{'--'});
        i=i+1;legend_str{i}=[' poly',num2str(stats.parameters.aperiodic_order),' ',str.show(stats.parameters.seg_lengths(j))];
      end
       stats.out.plot('column',1);i=i+1;legend_str{i}='out';
      stats.rest.plot('column',1);i=i+1;legend_str{i}='rest';
      stats.orig.plot('column',1);i=i+1;legend_str{i}='original';
      plotting.legend('plot_legend',legend_str,'plot_legend_location','westoutside');
      plotting.font_size;
      plotting.line_color;
      plotting.line_width;
    end
    function f=component_split_fileroot(obj,dir,ext,search,lower_durat,upper_durat,nr_seg,t_extrapolate,t_crop,aperiodic_order)
      %build filename in case it needs to be saved outside
      f=fullfile(dir,[...;
        str.show({...
          obj.descriptor,...
          str.show({'search',            search         },'',          '_'),...
          strrep(str.show({'lower_durat',lower_durat    },'',          '_'),' ',''),...
          strrep(str.show({'upper_durat',upper_durat    },'',          '_'),' ',''),...
          str.show({'nr_seg',            nr_seg         },'',          '_'),...
          str.show({'t_extrapolate',     t_extrapolate  },'yyyy-mm-dd','_'),...
          str.show({'t_crop',            t_crop         },'yyyy-mm-dd','_'),...
          str.show({'aperiodic_order',   aperiodic_order},'',          '_'),...
          strjoin(obj.labels,'_')...
        },'join_char','.'),...
      ext]);
    end
    function J=component_split_fitness(obj,seg_length,tf,varargin)
      [~,~,J]=obj.component_ampl(seg_length,tf,varargin{:});
    end
    function [obj,stats]=component_split(obj,varargin)
      v=varargs.wrap('sources',{{...
        'search',         true,       @islogical;...
        'lower_durat',    obj.step*4, @isduration;...
        'upper_durat',    obj.span/2, @isduration;...
        'nr_seg',         3,          @num.isscalar;...
        't_extrapolate',  obj.stop,   @isdatetime;...
        't_crop' ,        obj.stop,   @(i) isdatetime(i) && i<=obj.stop;...
        'aperiodic_order',3,          @num.isscalar;...
        'plot_dir',       '',         @ischar;...
        'data_dir',       '',         @ischar;...
        'time_scale',     @years,     @(i) isa(i,'function_handle') && isduration(i(1));...
      }},varargin{:});
      v=varargs.wrap('sources',{{...
        'seg_lengths',    linspace(v.lower_durat,v.upper_durat,v.nr_seg), @isduration;...
      },v},varargin{:});
      %update some parameters, in case seg_lengths was given explicitly
      v.nr_seg=numel(v.seg_lengths);
      %check if data is already available (need data_dir non-empty)
      recompute=true;
      if ~isempty(v.data_dir)
        datafile=obj.component_split_fileroot(v.data_dir,'.mat',...
          v.search,v.lower_durat,v.upper_durat,v.nr_seg,v.t_extrapolate,v.t_crop,v.aperiodic_order);
        if exist(datafile,'file')
          str.say('Loading data file ',datafile)
          load(datafile,'obj','obj_orig','stats')
          recompute=false;
        end
      end
      %maybe don't need to recompute stuff
      if recompute
        %enforce t_crop
        obj=obj.trim(obj.start,v.t_crop);
        %save original object
        obj_orig=obj;
        std_orig=obj.stats('mode','std');
        %inits
        periodic=cell(1,v.nr_seg);
        aperiodic=cell(1,v.nr_seg);
        S=cell(1,v.nr_seg);
        %init rest
        rest=obj;
        %inform user
        str.say(obj.descriptor,' is going to be decomposed with the following options:');disp(v.show)
        %iterate
        for i=1:v.nr_seg
          if v.search
            str.say('For',obj.descriptor,'estimating',str.th(i),'out of',v.nr_seg,'periodic components')
          else
            str.say('For',obj.descriptor,'computing',v.nr_seg,'periodic components with lengths',str.show({v.seg_lengths(:)},'join_char',', '))
          end
          %compute polifit and remove from rest (unextended time domain)
          warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
          [aperiodic{i},S{i}]=rest.polyfit(v.aperiodic_order,rest.t);
          warning('on','MATLAB:polyfit:RepeatedPointsOrRescale')
          %decumulate aperiodic component (unextended time domain)
          rest=rest-aperiodic{i}.interp(obj.t);
          %search
          if v.search
            %define function to minimze
            fun=@(i) rest.component_split_fitness(v.time_scale(i),v.t_extrapolate,varargin{:});
            %get segment length for this iter
            v.seg_lengths(i)=v.time_scale(fminbnd(fun,...
              v.time_scale(v.lower_durat),...
              v.time_scale(v.upper_durat),...
              optimset('Display','off')...
            ));
            str.say(str.th(i),'periodic component has length',v.seg_lengths(i))
          end
          %get this periodic component (extended time domain)
          periodic{i}=rest.component_ampl(v.seg_lengths(i),v.t_extrapolate,varargin{:});
          %decumulate (unextended time domain)
          rest=rest-periodic{i}.interp(obj.t);
        end
        %accumulate (extended time domain)
        out=periodic{1}+obj.polyfit(S{1},periodic{1}.t);
        for i=2:v.nr_seg
          out=out+periodic{i}+obj.polyfit(S{i},periodic{1}.t);
        end
        %update output
        obj=out;
        %save data and std ratios
        stats=struct(...
          'parameters', v,...
          'periodic',   {periodic},...
          'aperiodic',  {aperiodic},...
          'rest',       rest,...
          'out',        out,...
          'orig',       obj_orig,...
          'stdr',struct(...
            'periodic', cellfun(@(i) i.stats('mode','std')./std_orig, periodic),...
            'aperiodic',cellfun(@(i) i.stats('mode','std')./std_orig,aperiodic),...
            'rest',               rest.stats('mode','std')./std_orig,...
            'out',                 out.stats('mode','std')./std_orig...
          )...
        );
        %inform
        str.say('periods:',stats.parameters.seg_lengths,'std ratios:',...
          'periodic', sqrt(sum(stats.stdr.periodic.^2)),...
          'aperiodic',sqrt(sum(stats.stdr.aperiodic.^2)),...
          'rest',     stats.stdr.rest,...
          'out',      stats.stdr.out...
        );
        %save it
        if ~isempty(v.data_dir)
          file.ensuredir(datafile);
          str.say('Saving data file ',datafile)
          save(datafile,'obj','obj_orig','stats')
        end
      end
      %plot it
      if ~isempty(v.plot_dir)
        obj.component_split_plot(stats);
        plotfile=obj_orig.component_split_fileroot(v.data_dir,'.png',...
          v.search,v.lower_durat,v.upper_durat,v.nr_seg,v.t_extrapolate,v.t_crop,v.aperiodic_order);
        file.ensuredir(plotfile);
        str.say('Saving plot ',plotfile)
        saveas(gcf,plotfile)
        close(gcf)
      end
    end
    %% parametric decomposition
    function [obj,pd_set]=pardecomp(obj,varargin)
       v=varargs.wrap('sources',{...
        {...
          %these arguments change the defaults values in common_ops_args
          'datafile',   '', @ischar;...
          'force',                        false ,@islogical;...
        }...
      },varargin{:});
      %check if data file is availble
      if ~file.exist(v.datafile) || v.force
        %get the pd-set
        out.pd_set=pardecomp.split(obj,varargin{:});
        %reconstruct the time series
        out.obj   =pardecomp.join(out.pd_set,varargin{:});
        %save the mat data
        file.save_mat(out,v.datafile)
      else
        %load the mat data
        [out,loaded_flag]=file.load_mat(v.datafile);
        assert(loaded_flag,['Problem loading ',v.datafile])
      end
      %unwrap data
      obj=out.obj;
      pd_set=out.pd_set;
    end
    function [obj,pd_set]=pardecomp_search(obj,varargin)
      %need some stuf
      v=varargs.wrap('sources',{....
        {...
          'T'                  [], @isnumeric   ;...
          'timescale',obj.x_units, @ischar      ;...
          'max_iter',          20, @num.isscalar;...
          'neg_delta',       1e-3, @num.isscalar;...
          'pos_delta',       5e-1, @num.isscalar;...
        },...
      },varargin{:});
      %init loop
      c=0;norm_now=inf;tab_len=20;
      disp(str.tablify(tab_len,'iter','norm prev','norm now','delta','new period'));
      %search
      while c<v.max_iter
        c=c+1;
        %save previous norm
        norm_prev=norm_now;
        %get pd-set
        [~,pd_set]=obj.pardecomp(v.varargin{:},'quiet',true,'T',v.T);
        %get PSD of residuals
        rf=simplefreqseries.transmute(pd_set.res).psd;
        %get period where PSD is maximum
        max_idx=find(rf.y(:,1)==max(rf.y(:,1)),1,'first');
        simpletimeseries.time2num(obj.step,0,obj.x_units);
        T_now=simpletimeseries.num2time(1/rf.x(max_idx),0,pd_set.timescale);
        %update norm and delta
        norm_now=norm(pd_set.norm.y);
        delta=(norm_now-norm_prev)/norm_now;
        %inform user
        disp(str.tablify(tab_len,c,norm_prev,norm_now,(norm_now-norm_prev)/norm_now,time.num2duration(T_now,v.timescale)))
        %do not add already existing periods
        if any(v.T==T_now)
          disp('The new period is already in T, stop searching.')
          break
        end
        if delta<0 && -delta<v.neg_delta
          disp('The norm of the residuals decreased too little, stop searching.')
          break
        end
        if delta>0 && delta>v.pos_delta
          disp('The norm of the residuals increased too much, stop searching.')
          break
        end
        %increment T
        v.T(end+1)=simpletimeseries.time2num(T_now,0,pd_set.timescale);
      end
      obj=pardecomp.join(pd_set,'time',obj.t);
    end
    %% plot methots
    function out=plot(obj,varargin)
      %call superclass
      out=plot@simpledata(obj,varargin{:});
      %annotate
      out.xlabel='time';
      xlabel(out.xlabel)
      %outputs
      if nargout == 0
        clear out
      end
    end
    %% export methods
    function out=pluck(obj,t_now)
      assert(isscalar(t_now),'Input T_now must be a scalar.')
      if isdatetime(t_now)
        x_now=obj.t2x(t_now);
      else
        x_now=t_now;
      end
      %call mother routine
      out=pluck@simpledata(obj,x_now);
    end
    function export(obj,filename,filetype,varargin)
      %NOTICE: GRACE-specific formats have been moved to grace.m
      p=machinery.inputParser;
      p.addRequired( 'filename',             @ischar);
      p.addRequired( 'filetype',             @ischar);
      v=varargs.wrap('parser',p,'sources',{{...
        'header',  'default',   @ischar;...
        'columns', 1:obj.width, @isnumeric;...
        'sat_name','',          @ischar;...
        'force',   false,       @islogical;...
      }},'mandatory',{filename,filetype},varargin{:});
      if ~exist(filename,'file') || v.force
        disp([datestr(now),': start exporting ',filename])
        %make sure this directory exists
        assert(file.ensuredir(filename),['Error creating directory of file ',filename,'.'])
        %open the file (sanity done inside)
        fid=file.open(filename,'w');
        %translate legacy usage
        if isempty(v.header); v.header='default';end
        %branch on type of file
        switch filetype
        case 'ascii'
          %enforce requested header type/value
          switch lower(v.header)
          case 'default'
            dh=[...
'# Column 1:    Date (yyyy-mm-dd)',10,...
'# Column 2:    Time (hh:mm:ss.sss)',10,...
'# Column 3:    Time system (',obj.timesystem,')',10,...
'# Column 4:    Modified Julian Day (including fraction of day)',10];
            %build rest of the default header
            for i=1:numel(v.columns)
              dh=[dh,...
                '# Column ',num2str(i+4),':    ',...
                  obj.labels{v.columns(i)},' (',...
                  obj.units{v.columns(i)},')',10];  %#ok<AGROW>
            end
            fprintf(fid,'%s',dh);
          case 'none'
            %do nothing
          otherwise
            fprintf(fid,'%s',v.header);
          end
          %build time vectors
          time_str=datestr(obj.t_masked,'yyyy-mm-dd HH:MM:SS.FFF');
          mjd=obj.mjd(obj.mask);
          %build format string
          fmt=['%s UTC %14.8f',repmat(' %16.8e',1,numel(v.columns)),'\n'];
          %build output data
          y=obj.y_masked([],v.columns);
          %sanity
          if size(time_str,1)~=size(y,1)
            error('BUG TRAP: discrepancy in the sizes of time_str and y. Debug needed.')
          end
          %save the data
          s.msg=['exporting ',obj.descriptor];s.n=size(time_str,1);
          for i=1:size(time_str,1)
            fprintf(fid,fmt,time_str(i,:),mjd(i),y(i,:));
            s=time.progress(s,i);
          end
        otherwise
          error(['Cannot handle exporting time series to files of type ''',filetype,'''.'])
        end
        fclose(fid);
      end
    end
  end
end