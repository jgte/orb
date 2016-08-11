classdef simpletimeseries < simpledata
  %static
  properties(Constant,GetAccess=private)
    valid_formats=struct(...
      'char',{{...
        'yyyy-MM-dd hh:mm:ss.sss',...
        'yyyyMMdd''T''hhmmss',...
        'yyyyMMddhhmmss.sss',...
        'yyyy MM dd hh mm ss.sss',...
        'yyyy-MM-dd',...
        'yyyyMMdd'...
      }},...
      'double',{{...
        'datenum',...
        'excel',...
        'excel1904',...
        'juliandate',...
        'modifiedjuliandate',...
        'posixtime',...
        'yyyymmdd',...
        'gpstime',...
        'gpsweeksecond',...
        'doy'...
      }},...
      'datetime',{{...
        'datetime'...
      }}...
    );
    parameter_list=struct(...
      'format',struct('default','modifiedjuliandate','validation',@(i) ischar(i)),...
      'units', struct('default',{{''}},              'validation',@(i) iscellstr(i)),...
      'debug', struct('default',false,               'validation',@(i) islogical(i) && iscalar(i))...
    );
    gps_zero_epoch='1980-01-06';
  end
  properties(Constant)
    % table of leap seconds since 6 Jan 1980:
    leap_seconds=[...
      datetime('1981-07-01'),... 1981  Jul.   1  - 1s
      datetime('1982-07-01'),... 1982  Jul.   1  - 1s
      datetime('1983-07-01'),... 1983  Jul.   1  - 1s
      datetime('1985-07-01'),... 1985  Jul.   1  - 1s
      datetime('1988-01-01'),... 1988  Jan.   1  - 1s
      datetime('1990-01-01'),... 1990  Jan.   1  - 1s
      datetime('1991-01-01'),... 1991  Jan.   1  - 1s
      datetime('1992-07-01'),... 1992  Jul.   1  - 1s
      datetime('1993-07-01'),... 1993  Jul.   1  - 1s
      datetime('1994-07-01'),... 1994  Jul.   1  - 1s
      datetime('1996-01-01'),... 1996  Jan.   1  - 1s
      datetime('1997-07-01'),... 1997  Jul.   1  - 1s
      datetime('1999-01-01'),... 1999  Jan.   1  - 1s
      datetime('2006-01-01'),... 2006  Jan.   1  - 1s
      datetime('2009-01-01'),... 2009  Jan.   1  - 1s
      datetime('2012-07-01'),... 2012  Jul.   1  - 1s 
      datetime('2015-07-01')...  2015  Jul.   1  - 1s 
    ];
  end
  %read only
  properties(SetAccess=private)
    step
    debug
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
  end
  %These parameters should not modify the data in any way; they should
  %only describe the data or the input/output format of it.
  %NOTE: if you add something here, update simpletimeseries.parameter_list
  properties(GetAccess=public,SetAccess=public)
    format
  end
  methods(Static)
    function out=timescale(in)
      out=seconds(in);
    end
    function out=valid_t(in)
      out=isdatetime(in);
    end
    function out=valid_epoch(in)
      out=isdatetime(in) && isscalar(in);
    end
    function out=time2num(in,epoch)
      if ~exist('epoch','var') || isempty(epoch)
        epoch=in(1);
      end
      out=simpletimeseries.timescale(in-epoch);
    end
    function out=num2time(in,epoch)
      if ~exist('epoch','var') || isempty(epoch)
        error([mfilename,': need input ''epoch''.'])
      end
      out=epoch+simpletimeseries.timescale(in);
    end
    function out=parameters
      out=fieldnames(simpletimeseries.parameter_list);
    end
    function out=timestep(in,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'in',              @(i) isdatetime(i));
      p.addParameter('nsigma',    4,    @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('max_iter',  10,   @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('sigma_iter',2,    @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('sigma_crit',1e-9, @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('curr_iter', 0,    @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('disp_flag', false,@(i) islogical(i));
      % parse it
      p.parse(in,varargin{:});
      %handle singularities
      switch numel(in)
        case 0
          error([mfilename,': cannot handle empty time stamps'])
        case 1
          out=NaN;
          return
      end
      %get numeric diff of time
      tdiff=simpletimeseries.timescale(diff(in));
      %get diff of time domain without jumps
      outdiff=simpledata.rm_outliers(tdiff,varargin{:});
      %get rid of nans
      outdiff=outdiff(~isnan(outdiff));
      if std(outdiff)>p.Results.sigma_crit*mean(outdiff) && p.Results.curr_iter < p.Results.max_iter
        %reduce sigma
        nsigma_new=p.Results.nsigma/p.Results.sigma_iter;
        %send feedback
        if p.Results.disp_flag
          disp([mfilename,': failed to determine the timestep, since std(delta t) is ',num2str(std(outdiff)),...
            '. Reducing NSIGMA from ',num2str(p.Results.nsigma),' to ',num2str(nsigma_new),'.'])
        end
        %recursive call
        vararginnow=simpledata.vararginclean(varargin,{'nsigma','curr_iter','disp_flag'});
        out=simpletimeseries.timestep(in,...
          'nsigma',nsigma_new,...
          'curr_iter',p.Results.curr_iter+1,...
          'disp_flag',false,...
          vararginnow{:});
      elseif isempty(outdiff)
        %dead end, sigma was reduced too much and all data is flagged as
        %outliers: nothing to do but to give some estimated of the previous
        %sigma (rounded to micro-seconds to avoid round off errors)
        vararginnow=simpledata.vararginclean(varargin,{'nsigma'});
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
        disp([mfilename,': final timestep is ',char(out),'.'])
      end
    end
    function v=fix_interp_over_gaps_narrower_than(v)
      if ~iscell(v)
        error([mfilename,': expecting input ''v'' to be a cell array, not a ',class(v),'.'])
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
    %this function converts from many forms of date/time representations to
    %matlab's 'datetime' class.
    function [out,format_out]=ToDateTime(in,format_in,debug)
      if ~exist('debug','var')
        debug=simpletimeseries.parameter_list.debug.default;
      end
      if ~exist('format_in','var')
        if ~isfield(simpletimeseries.valid_formats, class(in))
          error([mfilename,': there is no default format for inputs of class ',class(in),'; optional input ''format'' must be given.'])
        end
        format_in='';
      end
      switch class(in)
      case 'datetime'
        out=in;
        format_out='datetime'; %This is assumed to be UTC (no exceptions!)
      case 'char'
        if isempty(format_in)
          out=NaN;
          for i=simpletimeseries.valid_formats.(class(in))
            try
              out=datetime(in,'InputFormat',i{1});
            catch
              continue
            end
            %format is determined automatically, so set it
            format_out=i{1};
            break
          end
          if ~isdatetime(out)
            error([mfilename,': can not understand time ',class(in),': ',in])
          end
        else
          out=datetime(in,'InputFormat',format_in);
          %keep format, it was not attributed automatically
          format_out=format_in;
        end
      case 'double'
        if isempty(format_in)
          %assume no format is valid
          out=NaN;
          %loop over all known formats
          for i=simpletimeseries.valid_formats.(class(in))
            try
              out=datetime(in,'ConvertFrom',i{1});
            catch
              continue
            end
            %format is determined automatically, so set it
            format_out=i{1};
            break
          end
          %catch when no format if found
          if ~isdatetime(out)
            error([mfilename,': can not understand time ',class(in),': ',num2str(in)])
          end
        else
          switch format_in
          case 'gpstime'
            out=datetime(in,'convertfrom','epochtime','epoch',simpletimeseries.gps_zero_epoch);
            for i=1:numel(simpletimeseries.leap_seconds)
              out=out-seconds(out>simpletimeseries.leap_seconds(i));
            end
          case 'datevector'
            out=datetime(in);
          case 'gpsweeksecond'
            if size(in,2)~=2
              error([mfilename,': when format is ''',format_in,''', need input to have ',num2str(cols),' columns, not ',num2str(size(in,2)),'.'])
            end
            if any(floor(in(:,1))~=in(:,1))
              error([mfilename,': when format is ''',format_in,''', the first column must only contain integers.'])
            end
            out=datetime(gps2date(in(:,1),in(:,2)));
          case 'doy'
            if size(in,2)~=3
              error([mfilename,': when format is ''',format_in,''', need input to have ',num2str(cols),' columns, not ',num2str(size(in,2)),'.'])
            end
            if any(floor(in(:,1))~=in(:,1))
              error([mfilename,': when format is ''',format_in,''', the first column must only contain integers.'])
            end
            if any(floor(in(:,2))~=in(:,2))
              error([mfilename,': when format is ''',format_in,''', the second column must only contain integers.'])
            end
            tmp=datevec(datenum(in(:,1),1,1)+in(:,2)-1);       %year, month and day
            tmp(:,4) = floor(in(:,3)/3600);                    %hours
            tmp(:,5) = floor(in(:,3)/60 - tmp(:,4)*60);        %minutes
            tmp(:,6) = in(:,3) - tmp(:,4)*3600 - tmp(:,5)*60;  %seconds
            out=datetime(tmp);
          otherwise
            out=datetime(in,'ConvertFrom',format_in);
          end
          %keep format, it was not attributed automatically
          format_out=format_in;
        end          
      otherwise
        out=datetime(in);
        format_out='default';
      end
      if ~isempty(format_in) && ~strcmp(format_out,format_in) && debug
        disp(['WARNING: format changed from ''',format_in,''' to ''',format_out,'''.'])
      end
      if ~isdatetime(out)
        error([mfilename,': output must be datetime, not ',class(out),'. Debug needed!'])
      end
    end
    %this function performs the inverse convertion as 'ToDateTime'
    function out=FromDateTime(in,format)
      if ~isdatetime(in)
        error([mfilename,': input must be datetime, not ',class(in),'. Debug needed!'])
      end
      if ~exist('format','var') || isempty(format)
        format='default';
      end
      switch format
      case {'datetime','default'};
        out=in;
      case 'datenum'
        out=datenum(in);
      case 'datevec'
        out=datevec(in);
      case 'excel'
        out=exceltime(in,'1900');
      case 'excel1904'
        out=exceltime(in,'1904');
      case {'juliandate','modifiedjuliandate'}
        out=juliandate(in,format);
      case 'posixtime'
        out=posixtime(in);
      case 'yyyymmdd'
        out=yyyymmdd(in);
      case 'gpstime'
        out=seconds(in-datetime(simpletimeseries.gps_zero_epoch));
        for i=1:numel(simpletimeseries.leap_seconds)
          out=out+(in>simpletimeseries.leap_seconds(i));
        end
      case 'gpsweeksecond'
        [gps_week, gps_sow] = date2gps(simpletimeseries.FromDateTime(in,'datevec'));
        out=[gps_week,gps_sow];
      case 'doy'
        gps_week_sow=simpletimeseries.FromDateTime(in,'gpsweeksecond');
        [date, doy] = gps2date(gps_week_sow(:,1),gps_week_sow(:,2));
        out=[date(:,1),doy,date(:,4)*3600+date(:,5)*60+date(:,6)];
      otherwise
        out=char(datetime(in,'Format',format));
      end
    end
    %this function tests the reciprocity of From/ToDateTime 
    function test_time(n,max_date,col_width)
      if ~exist('n','var') || isempty(n)
        n=100;
      end
      if ~exist('max_date','var') || isempty(max_date)
        max_date=datetime([2100,12,31,23,59,59]);
      end
      if ~exist('tab_len','var') || isempty(col_width)
        col_width=[0,24,0,10,0];
      end
      for i=simpletimeseries.valid_formats.double
        switch i{1}
        case 'yyyymmdd'
          year_list=round(rand(n,1)*year( max_date));
          month_list=ceil(rand(n,1)*month(max_date));
          day_list=ceil(rand(n,1).*eomday(year_list,month_list));
          in=year_list*10000+month_list*100+day_list;
        case 'gpsweeksecond'
          max_gpsweeksecond=simpletimeseries.FromDateTime(max_date,i{1});
          in=[round(rand(n,1)*max_gpsweeksecond(1)),rand(n,1)*max_gpsweeksecond(2)];
        case 'doy'
          max_doy=simpletimeseries.FromDateTime(max_date,i{1});
          in=[round(rand(n,1)*max_doy(1)),round(rand(n,1)*max_doy(2)),rand(n,1)*max_doy(3)];
          in(in(:,2)==0,2)=1;
          in(in(:,3)==0,3)=1;
        otherwise
          in=rand(n,1)*simpletimeseries.FromDateTime(max_date,i{1});
        end
        tic
        [out,format_here]=simpletimeseries.ToDateTime(in,i{1});
        in_check=simpletimeseries.FromDateTime(out,format_here);
        dt=toc;
        switch i{1}
          case 'datenum'
            crit=1e-9;
          case 'excel'
            crit=1e-10;
          case 'excel1904'
            crit=1e-10;
          case 'juliandate'
            crit=1e-9;
          case 'modifiedjuliandate'
            crit=1e-11;
          case 'posixtime'
            crit=1e-6;
          case 'gpstime'
            crit=1e-6;
          otherwise
            crit=0;
        end
        if any(abs(in-in_check)>crit)
          idx=find(abs(in-in_check)>crit,1,'first');
          error([mfilename,': test failed for format ',i{1},':',10,...
            i{1},'  (in): ',num2str(in(idx,:)),10,...
            i{1},' (out): ',num2str(in_check(idx,:)),10,...
            'diff: ',num2str(abs(in(idx)-in_check(idx))),10,...
            'date:',datestr(simpletimeseries.ToDateTime(in(idx,:),i{1}))...
          ])
        else
          out={'Format',i{1},'ok',num2str(round(n/dt)),'ops/sec'};
          j=2;out{j}=[out{j},repmat(' ',1,col_width(j)-length(out{j}))];
          j=4;out{j}=[repmat(' ',1,col_width(j)-length(out{j})),out{j}];
          disp(strjoin(out,' '))
        end
      end
      
      for i=simpletimeseries.valid_formats.char
        year_list=round(rand(n,1)*year( max_date));
        month_list=ceil(rand(n,1)*month(max_date));
        day_list=ceil(rand(n,1).*eomday(year_list,month_list));
        hour_list=floor(rand(n,1)*24);
        minute_list=floor(rand(n,1)*60);
        second_list=floor(rand(n,1)*60);
        in=char(datetime('now','Format',i{1})+[...
          year_list,...
          month_list,...
          day_list,...
          hour_list,...
          minute_list,...
          second_list...
        ]);
        tic
        [out,format_here]=simpletimeseries.ToDateTime(in,i{1});
        in_check=simpletimeseries.FromDateTime(out,format_here);
        dt=toc;
        if ~strcmp(in,in_check)
          error([mfilename,': test failed for format ',i{1},10,in(1,:),10,in_check(1,:),10,'debug needed'])
        else
          out={'Format',i{1},'ok',num2str(round(n/dt)),'ops/sec'};
          j=2;out{j}=[out{j},repmat(' ',1,col_width(j)-length(out{j}))];
          j=4;out{j}=[repmat(' ',1,col_width(j)-length(out{j})),out{j}];
          disp(strjoin(out,' '))
        end
      end
    end
    %general test for the current object
    function test(l,w)
      if ~exist('l','var') || isempty(l)
        l=1000;
      end
      if ~exist('w','var') || isempty(w)
        w=1;
      end
%       %test the time conversions
%       simpletimeseries.test_time(l)
      %test current object
      args=simpledata.test_parameters('args',l,w);
      now=juliandate(datetime('now'),'modifiedjuliandate');
     
      
      i=0;
      a=simpletimeseries(...
        now+[1:round(l/3),round(l*2/3):round(l*4/3)],...  % t
        simpledata.test_parameters('y',l,w),...           % y
        'mask',simpledata.test_parameters('mask',l),...
        args{:},...
        'format','modifiedjuliandate'...
      );
      i=i+1;h{i}=figure('visible','off'); a.plot('title', 'original')
      
      lines1=cell(w,1);lines1(:)={'-o'};
      lines2=cell(w,1);lines2(:)={'-x'};
      lines3=cell(w,1);lines3(:)={'-+'};
      i=i+1;h{i}=figure('visible','on');
      a.plot('line',lines1)
      a.median(10).plot('line',lines2);
      a.medfilt(10).plot('line',lines3);
      legend('origina','median','medfilt')
      title('median (operation not saved)');

      b=a.resample;
      a=a.fill;
      i=i+1;h{i}=figure('visible','off');
      a.plot('line',lines1); hold on; b.plot('title','fill','line',lines2)
      legend('fill','resample')
      
      a=a.append(...
        simpletimeseries(...
          a.stop+(round(l/3):round(4*l/3)-1),...
          simpledata.test_parameters('y',l,w),...
          'mask',simpledata.test_parameters('mask',l,w),...
          args{:}...
        )...
      );
      i=i+1;h{i}=figure('visible','off'); a.plot('title','append')
      
      a=a.trim(...
        datetime(now+round(-l/2),'convertfrom','modifiedjuliandate'),...
        datetime(now+round( l/2),'convertfrom','modifiedjuliandate')...
      );
      i=i+1;h{i}=figure('visible','off'); a.plot('title','trim')

      b=a.resample(...
        days(0.8) ...
      );
      i=i+1;h{i}=figure('visible','off'); 
      a.plot('line',lines1); hold on; b.plot('title','resampled','line',lines2)
      legend('original','resampled')
      
      a=a.extend(...
        100 ...
      ).extend(...
        -100 ...
      );
      i=i+1;h{i}=figure('visible','off'); a.plot('title','extend')
      
      a=a.slice(...
        datetime(now+round(-l/5),'convertfrom','modifiedjuliandate'),...
        datetime(now+round( l/5),'convertfrom','modifiedjuliandate')...
      );
      i=i+1;h{i}=figure('visible','off'); a.plot('title','delete')
      
      for i=numel(h):-1:1
        set(h{i},'visible','on')
      end
    end
    %% import methods
    function obj=import(filename)
      %if argument filename has a wild card, then load all those files
      if ~isempty(strfind(filename,'*'))
        p=fileparts(filename);
        file_list=dir(filename);
        for i=1:numel(file_list)
          file_now=fullfile(p,file_list(i).name);
          disp([mfilename,': reading data from file ',file_now])
          %read the data
          obj_now=simpletimeseries.import(file_now);
          %determine current day
          day_now=datetime(yyyymmdd(obj_now.t(round(obj_now.length/2))),'ConvertFrom','yyyymmdd');
          %get rid of overlaps
          obj_now=obj_now.trim(day_now,day_now+hours(24)-obj_now.step);
          %append or initialize
          if i==1
            obj=obj_now;
          else
            obj=obj.append(obj_now);
          end
        end
        return
      end
      %determine file type
      [~,~,e]=fileparts(filename);
      %branch on extension
      switch e
      case '.sigma'
        fid=fopen(filename);
        raw = textscan(fid,'%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne',1);
        fclose(fid);
        %building time domain
        t=datetime([double([raw{1:5}]),raw{6}]);
        %building data domain
        y=[raw{7:end}];
        %building object
        obj=simpletimeseries(t,y,...
          'format','datetime',...
          'y_units',{'m','m','m','s','m^2','m^2','m^2','s^2','m^2','m^2','ms','m^s','ms','ms'},...
          'labels', {'x','y','z','t','xx', 'yy', 'zz', 'tt', 'xy', 'xz', 'xt','yz', 'yt','zt'},...
          'descriptor',['kinematic orbit from file ',filename]...
         );
      otherwise
        error([mfilename,': cannot handle files of type ''',e,'''.'])
      end
    end
    %% utilities
    function out=list(start,stop,period)
      p=inputParser;
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
  end
  methods
    %% constructor
    function obj=simpletimeseries(t,y,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't'); %this can be char, double or datetime
      p.addRequired( 'y',      @(i) simpledata.valid_y(i));
      %declare parameters
      for j=1:numel(simpletimeseries.parameters)
        %shorter names
        pn=simpletimeseries.parameters{j};
        %declare parameters
        p.addParameter(pn,simpletimeseries.parameter_list.(pn).default,simpletimeseries.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(t,y,varargin{:});
      % get datetime 
      [t,f]=simpletimeseries.ToDateTime(t,p.Results.format);
      %call superclass (create empty object, assignment comes later)
      obj=obj@simpledata(simpletimeseries.time2num(t),y,...
        'epoch', t(1),...
        'x_units','time',...
        'y_units',p.Results.units,...
        varargin{:}...
      );
      %save input format
      obj.format=f;
      % save parameters
      for i=1:numel(simpletimeseries.parameters)
        %shorter names
        pn=simpletimeseries.parameters{i};
        %parameter 'units' has already been handled when calling simpledata
        %constructor, so skip it
        if strcmp(pn,'units')
          continue
        end
        if ~isscalar(p.Results.(pn))
          %vectors are always lines (easier to handle strings)
          obj.(pn)=transpose(p.Results.(pn)(:));
        else
          obj.(pn)=p.Results.(pn);
        end
      end
    end
    function obj=assign(obj,y,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'y'      ,          @(i) simpledata.valid_y(i));
      p.addParameter('t'      ,obj.t,    @(i) simpletimeseries.valid_t(i));
      p.addParameter('epoch'  ,obj.epoch,@(i) simpletimeseries.valid_epoch(i));
      % parse it
      p.parse(y,varargin{:});
      % simpler names
      x_present=isfield(p.Unmatched,'x');
      t_present=~any(strcmp(p.UsingDefaults,'t'));
      %cannot have both 't' and 'x'
      if x_present && t_present
        error([mfilename,': cannot handle both inputs ''x'' and ''t''.'])
      end
      %if 't' is not present, then pass it on to simple data
      if ~t_present
        obj=assign@simpledata(obj,y,varargin{:});
        %if there is no 'x', then this is a simple assignment of y
        if ~x_present; return; end
      end
      %if 't' is present, assign 'x'
      if t_present
        obj=assign@simpledata(obj,y,'x',simpletimeseries.time2num(p.Results.t),varargin{:});
      end
      %update epoch (needed to derive obj.t from obj.x)
      if ~isempty(p.Results.epoch)
        obj.epochi=p.Results.epoch;
      elseif t_present
        obj.epochi=p.Results.t(1);
      else 
        error([mfilename,': cannot derive epoch without either input ''epoch'' or ''t''.'])
      end
      %update local records
      obj.step=simpletimeseries.timestep(obj.t);
      %sanitize (don't pass t, since it can be deliberatly changed)
      obj.check_st
    end
    function obj=copy_metadata(obj,obj_in)
      %call superclass
      obj=copy_metadata@simpledata(obj,obj_in);
      %propagate parameters of this object
      parameters=simpletimeseries.parameters;
      for i=1:numel(parameters)
        if isprop(obj,parameters{i}) && isprop(obj_in,parameters{i})
          obj.(parameters{i})=obj_in.(parameters{i});
        end
      end
    end
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      relevant_parameters={'step','format','epoch','start','stop'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print superclass
      print@simpledata(obj,tab)
    end
    %% t methods
    function x_out=t2x(obj,t_now)
      if simpletimeseries.valid_t(t_now)
        x_out=simpletimeseries.time2num(t_now,obj.epoch);
      else
        x_out=t_now;
      end
    end
    function t_out=x2t(obj,x_now)
      if simpledata.valid_x(x_now)
        t_out=simpletimeseries.num2time(x_now,obj.epoch);
      else
        t_out=x_now;
      end
    end
    function obj=set.t(obj,t_now)
      %NOTICE: this blindly changes the time domain!
      obj=obj.assign(obj.y,'t',t_now);
    end
    function out=get.t(obj)
      if isempty(obj.x)
        out=[];
      else
        out=obj.x2t(obj.x);
      end
    end
    function out=isfilled(obj)
      out=obj.length==(obj.stop-obj.start+obj.step)/obj.step;
    end
    function out=isx1zero(obj)
      %handle empty object
      if isempty(obj.x)
        out=true;
        return 
      end
      %this function checks that if obj.x(1) is zero
      %it also ensures that obj.epoch and obj.t(1) agree with that
      test=[obj.x(1)==0,obj.start==obj.epoch];
      %sanity
      if test(1)~=test(2)
        error([mfilename,':',10,...
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
    function out=start(obj)
      out=obj.t(1);
    end
    function out=stop(obj)
      out=obj.t(obj.length);
    end
    function out=span(obj)
      out=obj.stop-obj.start;
    end
    function out=t_domain(obj,step_now)
      if ~exist('step_now','var') || isempty(step_now)
        step_now=obj.step;
      end
      out=transpose(obj.start:step_now:obj.stop);
    end
    function obj=set.t_formatted(obj,t_now)
      [obj.t,format_now]=simpletimeseries.ToDateTime(t_now,obj.format);
      if ~strcmp(format_now,format_in)
        obj.format=format_now;
      end
      %sanitize
      obj.check_st(t_now)
    end
    function out=get.t_formatted(obj)
      out=simpletimeseries.FromDateTime(obj.t,obj.format);
    end
    function out=t_masked(obj,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      end
      out=obj.t(mask);
    end
    function out=idx(obj,t_now,varargin)
      out=idx@simpledata(obj,obj.t2x(t_now),varargin{:});
    end
    %% step methods
    function out=step_num(obj)
      out=simpletimeseries.timescale(obj.step);
    end
    function out=step_get(obj)
      out=simpletimeseries.timestep(obj.t);
    end
    function obj=step_update(obj)
      obj.step=simpletimeseries.timestep(obj.t);
    end
    %% epoch methods
    function obj=set.epoch(obj,epoch)
      if ~simpletimeseries.valid_epoch(epoch)
        error([mfilename,': invalid input ''epoch''.'])
      end
      %get current time domain
      t_old=obj.t;
      %set epoch
      obj.epochi=epoch;
      %shift x
      obj=obj.x_set(simpletimeseries.time2num(t_old,epoch));
      %sanity
      if any(seconds(t_old-obj.t).^2>1e-20)
        error([mfilename,': changing epoch cause the time domain to also change.'])
      end
    end
    function out=get.epoch(obj)
      out=obj.epochi;
    end
    function obj=epoch_update(obj)
      obj.epoch=obj.t(1);
    end
    %% management methods
    function check_st(obj,t_now)
      %check consistency in the values of obj.start and obj.epoch
      obj.isx1zero;
      %check for monotonously increasing time domain
      if any(diff(obj.x)<=0)
        error([mfilename,': the time domain is not monotonously increasing.'])
      end
      if exist('t_now','var') && ~isempty(t_now)
        %check for consistency in the time domain
        if numel(obj.t) ~= numel(t_now) || any(seconds(obj.t - t_now(:)).^2>1e-18)
          error([mfilename,': the time domain is not consistent with input ''t_now''.'])
        end
      end
    end
    %% edit methods
    %the add method doesn't make sense here
    %the remove method can be called directly
    function obj=trim(obj,start,stop)      
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
    function obj=resample(obj,step_now)
      % this function is a special case of interpolation
      if ~exist('step_now','var') || isempty(step_now)
        step_now=obj.step_get;
      end
      if ~isduration(step_now)
        error([mfilename,': expecting input ''step_now'' to be duration, not ',class(step_now),'.'])
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
    %the detrend method can be called directly
    %the outlier method can be called directly
    %the medfilt method can be called directly
    function obj=median(obj,n)
      %call superclass
      obj=median@simpledata(obj,n);
      %resample (if needed, which is checked inside resample)
      obj=obj.resample;
    end
    function obj=extend(obj,nr_epochs)
      switch class(nr_epochs)
      case 'double'
        if nr_epochs==0
          return
        end
        if (nr_epochs~=round(nr_epochs))
          error([mfilename,': input ''nr_epochs'' must be an integer, not ',num2str(nr_epochs),'.'])
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
          error([mfilename,': input ''t'' (',datestr(t_now),') ',...
            'must be larger than obj.stop (',datestr(obj.stop),') ',...
            'or smaller than than obj.start (',datestr(obj.start),').'...
          ]);
        end
        obj=extend(obj,(t_now-t_ref)/obj.step);
      otherwise
        error([mfilename,': cannot handle input ''nr_epochs'' of class ',class(nr_epochs),'.'])
      end
    end
    %TODO: check if this routine is redundant
    function obj=fill(obj)
      %trivial call
      if obj.isfilled
        return
      end
      %build complete time domain
      t_new=obj.t_domain;
      t_now=obj.t;
      disp(['Need to fill in missing epochs: ',num2str(numel(t_new)-obj.length),' ('...
        num2str((numel(t_new)-obj.length)/numel(t_new)*1e2),'%).'])
      %find out where there are gaps
      gap_idx=find(diff(t_now)~=obj.step);
      %loop over all implicit gaps (i.e. missing epochs)
      s.msg=[mfilename,': populating missing epochs'];s.n=numel(gap_idx);
      for i=1:numel(gap_idx)
        %create patch
        t_patch=transpose((t_now(gap_idx(i))+obj.step):obj.step:(t_now(gap_idx(i)+1)-obj.step));
        %save data with patch (it gets deleted when assigning to x)
        y_patched=[obj.y(1:gap_idx(i),:);...
                   nan(numel(t_patch),obj.width);...
                   obj.y(gap_idx(i)+1:end,:)];
        %create patched t
        t_patched=[t_now(1:gap_idx(i));...
                  t_patch;...
                  t_now(gap_idx(i)+1:end)];
        %propagate y
        obj=obj.assign(y_patched,'t',t_patched);
        %user feedback
        s=time.progress(s,i);
      end
      %sanitize
      obj.check_st(t_new);
    end
    %% multiple object manipulation
    function out=istequal(obj1,obj2)
      %make sure things are up to date
      obj1=obj1.t_reset;
      obj2=obj2.t_reset;
      out= obj1.length == obj2.length && ...
           obj1.start  == obj2.start  && ...
           obj1.stop   == obj2.stop   && ...
           obj1.step   == obj2.step   && ...
           obj1.epoch  == obj2.epoch;
    end
    function [obj1,obj2]=consolidate(obj1,obj2)
      %extends the time domain of both objects to be in agreement
      %with the each other
      compatible(obj1,obj2)
      %trivial call
      if istequal(obj1,obj2)
        return
      end
      %build extended time domain, with gcd timesteo, rounded to the nearest second
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
    function out=step_gcd(obj1,obj2)
      out=simpletimeseries.timescale(gcd(...
        simpletimeseries.timescale(obj1.step),...
        simpletimeseries.timescale(obj2.step)...
      ));
    end
    function out=step_lcm(obj1,obj2)
      out=simpletimeseries.timescale(lcm(...
        simpletimeseries.timescale(obj1.step),...
        simpletimeseries.timescale(obj2.step)...
      ));
    end
    function [obj1,obj2]=matchstep(obj1,obj2)
      %trivial call
      if obj1.step==obj2.step
        return
      end
      %new timestep is the greatest common divisor
      step_now=step_gcd(obj1,obj2);
      %resample to the common step size
      obj1=obj1.resample(step_now);
      obj2=obj2.resample(step_now);
    end
    function [obj,idx1,idx2]=append(obj1,obj2)
      if isa(obj1,'simpletimeseries') && isa(obj2,'simpletimeseries')
        %need step to be the same
        if obj1.step~=obj2.step
          [obj1,obj2]=matchstep(obj1,obj2);
        end
        %match epochs
        obj2.epoch=obj1.epoch;
      end
      %call upstream method
      [obj,idx1,idx2]=append@simpledata(obj1,obj2);
    end
    %% plot methots
    function out=plot(obj,varargin)
      %call superclass
      out=plot@simpledata(obj,varargin{:});
      %using internal Matlab representation for dates
      lines_now=get(gca,'children');
      for i=1:numel(lines_now)
        if numel(lines_now(i).XData) == obj.length
          lines_now(i).XData=datenum(obj.t);
        end
      end
      set(gca,'XTick',datenum(obj.t));
      datetick('x',time.format(seconds(obj.span)))
      %annotate
      out.xlabel='time';
      xlabel(out.xlabel)
      %outputs
      if nargout == 0
        clear out
      end
    end
    %% export methods
    function ascii(obj,filename,varargin)
      default_header=[...
'# Column 1:    Date (yyyy-mm-dd)',10,...
'# Column 2:    Time (hh:mm:ss.sss)',10,...
'# Column 3:    Time system (UTC)',10,...
'# Column 4:    Modified Julian Day (including fraction of day)',10];
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'filename',             @(i) ischar(i));
      p.addParameter('header',  '',          @(i) ischar(i));
      p.addParameter('columns', 1:obj.width, @(i)isnumeric(i));
      % parse it
      p.parse(filename,varargin{:});
      if isempty(dir(filename))
        disp([datestr(now),': start exporting ',filename])
        %open the file
        [fid,msg]=fopen(filename,'w');
        if fid <=0
          error([mfilename,': error opening ',filename,': ',msg])
        end
        %write the header
        if isempty(p.Results.header)
          %use default header, none was specified
          header=default_header;
          %build rest of the default header
          for i=1:numel(p.Results.columns)
            header=[header,...
              '# Column ',num2str(i+4),':    ',...
                obj.labels{p.Results.columns(i)},' (',...
                obj.y_units{p.Results.columns(i)},')',10]; %#ok<AGROW>
          end
        else
          header=p.Results.header;
        end
        fprintf(fid,'%s',header);
        %build time vectors
        time_str=datestr(obj.t_idx(obj.mask),'yyyy-mm-dd HH:MM:SS.FFF');
        mjd=simpletimeseries.FromDateTime(obj.t_idx(obj.mask),'modifiedjuliandate');
        %build format string
        fmt='%s UTC %14.8f';
        for j=1:numel(p.Results.columns)
          fmt=[fmt,' %16.8e']; %#ok<AGROW>
        end
        fmt=[fmt,'\n'];
        %build output data
        y=obj.y(obj.mask,p.Results.columns);
        %sanity
        if size(time_str,1)~=size(y,1)
          error([mfilename,': discrepancy in the sizes of time_str and y. Debug needed.'])
        end
        %save the data
        s.msg=['exporting ',obj.descriptor];s.n=size(time_str,1);
        for i=1:size(time_str,1)
          fprintf(fid,fmt,time_str(i,:),mjd(i),y(i,:));
          s=time.progress(s,i);
        end
        fclose(fid);
      end
    
    end
  end
end

function [date, doy, dow] = gps2date(gps_week, gps_sow)
  % SYNTAX:
  %   [date, doy, dow] = gps2date(gps_week, gps_sow);
  %
  % INPUT:
  %   gps_week = GPS week
  %   gps_sow  = GPS seconds of week
  %
  % OUTPUT:
  %   date = date [year month day hour min sec]
  %   doy  = day of year
  %   dow  = day of week
  %
  % DESCRIPTION:
  %   Conversion from GPS time to calendar date and day of year (DOY).

  %----------------------------------------------------------------------------------------------
  %                           goGPS v0.4.3
  %
  % Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
  %----------------------------------------------------------------------------------------------
  %
  %    This program is free software: you can redistribute it and/or modify
  %    it under the terms of the GNU General Public License as published by
  %    the Free Software Foundation, either version 3 of the License, or
  %    (at your option) any later version.
  %
  %    This program is distributed in the hope that it will be useful,
  %    but WITHOUT ANY WARRANTY; without even the implied warranty of
  %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  %    GNU General Public License for more details.
  %
  %    You should have received a copy of the GNU General Public License
  %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  %----------------------------------------------------------------------------------------------

  gps_start_datenum = 723186; %This is datenum([1980,1,6,0,0,0])

  gps_dow = fix(gps_sow/86400);                             %day of week
  date = datevec(gps_start_datenum + 7*gps_week + gps_dow); %calendar date up to days
  gps_sod = gps_sow - gps_dow*86400;                        %seconds of day
  date(:,4) = floor(gps_sod/3600);                          %hours
  date(:,5) = floor(gps_sod/60 - date(:,4)*60);             %minutes
  date(:,6) = gps_sod - date(:,4)*3600 - date(:,5)*60;      %seconds

  %day of year (DOY)
  if (nargout > 1)
      doy = date2doy(datenum(date));
      doy = floor(doy);
  end

  %day of week
  if (nargout > 2)
      dow = gps_dow;
  end
end

function [gps_week, gps_sow, gps_dow] = date2gps(date)
  % SYNTAX:
  %   [gps_week, gps_sow, gps_dow] = date2gps(date);
  %
  % INPUT:
  %   date = date [year, month, day, hour, min, sec]
  %
  % OUTPUT:
  %   gps_week = GPS week
  %   gps_sow  = GPS seconds of week
  %   gps_dow  = GPS day of week
  %
  % DESCRIPTION:
  %   Conversion from calendar date to GPS time.

  %----------------------------------------------------------------------------------------------
  %                           goGPS v0.4.3
  %
  % Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
  %----------------------------------------------------------------------------------------------
  %
  %    This program is free software: you can redistribute it and/or modify
  %    it under the terms of the GNU General Public License as published by
  %    the Free Software Foundation, either version 3 of the License, or
  %    (at your option) any later version.
  %
  %    This program is distributed in the hope that it will be useful,
  %    but WITHOUT ANY WARRANTY; without even the implied warranty of
  %    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  %    GNU General Public License for more details.
  %
  %    You should have received a copy of the GNU General Public License
  %    along with this program.  If not, see <http://www.gnu.org/licenses/>.
  %----------------------------------------------------------------------------------------------

  gps_start_datenum = 723186; %This is datenum([1980,1,6,0,0,0])

  %number of days since the beginning of GPS time
  deltat   = (datenum([date(:,1), date(:,2), date(:,3)]) - gps_start_datenum);

  gps_week = floor(deltat/7);            %GPS week
  gps_dow  = floor(deltat - gps_week*7); %GPS day of week
  gps_sow  = (deltat - gps_week*7)*86400; 
  gps_sow = gps_sow + date(:,4)*3600 + date(:,5)*60 + date(:,6); %GPS seconds of week
end
  
function [doy,fraction] = date2doy(inputDate)
%DATE2DOY  Converts the date to decimal day of the year.
  %   [doy,fraction] = date2doy(inputDate)
  %
  %   Descriptions of Input Variables:
  %   inputDate:  Input date as a MATLAB serial datenumber
  %
  %   Descriptions of Output Variables:
  %   doy: Decimal day of the year. For instance, an output of 1.5 would 
  %       indicate that the time is noon on January 1.
  %   fraction: Outputs the fraction of a year that has elapsed by the input
  %       date.
  %
  %   Example(s):
  %   >> [doy,frac] = date2doy(datenum('1/25/2004'));
  %
  %   See also:

  % Author: Anthony Kendall
  % Contact: anthony [dot] kendall [at] gmail [dot] com
  % Created: 2008-03-11
  % Copyright 2008 Michigan State University.

  %Want inputs in rowwise format
  [doy,fraction] = deal(zeros(size(inputDate)));
  inputDate = inputDate(:);

  %Parse the inputDate
  [dateVector] = datevec(inputDate);

  %Set everything in the date vector to 0 except for the year
  dateVector(:,2:end) = 0;
  dateYearBegin = datenum(dateVector);

  %Calculate the day of the year
  doyRow = inputDate - dateYearBegin;

  %Optionally, calculate the fraction of the year that has elapsed
  flagFrac = (nargout==2);
  if flagFrac
      %Set the date as the end of the year
      dateVector(:,1) = dateVector(:,1) + 1;
      dateYearEnd = datenum(dateVector);
      fracRow = (doyRow - 1) ./ (dateYearEnd - dateYearBegin);
  end

  %Fill appropriately-sized output array
  doy(:) = doyRow;
  if flagFrac
      fraction(:) = fracRow;
  end
end