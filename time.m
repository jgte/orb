classdef time
  properties(Constant,GetAccess=private)
    formt_list={'S.FFF' ,'SS.FFF','HH:MM:SS','HH:MM','dd/mm HH:MM','dd/mm/yy','mm/yyyy' ,'yyyy'}
    units_list={'ms'    ,'s'     ,'min'     ,'hrs'  ,'days'       ,'mon'     ,'yrs'     ,'cent'};
    magnt_list=[0.01    ,60      ,3600      ,86400  ,2678400      ,31536000  ,3153600000,inf   ];
    funct_list={...
      @(i) seconds(i*1e-3),...
      @(i) seconds(i),...
      @(i) minutes(i),...
      @(i) hours(i),...
      @(i) days(i),...
      @(i) days(365.25*i/12),...
      @(i) years(i),...
      @(i) years(i*100)...
    };
    millennium=2000;
    valid_formats=struct(...
      'char',{{...
        'yyyy-MM-dd HH:mm:ss.sss',...
        'yyyyMMdd''T''HHmmss',...
        'yyyyMMddHHmmss.sss',...
        'yyyy MM dd HH mm ss.sss',...
        'yyyy-MM-dd',...
        'yyyyMMdd'...
      }},...                     %needs revision when adding new numeric format
      'double',{{...             https://www.mathworks.com/help/matlab/ref/datetime.html#namevaluepairarguments
        'datenum'           ,... Number of days since 0-Jan-0000 (proleptic ISO calendar).
        'excel'             ,... Number of days since 0-Jan-1900.
        'excel1904'         ,... Number of days since 0-Jan-1904.
        'juliandate'        ,... Number of days since noon UTC 24-Nov-4714 BCE (proleptic Gregorian calendar).
        'modifiedjuliandate',... Number of days since midnight UTC 17-Nov-1858.
        'posixtime'         ,... Number of seconds since 1-Jan-1970 00:00:00 UTC, not counting leap seconds.
        'yyyymmdd'          ,... Dates as YYYYMMDD numeric values. For example, 20140402 represents April 2, 2014.
        'J2000sec'          ,... Number of seconds since 1-Jan-2000 00:00:00 UTC, not counting leap seconds.
        'Week0sec'          ,... Number of seconds since 6-Jan-1980 00:00:00 UTC, not counting leap seconds.
        'gpstime'           ,...
        'gpsweeksecond'     ,...
        'yeardoysec'        ...
      }},...
      'uint64',{{...
        'ntp'               ,... Number of "clock ticks" since 1-Jan-1900 00:00:00 UTC, where each clock tick is 2^-32 seconds.
        '.net'              ,... Number of "clock ticks" since 1-Jan-0001 00:00:00 UTC, where each clock tick is 100 ns.
        'ntfs'              ,... Number of "clock ticks" since 1-Jan-1601 00:00:00 UTC, where each clock tick is 100 ns.
      }},...
      'datetime',{{...
        'datetime'...
      }}...
    );
    gps_zero_epoch='1980-01-06';
  end
  properties(Constant)
    zero_date=datetime(0,  'ConvertFrom','datenum');
     inf_date=datetime(Inf,'ConvertFrom','datenum');
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
      datetime('2017-01-01')...  2017  Jan.   1  - 1s
    ];
  end
  methods(Static)
    function test
      disp(' - time.{ToDateTime,FromDateTime} -')
      d0=datetime('now');
      fmt=[time.valid_formats.double(:);time.valid_formats.uint64(:);time.valid_formats.char(:)];
      for i=1:numel(fmt)
        d1=time.FromDateTime(d0,fmt{i});
        d2=time.ToDateTime(d1,fmt{i});
        switch fmt{i}
        case {'yyyymmdd','yyyy-MM-dd','yyyyMMdd'}
          dtest=dateshift(d0, 'start', 'day');
        otherwise
          dtest=d0;
        end
        str.say(fmt{i},':',d1,' +/- ',seconds(abs(dtest-d2)),' sec')
      end
      disp(' - time.{scale,units,format,str} -')
      fmt='%012.12g';
      for i=-4:11
        x=10^i;
        disp([...
          'magnitd: ',num2str(i,'%02i'),'; ',...
          'seconds: ',num2str(x,fmt),'; ',...
          'scale  : ',num2str(time.scale(x)),'; ',...
          'units  : ',time.units(x),'; ',...
          'fmt    : ',time.format(x),': ',...
          'str    : ',time.str(x),...
        ])
      end
    end
    %% tests for mutated inputs
    function out=mutate(in)
      switch class(in)
      case 'cell'
        out=cellfun(@time.mutate,in);
      case 'datetime'
        out=in;
      case {'double','single','uint8','uint16','uint32','uint64','int8','int16','int32','int64'}
        out=datenum(in);
      otherwise
        error(['Cannot handle class ''',class(in),'''.'])
      end
    end
    function out=iszero(in)
      switch class(in)
      case 'cell'
        out=cellfun(@time.iszero,in);
      otherwise
        out=seconds(time.mutate(in)-time.zero_date)<2e-6;
      end
    end
    function out=isinf(in)
      switch class(in)
      case 'cell'
        out=cellfun(@time.isinf,in);
      otherwise
        out=any(isinf(time.mutate(in)));
      end
    end
    function out=isfinite(in)
      if isempty(in)
        out=false;
      else
        switch class(in)
        case 'cell'
          out=cellfun(@time.isfinite,in);
        otherwise
          out=~time.iszero(in)&&~time.isinf(in);
        end
      end
    end
    %% constructors
    function out=rand(n)
      out=datetime(sort(datenum(round([...
        rand(n,1),...
        rand(n,1)*12,...
        rand(n,1)*31,...
        rand(n,3)*60 ...
      ]))),'ConvertFrom','datenum')+years(2000);
    end
    function out=domain(start,stop,step)
      out=transpose(start:step:stop);
      %make sure out includes obj.stop (hopefully patched with NaNs later)
      if out(end)<stop
        out(end+1)=out(end)+step;
      end
    end
    %% auto time/duration/strings
    function out=scale(in)
      out=find(in<=time.magnt_list,1,'first');
    end
    function out=scale_from_units(units)
      units=time.translate_units(units);
      for i=1:numel(time.units_list)
        if strcmp(units,time.units_list{i})
          out=time.magnt_list(i);
          return
        end
      end
    end
    function out=translate_units(in)
      switch lower(in)
      case {'milisecond','miliseconds','ms'}
        out = 'ms';
      case {'second','seconds','sec','s'}
        out = 's';
      case {'minute','minutes','min','m'}
        out = 'min';
      case {'hour','hours','hr','hrs','h'}
        out = 'hrs';
      case {'day','days','d'}
        out = 'days';
      case {'week','weeks','w'}
        out = 'weeks';
      case {'month','months','mon'}
        out = 'mon';
      case {'year','years','yrs','y'}
        out = 'yrs';
      case {'centuries','century','cent','c'}
        out='cent';
      otherwise
        warning(['unknown time units ''',in,'''. Using seconds.'])
        out = 's';
      end
    end
    function out=units(in)
      out=time.units_list{time.scale(in)};
    end
    function out=format(in)
      out=time.formt_list{time.scale(in)};
    end
    function [out,units]=scale_domain(in,units)
      % [OUT,UNITS,SCALE]=SCALE_DOMAIN(IN,UNITS) gets the units from the
      % input time vector and returns a new scaled time vector appropriate
      % for plotting.
      if ~exist('units','var') || isempty(units)
        units=time.units(in);
      end
      out=in./time.scale_from_units(units);
    end
    function out = str(seconds,element_max)
      if ~exist('element_max','var') || isempty(element_max)
        element_max=3;
      end
      if ~isscalar(seconds)
        error('input must be scalar')
      end
      if isnan(seconds)
        out='NaN';
      elseif ~isfinite(seconds)
        out='Inf';
      else
        steps = [1,60,60*60,60*60*24,60*60*24*7,60*60*24*31,60*60*24*365];
        units = ['s','m','h','d','w','m','y'];
        out=cell(size(units));
        seconds = abs(seconds);
        for i=numel(steps):-1:1
          if (seconds > steps(i))
            if i~=1
              out{i} = [num2str(floor(seconds/steps(i))),units(i)];
            else
              out{i} = [num2str(seconds,3),units(i)];
            end
            seconds = rem(seconds,steps(i));
          elseif (i==1)
            out{i} = [num2str(seconds,3),units(i)];
          end
        end
        %remove empty
        out=cells.rm_empty(out);
        %enforce max elements
        out=out(numel(out):-1:numel(out)-min([element_max,numel(out)])+1);
        %make a string out of it
        out=strjoin(out,' ');
      end
    end
    function out=num2duration(in,units)
      units=time.translate_units(units);
      for i=1:numel(time.units_list)
        if strcmp(units,time.units_list{i})
          out=time.funct_list{i}(in);
          return
        end
      end
      error(['Cannot understand units ''',units,'''; options are: ',strjoin(time.units_list,','),'.'])
    end
    function out=duration2num(in,units)
      out=time.num2duration(in,units);
    end
    %% monitor lengthy iterations
    function s=progress(s,i)
      % %Example:
      % s.msg='something'; s.n=n;
      % for i=1:n
      %     <do something>
      %     s=time.progress(s,i);
      % end
      % clear s
      % %or:
      % s.msg='something'; s.n=n;
      % while criteria
      %     <do something>
      %     criteria=<some check>
      %     s=time.progress(s);
      % end
      if ~exist('i','var') || isempty(i)
        if isfield(s,'c')
          i=s.c+1;
        else
          i=1;
        end
      end
      if ~exist('s','var') || isempty(s) || ~isfield(s,'n')
        error('need structure in argument ''s'' with field ''n'' (the loop total nr of iterations).');
      end
      %defaults
      if ~isfield(s,'step')
        s.step=1; %seconds
      end
      if ~isfield(s,'msg')
        s.msg='progress';
      end
      %initialize
      if ~isfield(s,'l')
        s.l=0;
        s.l_max=0;
        s.e_show=0;
        s.dt=zeros(1,s.n);
        s.c=0;
      end
      if ~isfield(s,'stationary')
        s.stationary=true;
      end
      %finalize
      if i>=s.n
        if s.c>0
          fprintf('%s',repmat(8,1,s.l_max))
          out=sprintf('%s',[s.msg,': ',num2str(i,'%03d'),' iter done in ',time.str(sum(s.dt))]);
          fprintf('%s\n',[out,repmat(' ',1,s.l_max-length(out))])
        end
        return
      end
      %gather measurements
      if s.c>0
        dt=toc;
      else
        dt=0;
      end
      s.c=s.c+1;
      s.dt(s.c)=dt;
      s.dt_std=std(s.dt(1:s.c));
      s.e=sum(s.dt);
      s.f=(s.n-i)*s.e/s.c;
      s.f_err=s.dt_std*(s.n-i);
      %give feedback during iterations
      if (s.e-s.e_show > s.step)
        fprintf('%s',repmat(8,1,s.l_max))
        out=sprintf('%s',[s.msg,': ',num2str(round(i/s.n*100)),'% (',...
          'finished in ',time.str( s.f ),' +/-',time.str(s.f_err),...
          ', elapsed ',  time.str( s.e ),...
          ', iteration ',time.str( s.e/s.c ),' +/-',time.str(s.dt_std),...
          ', total ',    time.str( s.e+s.f-s.f_err ),' - ',time.str( s.e+s.f+s.f_err )...
          ')']);
        if s.stationary
          fprintf('%s',[out,repmat(' ',1,s.l_max-length(out))])
        else
          fprintf('%s\n',out)
        end
        s.l=length(out);      %nr of chars to be deleted
        s.l_max=max([s.l_max,s.l]);
        s.e_show=s.e;   %next user feedback only after s.step has elapsed after this
      end
      tic
    end
    %% lists
    function [startlist,stoplist]=day_list(start,stop)
      p=machinery.inputParser;
      p.addRequired( 'start',   @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'stop',    @(i) isscalar(i) && isdatetime(i));
      p.parse(start,stop)
      round_start=dateshift(start,'start','day');
        round_end=dateshift(stop, 'start','day');
      startlist=round_start:days(1):round_end;
      if nargout>1
        stoplist=startlist+days(1);
      end
    end
    function [startlist,stoplist]=month_list(start,stop)
      p=machinery.inputParser;
      p.addRequired( 'start',   @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'stop',    @(i) isscalar(i) && isdatetime(i));
      p.parse(start,stop)
      %https://www.mathworks.com/matlabcentral/answers/73749-how-do-i-create-a-vector-with-the-first-day-of-each-month#answer_83695
      dvec      = datevec(datenum(start):datenum(stop));
      duniq     = unique(dvec(:, 1:2), 'rows');
      startlist = datetime(datenum(duniq(:,1), duniq(:,2), 1),'ConvertFrom','datenum');
      if nargout>1
        if numel(startlist)==1
          stoplist = dateshift(startlist,'end','month')+days(1);
        else
          stoplist = [startlist(2:end);dateshift(startlist(end),'end','month')+days(1)];
        end
      end
    end
    function [startlist,stoplist]=year_list(start,stop)
      p=machinery.inputParser;
      p.addRequired( 'start',   @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'stop',    @(i) isscalar(i) && isdatetime(i));
      p.parse(start,stop)
      %https://www.mathworks.com/matlabcentral/answers/73749-how-do-i-create-a-vector-with-the-first-day-of-each-month#answer_83695
      dvec      = datevec(datenum(start):datenum(stop));
      duniq     = unique(dvec(:, 1), 'rows');
      startlist = datetime(datenum(duniq(:,1), 1, 1),'ConvertFrom','datenum');
      if nargout>1
        if numel(startlist)>1
          stoplist = [startlist(2:end);dateshift(startlist(end),'end','year')+days(1)];
        else
          stoplist = dateshift(startlist(end),'end','year')+days(1);
        end
      end
    end
    function [startlist,stoplist]=list(start,stop,period)
      p=machinery.inputParser;
      p.addRequired( 'start',   @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'stop',    @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'period',  @isduration);
      p.parse(start,stop,period)
      switch period
      case days(1)
        [startlist,stoplist]=time.day_list(start,stop);
      case years(1)/12
        [startlist,stoplist]=time.month_list(start,stop);
      case years(1)
        [startlist,stoplist]=time.year_list(start,stop);
      otherwise
        startlist=start:period:stop;
        stoplist=startlist+period;
      end
    end
    %% time representations
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
          doy = time.date2doy(datenum(date));
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
    function date=doy2datetime(year,doy)
      date=datetime(datenum([year 0 0 0 0 0])+doy,'ConvertFrom','datenum');
    end
    function out=gpssec2datetime(in,zero_epoch)
      if ~exist('zero_epoch','var')
        zero_epoch=time.gps_zero_epoch;
      end
      out=datetime(in,...
        'convertfrom','epochtime',...
        'epoch',zero_epoch...
      );
    end
    function out=datetime2gpssec(in,zero_epoch)
      if ~exist('zero_epoch','var')
        zero_epoch=time.gps_zero_epoch;
      end
      out=seconds(in-datetime(zero_epoch));
    end
    function utc=gps2utc(gps)
      utc=gps;
      for i=1:numel(time.leap_seconds)
        utc=utc-seconds(utc>time.leap_seconds(i));
      end
    end
    function gps=utc2gps(utc)
      gps=utc;
      for i=1:numel(time.leap_seconds)
        gps=gps+seconds(gps>time.leap_seconds(i));
      end
    end
    %this function converts from many forms of date/time representations to
    %matlab's 'datetime' class.
    function [out,format_out]=ToDateTime(in,format_in,debug)
      if ~exist('debug','var')
        debug=false;
      end
      if ~exist('format_in','var')
        assert(isfield(time.valid_formats, class(in)),...
          ['there is no default format for inputs of class ',class(in),...
          '; optional input ''format'' must be given.'])
        format_in='';
      end
      switch class(in)
      case 'datetime'
        out=in;
        format_out='datetime'; %This is assumed to be UTC (no exceptions!)
      case 'char'
        if isempty(format_in)
          out=NaN;
          for i=time.valid_formats.(class(in))
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
            error(['can not understand time ',class(in),': ',in])
          end
        else
          out=datetime(in,'InputFormat',format_in);
          %keep format, it was not attributed automatically
          format_out=format_in;
        end
      case {'double','uint64'}
        if isempty(format_in)
          %assume no format is valid
          out=NaN;
          %loop over all known formats
          for i=time.valid_formats.(class(in))
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
            error(['can not understand time ',class(in),': ',num2str(in)])
          end
        else
          switch format_in %needs revision when adding new numeric format
          case 'gpstime'
            out=time.gps2utc(...
              time.gpssec2datetime(in)...
            );
          case 'J2000sec'
            out=datetime(in,...
              'convertfrom','epochtime',...
              'epoch','2000-01-01'...
            );
          case 'Week0sec'
            out=datetime(in,...
              'convertfrom','epochtime',...
              'epoch','1980-01-06'...
            );
          case 'datevector'
            out=datetime(in);
          case 'gpsweeksecond'
            cols=2;
            if size(in,2)~=cols
              error(['when format is ''',format_in,''', need input to have ',num2str(cols),' columns, not ',num2str(size(in,2)),'.'])
            end
            if any(floor(in(:,1))~=in(:,1))
              error(['when format is ''',format_in,''', the first column must only contain integers.'])
            end
            out=datetime(time.gps2date(in(:,1),in(:,2)));
          case 'yeardoysec'
            cols=3;
            if size(in,2)~=cols
              error(['when format is ''',format_in,''', need input to have ',num2str(cols),' columns, not ',num2str(size(in,2)),'.'])
            end
            if any(floor(in(:,1))~=in(:,1))
              error(['when format is ''',format_in,''', the first column must only contain integers.'])
            end
            if any(floor(in(:,2))~=in(:,2))
              error(['when format is ''',format_in,''', the second column must only contain integers.'])
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
        error(['output must be datetime, not ',class(out),'. Debug needed!'])
      end
    end
    %this function performs the inverse convertion as 'ToDateTime'
    function out=FromDateTime(in,format)
      if ~isdatetime(in)
        error(['input must be datetime, not ',class(in),'. Debug needed!'])
      end
      if ~exist('format','var') || isempty(format)
        format='default';
      end
      switch format %needs revision when adding new numeric format
      case {'datetime','default'}
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
      case 'ntp'
        out= uint64((in-datetime(1900,1,1,0,0,0))/seconds(2^-32));
      case '.net'
        out= uint64((in-datetime(0001,1,1,0,0,0))/seconds(100e-9));
      case 'ntfs'
        out= uint64((in-datetime(1601,1,1,0,0,0))/seconds(100e-9));
      case 'J2000sec'
        out=seconds(in-datetime(2000,1,1,0,0,0));
      case 'Week0sec'
        out=seconds(in-datetime(1980,1,6,0,0,0));
      case 'gpstime'
        out=time.datetime2gpssec(...
          time.utc2gps(in)...
        );
      case 'gpsweeksecond'
        [gps_week, gps_sow] = time.date2gps(time.FromDateTime(in,'datevec'));
        out=[gps_week,gps_sow];
      case 'yeardoysec'
        gps_week_sow=time.FromDateTime(in,'gpsweeksecond');
        [date, doy] = time.gps2date(gps_week_sow(:,1),gps_week_sow(:,2));
        out=[date(:,1),doy,date(:,4)*3600+date(:,5)*60+date(:,6)];
      otherwise
        out=char(datetime(in,'Format',format));
      end
    end
    function out=sod(in)
      % Converts the date to second of day.
      out=seconds(in-dateshift(in,'start','day'));
    end
    function out=fod(in)
      % Converts the date to fraction of day.
      out=time.sod(in)/seconds(days(1));
    end
    function out=doy(in)
      % Converts the date to second of day.
      out=days(in-dateshift(in,'start','year')+1);
    end
    %shortcuts
    function out=mjd(in)
      out=time.FromDateTime(in,'modifiedjuliandate');
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
      for i=time.valid_formats.double
        switch i{1}
        case 'yyyymmdd'
          year_list=round(rand(n,1)*year( max_date));
          month_list=ceil(rand(n,1)*month(max_date));
          day_list=ceil(rand(n,1).*eomday(year_list,month_list));
          in=year_list*10000+month_list*100+day_list;
        case 'gpsweeksecond'
          max_gpsweeksecond=time.FromDateTime(max_date,i{1});
          in=[round(rand(n,1)*max_gpsweeksecond(1)),rand(n,1)*max_gpsweeksecond(2)];
        case 'yeardoysec'
          max_doy=time.FromDateTime(max_date,i{1});
          in=[round(rand(n,1)*max_doy(1)),round(rand(n,1)*max_doy(2)),rand(n,1)*max_doy(3)];
          in(in(:,2)==0,2)=1;
          in(in(:,3)==0,3)=1;
        otherwise
          in=rand(n,1)*time.FromDateTime(max_date,i{1});
        end
        tic
        [out,format_here]=time.ToDateTime(in,i{1});
        in_check=time.FromDateTime(out,format_here);
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
          error(['test failed for format ',i{1},':',10,...
            i{1},'  (in): ',num2str(in(idx,:)),10,...
            i{1},' (out): ',num2str(in_check(idx,:)),10,...
            'diff: ',num2str(abs(in(idx)-in_check(idx))),10,...
            'date:',datestr(time.ToDateTime(in(idx,:),i{1}))...
          ])
        else
          out={'Format',i{1},'ok',num2str(round(n/dt)),'ops/sec'};
          j=2;out{j}=[out{j},repmat(' ',1,col_width(j)-length(out{j}))];
          j=4;out{j}=[repmat(' ',1,col_width(j)-length(out{j})),out{j}];
          disp(strjoin(out,' '))
        end
      end

      for i=time.valid_formats.char
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
        [out,format_here]=time.ToDateTime(in,i{1});
        in_check=time.FromDateTime(out,format_here);
        dt=toc;
        if ~strcmp(in,in_check)
          error(['test failed for format ',i{1},10,in(1,:),10,in_check(1,:),10,'debug needed'])
        else
          out={'Format',i{1},'ok',num2str(round(n/dt)),'ops/sec'};
          j=2;out{j}=[out{j},repmat(' ',1,col_width(j)-length(out{j}))];
          j=4;out{j}=[repmat(' ',1,col_width(j)-length(out{j})),out{j}];
          disp(strjoin(out,' '))
        end
      end
    end
    function out=round_seconds(in)
      out=seconds(round(seconds(in)));
    end
    function out=month2int(in)
      switch lower(in)
      case {'january',  'jan','1','01'};out=1;
      case {'february', 'feb','2','02'};out=2;
      case {'march',    'mar','3','03'};out=3;
      case {'april',    'apr','4','04'};out=4;
      case {'may',            '5','05'};out=5;
      case {'june',     'jun','6','06'};out=6;
      case {'july',     'jul','7','07'};out=7;
      case {'august',   'aug','8','08'};out=8;
      case {'september','sep','9','09'};out=9;
      case {'october',  'oct',    '10'};out=10;
      case {'november', 'nov',    '11'};out=11;
      case {'december', 'dec',    '12'};out=12;
      otherwise
        error(['Cannot understand month ''',in,'''.'])
      end
    end
    %% millennium fixing routines
    function year=millennium_add(year)
      assert(isnumeric(year),'Input ''year'' has to be numeric.')
      fix_idx=year<time.millennium;
      if any(fix_idx)
        year(fix_idx)=year(fix_idx)+time.millennium;
      end
    end
    function year=millennium_rm(year)
      assert(isnumeric(year),'Input ''year'' has to be numeric.')
      fix_idx=year>=time.millennium;
      if any(fix_idx)
        year(fix_idx)=year(fix_idx)-time.millennium;
      end
    end
    %% operations
    function out=union(time_domains,time_step)
      switch class(time_domains)
        case 'cell'
          mint=cellfun(@(i) i(1),  time_domains,'UniformOutput',false);mint=min([mint{:}]);
          maxt=cellfun(@(i) i(end),time_domains,'UniformOutput',false);maxt=max([maxt{:}]);
          out=mint:time_step:maxt;
        case 'datetime'
          out=time_domains(1):time_step:time_domains(end);
        otherwise
          error(['Cannot handle input ''time_domains'' of class ''',class(time_domains),'''.'])
      end
    end
  end
end