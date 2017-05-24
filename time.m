classdef time
  %static
  properties(Constant,GetAccess=private)
    formt_list={'S.FFF','SS.FFF','HH:MM:SS','HH:MM','dd/mm HH:MM','dd/mm/yy','mm/yyyy' ,'yyyy'}
    units_list={'ms'   ,'s'     ,'min'     ,'hrs'  ,'days'       ,'mon'     ,'yrs'     ,'cent'};
    magnt_list=[0.01   ,60      ,3600      ,86400  ,2678400      ,31536000  ,3153600000,inf   ];
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
  end    
  properties(Constant)
    zero_date=datetime(0,  'ConvertFrom','datenum');
     inf_date=datetime(Inf,'ConvertFrom','datenum');
  end
  methods(Static)
    function test
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
    function out=units(in)
      out=time.units_list{time.scale(in)};
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
        disp([mfilename,':WARNING: unknown time units ''',time_units,'''. Using seconds.'])
        out = 's';
      end
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
    function out = str(seconds)
      if ~isscalar(seconds)
        error([mfilename,': input must be scalar'])
      end
      if isnan(seconds)
        out='NaN';
      elseif ~isfinite(seconds)
        out='Inf';
      else
        steps = [1,60,60*60,60*60*24,60*60*24*7,60*60*24*31,60*60*24*365];
        units = ['s','m','h','d','w','m','y'];
        out=[];
        seconds = abs(seconds);
        for i=numel(steps):-1:1
          if (seconds > steps(i))
            if i~=1
              out = [out,num2str(floor(seconds/steps(i))),units(i),',']; %#ok<AGROW>
            else
              out = [out,num2str(seconds,3),units(i)]; %#ok<AGROW>
            end
            seconds = rem(seconds,steps(i));
          elseif (i==1)
            out = [out,num2str(seconds,3),units(i)]; %#ok<AGROW>
          end
        end
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
    end
    function s=progress(s,i)
      % %Example: 
      % s.msg='something'; s.n=n;
      % for i=1:n
      %     <do something>
      %     s=time.progress(s,i);
      % end
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
        error([mfilename,': need structure in argument ''s'' with field ''n'' (the loop total nr of iterations).']);
      end
      %defaults
      if ~isfield(s,'step')
        s.step=0.01;
      end
      if ~isfield(s,'msg')
        s.msg='';
      end
      %initialize
      if ~isfield(s,'l')
        s.l=0;
        s.l_max=0;
        s.step_c=0;
        tic
        s.t=zeros(1,ceil(1/s.step));
        s.i=zeros(1,ceil(1/s.step));
        s.c=0;
      end
      %finalize
      if i>=s.n
        if s.c>0
          fprintf('%s\n',sprintf('%s',[...
            repmat(8,1,s.l_max),s.msg,': ',...
            num2str(i,'%03d'),' iter done in ',time.str(s.e(end)),...
            '                            '...
          ]))
        end
        return
      end
      %give feedback during iterations
      if (i > s.step_c)
        dt=toc;
        s.c=s.c+1;
        s.i(s.c)=i;
        s.t(s.c)=dt;
        s.e(s.c)=sum(s.t(1:s.c));
        out=sprintf('%s',[repmat(8,1,s.l),s.msg,': ',num2str(round(i/s.n*100)),'% (',...
          'finished in ',time.str( (1/s.step-s.c)*mean(s.t(1:s.c)) ),...
          ', elapsed ',time.str( s.e(s.c) ),...
          ')               ']);
        fprintf('%s',out)
        s.l=numel(out);      %nr of chars to be deleted
        s.l_max=max([s.l_max,s.l]);
        s.step_c=i+s.step*s.n;   %next user feedback at this iter
        tic
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
    function date=doy2datetime(year,doy)
      date=datetime(datenum([year 0 0 0 0 0])+doy,'ConvertFrom','datenum');
    end
    function out=rand(n)
      out=datetime(sort(datenum(round([...
        rand(n,1),...
        rand(n,1)*12,...
        rand(n,1)*31,...
        rand(n,3)*60 ...
      ]))),'ConvertFrom','datenum')+years(2000);
    end
    function [startlist,stoplist]=day_list(start,stop)
      p=inputParser;
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
      p=inputParser;
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
      p=inputParser;
      p.addRequired( 'start',   @(i) isscalar(i) && isdatetime(i));
      p.addRequired( 'stop',    @(i) isscalar(i) && isdatetime(i));
      p.parse(start,stop)
      %https://www.mathworks.com/matlabcentral/answers/73749-how-do-i-create-a-vector-with-the-first-day-of-each-month#answer_83695
      dvec      = datevec(datenum(start):datenum(stop));
      duniq     = unique(dvec(:, 1), 'rows');
      startlist = datetime(datenum(duniq(:,1), 1, 1),'ConvertFrom','datenum');
      if nargout>1
        stoplist = [startlist(2:end);dateshift(startlist(end),'end','year')+days(1)];
      end
    end
    function out=isvalid(in)
      if isempty(in)
        out=false;
      else
        switch class(in)
        case 'cell'
          out=cellfun(@time.isvalid,in);
        case 'datetime'
          out= in~=time.zero_date;
        otherwise
          error([mfilename,': unsupported class ',class(in),'.'])
        end
      end
    end
  end
end