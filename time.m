classdef time
  %static
  properties(Constant,GetAccess=private)
    formt_list={'S.FFF','SS.FFF','HH:MM:SS','HH:MM','dd/mm HH:MM','dd/mm/yy','mm/yyyy' ,'yyyy'}
    units_list={'ms'   ,'s'     ,'min'     ,'hrs'  ,'days'       ,'mon'     ,'yrs'     ,'cent'};
    magnt_list=[0.01   ,60      ,3600      ,86400  ,2678400      ,31536000  ,3153600000,inf   ];
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
    function s=progress(s,i)
      % %Example: 
      % s.msg='something'; s.n=n;
      % for i=1:n
      %     <do something>
      %     s=simpledata.progress(s,i);
      % end
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
        fprintf('%s\n',sprintf('%s',[repmat(8,1,s.l_max),s.msg,': done in ',time.str(s.e(end)),'                            ']))
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
  end
end