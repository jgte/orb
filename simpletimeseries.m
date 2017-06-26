classdef simpletimeseries < simpledata
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list=struct(...
      'format',    struct('default','modifiedjuliandate','validation',@(i) ischar(i)),...
      'units',     struct('default',{{''}},              'validation',@(i) iscellstr(i)),...
      't_tol',     struct('default',1e-6,                'validation',@(i) isnumeric(i) && isscalar(i)),...
      'timesystem',struct('default','utc',               'validation',@(i) ischar(i)),...
      'debug',     struct('default',false,               'validation',@(i) islogical(i) && isscalar(i))...
    );
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'timesystem'};
  end
  properties(Constant)
    valid_timesystems={'utc','gps'};
  end
  %NOTE: edit this if you add a new parameter (if read only)
  properties(SetAccess=private)
    step
    debug
    timesystem
  end
 %These parameters should not modify the data in any way; they should
  %only describe the data or the input/output format of it.
  %NOTE: edit this if you add a new parameter (if read/write)
  properties(GetAccess=public,SetAccess=public)
    format
    t_tol
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
    function out=valid_timesystem(in)
      switch lower(in)
      case simpletimeseries.valid_timesystems
        out=true;
      otherwise
        out=false;
      end
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
    function out=ist(mode,t1,t2,tol)
      switch mode
      case {'=','==','equal'}
        if numel(t1)==numel(t2) 
          out=seconds(t1(:)-t2(:)).^2<tol.^2;
        elseif isscalar(t1)
          out=seconds(t1-t2(:)).^2<tol.^2;
        elseif isscalar(t2)
          out=seconds(t1(:)-t2).^2<tol.^2;
        else
          out=false;
        end
        return
      case {'<','less','smaller'}
        out=t1<t2;
        out(simpletimeseries.ist('==',t1,t2,tol))=false;
      case {'<=','lessorequal'}
        out=t1<t2;
        out(simpletimeseries.ist('==',t1,t2,tol))=true;
      case {'>','more','larger'}
        out=t1>t2;
        out(simpletimeseries.ist('==',t1,t2,tol))=false;
      case {'>=','moreorequal','largerorequal'}
        out=t1>t2;
        out(simpletimeseries.ist('==',t1,t2,tol))=true;
      otherwise
        error([mfilename,': unknown mode ''',mode,'''.'])
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
        error([mfilename,': input argument ''fields'' must be a cell array.'])
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
          error([mfilename,': cannot handle both inputs ''x'' and ''t''.'])
        end
      end
    end
    function out=parameters
      out=fieldnames(simpletimeseries.parameter_list);
    end
    function out=transmute(in)
      if isa(in,'simpletimeseries')
        %trivial call
        out=in;
      else
        %transmute into this object
        if isprop(in,'t')
          out=simpletimeseries(in.t,in.y,in.metadata{:});
        elseif isprop(in,'x')
          out=simpletimeseries(in.x,in.y,in.metadata{:});
        else
          error('Cannot find ''t'' or ''x''. Cannot continue.')
        end
      end
    end
    function out=timestep(in,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'in',                @(i) isdatetime(i));
      p.addParameter('nsigma',    4,      @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('max_iter',  10,     @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('sigma_iter',2,      @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('sigma_crit',1e-9,   @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('max_mean_ratio',1e3,@(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('curr_iter', 0,      @(i) isnumeric(i)  &&  isscalar(i));
      p.addParameter('disp_flag', false,  @(i) islogical(i));
      % parse it
      p.parse(in,varargin{:});
      %handle singularities
      switch numel(in)
        case 0
          error([mfilename,': cannot handle empty time stamps'])
        case 1
          out=0;
          return
      end
      %get numeric diff of time
      tdiff=simpletimeseries.timescale(diff(in));
      %large jumps produce erroneous results, so get rid of those first
      while std(tdiff)~=0 && max(tdiff)/mean(tdiff)>p.Results.max_mean_ratio
        %send feedback
        if p.Results.disp_flag
          disp([mfilename,': removing large gaps, since max(delta t) is ',num2str(max(diff(tdiff))),...
            ' and is ',num2str(max(diff(tdiff))/mean(diff(tdiff))),' times larger than mean(delta).'])
        end
        tdiff=simpledata.rm_outliers(tdiff,varargin{:});
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
    %general test for the current object
    function test(l,w)
      if ~exist('l','var') || isempty(l)
        l=1000;
      end
      if ~exist('w','var') || isempty(w)
        w=3;
      end
      %test current object
      args=simpledata.test_parameters('args',l,w);
      now=juliandate(datetime('now'),'modifiedjuliandate');
                
      i=0;      
            
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
      
      a=simpledata.test_parameters('all_T',l,w);
      i=i+1;h{i}=figure('visible','off'); a.plot('title', 'parametric decomposition','columns',1)
      b=a.parametric_decomposition(...
        'polynomial',ones(size(simpledata.test_parameters('y_poly_scale'))),...
        'sinusoidal',simpledata.test_parameters('T',l)...
      );
      fn=fields(b);tot=[];legend_str={};
      for i=1:numel(fn);
        if ~isempty(strfind(fn{i},'ts_'))
          legend_str{end+1}=fn{i}; %#ok<*AGROW>
          b.(fn{i}).plot('columns',1)
          if isempty(tot)
            tot=b.(fn{i});
          else
            tot=tot+b.(fn{i});
          end
        end
      end
      tot.plot('columns',1)
      legend_str=strrep(legend_str,'_','\_');
      legend('original',legend_str{:},'sum')
      
      for i=numel(h):-1:1
        set(h{i},'visible','on')
      end

    end
    %% import methods
    function filenames=unwrap_datafiles(in,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired(  'in',                                        @(i) ischar(i) || iscellstr(i));
      p.addParameter( 'start',       time.ToDateTime(0,'datenum'), @(i) isscalar(i) && isdatetime(i));
      p.addParameter( 'stop',        time.ToDateTime(0,'datenum'), @(i) isscalar(i) && isdatetime(i));
      p.addParameter( 'period',      days(1),                      @(i) isscalar(i) && isduration(i));
      p.addParameter( 'date_fmt',    'yyyy-mm-dd',                 @(i) ischar(i));
      p.addParameter( 'only_existing',true,                        @(i) islogical(i));
      p.parse(in,varargin{:})
      %loop over all inputs
      if iscellstr(in)
        filenames=cell(size(in));
        for i=1:numel(in);
          filenames{i}=simpletimeseries.unwrap_datafiles(in{i},varargin{:});
        end
        filenames=cells.flatten(filenames);
        %remove empty entries
        filenames=cells.rm_empty(filenames);
        return
      end
      %need fileparts
      [d,f,e]=fileparts(in);
      
      %if argument filename has a wild card, build a cell string with those names
      if ~isempty(strfind(in,'*'))
        file_list=dir(in);
        filenames=cell(size(file_list));
        for i=1:numel(file_list)
          filenames{i}=fullfile(d,file_list(i).name);
        end
        %unwrap these files again
        filenames=simpletimeseries.unwrap_datafiles(filenames,varargin{:});
        %done
        return
      end
      
      %if argument filename has a date place holder, build a cell string with those dates
      if ~isempty(strfind(in,'DATE_PLACE_HOLDER'))
        %sanity
        if p.Results.start>=p.Results.stop
          error([mfilename,': input ''start'' (',datestr(p.Results.start),...
            ') is not after input ''stop'' (',datestr(p.Results.stop),').'])
        end
        date_list=time.list(p.Results.start,p.Results.stop,p.Results.period);
        filenames=cell(size(date_list));
        for i=1:numel(date_list)
          filenames{i}=strrep(in,'DATE_PLACE_HOLDER',datestr(date_list(i),p.Results.date_fmt));
        end
        %unwrap these files again
        filenames=simpletimeseries.unwrap_datafiles(filenames,varargin{:});
        %done
        return
      end
      
      %if there's a .mat file along the argument filename, pass that on
      matfile=fullfile(d,[f,'.mat']);
      if ~isempty(dir(matfile))
        filenames=matfile;
        return
      end
      
      %now handling compressed files: prepend tar to extension if there
      if strcmp(f(end-3:end),'.tar')
        e=['.tar',e];
      end
      %try to uncompress archives
      try
        switch lower(e)
        case {'.z','.zip'}
          arch=true;
          filenames=unzip(in,d);
        case {'.tgz','.tar.gz','.tar'}
          arch=true;
          filenames=untar(in,d);
          %get rid of PaxHeaders
          filenames(cells.iscellstrfind(filenames,'PaxHeaders'))=[];
        case {'.gz','.gzip'}
          arch=true;
          filenames=gunzip(in,d);
        otherwise
          arch=false;  
        end
        if arch
          disp(['Extracted archive ''',in,'''.'])
        end
      catch
        %if the zip file is corrupted, assume data file is missing
        disp(['WARNING: error extracting archive ''',in,'''.'])
        return
      end
      %handle zipped files
      if arch
        %some sanity
        if ~iscell(filenames)
          error([mfilename,': expecting variable ''unzipped_filename'' to be a cellstr, not a ''',...
            class(filenames),'''.'])
        end
        if numel(filenames)~=1
          error([mfilename,': expecting zip archive ''',filenames,''' to contain one file only, not ',...
            num2str(numel(filenames)),':',10,strjoin(filenames,'\n')])
        end
        %and we're done (no recursive unwrapping!)
        return
      end
      
      %if none of the conditions above were met, this is the name of a single file (return char!)
      if isempty(dir(in)) && p.Results.only_existing
        filenames='';
      else
        filenames=in;
      end
    end
    function obj=import(filename,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter( 'save_mat', true,  @(i) isscalar(i) && islogical(i))
      p.addParameter( 'cut24hrs', true,  @(i) isscalar(i) && islogical(i))
      p.addParameter( 'del_arch', true,  @(i) isscalar(i) && islogical(i))
      p.parse(varargin{:})
      %use this flag to skip saving mat data (e.g. if input file is empty)
      skip_save_mat=false;
      %unwrap wildcards and place holders (output is always a cellstr)
      filename=simpletimeseries.unwrap_datafiles(filename,varargin{:});
      %if argument is a cell string, then load all those files
      if iscellstr(filename)
        for i=1:numel(filename)
          disp([mfilename,': reading data from file ',filename{i}])
          %read the data from a single file
          obj_now=simpletimeseries.import(filename{i});
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
        return
      end
      %split into parts
      [d,f,e]=fileparts(filename);
      %check if mat file is available
      datafile=fullfile(d,[f,'.mat']);
      if ~isempty(dir(datafile))
        load(datafile)
        %sanity on the loaded data
        if ~exist('obj','var')
          error([mfilename,': expecting to load variables ''obj'' from file ',datafile,'.'])
        end
        %we're done
        return
      end
      %some files have the format ID in front
      for i={'ACC1B','SCA1B','KBR1B','GNV1B','grc[AB]_gps_orb_.*\.acc'}
        if ~isempty(regexp(filename,i{1},'once'))
          e=i{1};
          break
        end
      end
      %branch on extension/format ID
      switch e
      case '.resid'
        fid=fopen(filename);
        raw = textscan(fid,'%f %f %f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne',1,'Headerlines',1);
        fclose(fid);
        %building time domain
        t=time.utc2gps(datetime(raw{1},...
          'convertfrom','epochtime',...
          'epoch','2000-01-01'...
        ));
        %building data domain
        y=[raw{5:7}];
        %determine coordinate
        coords={'AC0X','AC0Y','AC0Z'};
        idx=cells.iscellstrfind(filename,coords);
        %sanity
        assert(sum(idx)==1,[mfilename,': the name for .resid files must include (one of) AC0X, AC0Y or AC0Z, not ''',...
          filename,'''.'])
        %building object
        obj=simpletimeseries(t,y,...
          'format','datetime',...
          'y_units',{'m^2','m^2','m^2'},...
          'labels', {coords{idx},[coords{idx},'D'],[coords{idx},'Q']},...
          'timesystem','gps',...
          'descriptor',['model response from file ',filename]...
         );
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
          'timesystem','utc',...
          'descriptor',['kinematic orbit from file ',filename]...
         );
      case '.GraceAccCal'
        fmt='';
        if ~isempty(regexp(filename,'AC0[XYZ]\d?\.aak','once')) || ~isempty(regexp(filename,'AC0[XYZ]\d?\.accatt','once'))
          % 2002 04 05 2002.4.4. 23.59.47.00000000 1498260002 0.2784215319157E-07
          fmt='%d %d %d %s %s %d %f';
          units={'m/s^2',''};
          labels={str.clean(filename,{'file','grace','.'}),'Job ID','arc start'};
          time_fh=@(raw) time.utc2gps(...
            datetime(...
              strcat(...
                strrep(cellfun(@(x) [x(1:end-1),' '],raw{4},'UniformOutput',false),'.','/'),...
                strrep(strrep(raw{5},'.00000000',''),'.',':')...
              ),'InputFormat','yyyy/MM/dd HH:mm:ss'...
            )...
          );
          data_fh=@(raw) [raw{7},double(raw{6})];
          timesystem='gps';
          sanity_check=@(raw) true;
        end
        if ~isempty(regexp(filename,'AC0[XYZ][QD]\d?\.aak','once')) || ~isempty(regexp(filename,'AC0[XYZ][QD]\d?\.accatt','once'))
          % 2002 04 05 2002.4.4. 23.59.47.00000000 1498260002  0.1389481692269E-07 52368.99985
          fmt='%d %d %d %s %s %d %f %f';
          units={'m/s^2','','MJD days'};
          labels={str.clean(filename,{'file','grace','.'}),'Job ID','t_0','arc start'};
          time_fh=@(raw) time.utc2gps(...
            datetime(...
              strcat(...
                strrep(cellfun(@(x) [x(1:end-1),' '],raw{4},'UniformOutput',false),'.','/'),...
                strrep(strrep(raw{5},'.00000000',''),'.',':')...
              ),'InputFormat','yyyy/MM/dd HH:mm:ss'...
            )...
          );
          data_fh=@(raw) [raw{7},double(raw{6}),raw{8}];
          timesystem='gps';
          sanity_check=@(raw) true;
        end
        if ~isempty(regexp(filename,'AC0[XYZ]\d?\.estim','once')) || ~isempty(regexp(filename,'AC0[XYZ][DQ]\d?\.estim','once'))
          % 1    2  3  4  5  6  7     8 9   10      11       12                    13                     14                    15           16
          % 2002 04 05 04/05/02 52369 1 0.0 26400.0 1593715  3.774424464092000e-08 -3.585594302740665e-09 3.415865033817934e-08 2.82822E-09  71279987.
          fmt='%d %d %d %d/%d/%d %f %d %f %f %d %f %f %f %f %f';
          units={'m/s^2','','sec','sec','mjd','','m/s^2'};
          labels={str.clean(filename,{'file','grace','.'}),'Job ID','arc duration','arc start','arc t0','arc nr','TBD'};
          time_fh=@(raw) datetime(raw{7}+raw{9}/seconds(days(1)),...
            'ConvertFrom','modifiedjuliandate'...
          );
          data_fh=@(raw) [raw{14},double(raw{11}),raw{10},raw{9},time.mjd(time.utc2gps(time.ToDateTime(double(raw{16}),'J2000sec'))),double(raw{8}),double(raw{15})];
          timesystem='gps';
          sanity_check=@(raw) all(all([ raw{1}-2000==raw{6},raw{2}==raw{4},raw{3}==raw{5}]));
        end

        if isempty(fmt)
          error([mfilename,': cannot handle the GraceAccCal file ''',filename,'''.'])
        end
        %reading data
        fid = fopen(filename);
        raw = textscan(fid,fmt,'delimiter',' ','MultipleDelimsAsOne',1);
        fclose(fid);
        %keep some sanity
        assert(sanity_check(raw),['Failed sanity check on ',filename]);
        %building time domain
        t=time_fh(raw);
        %building data domain
        y=data_fh(raw);
        %sanity
        if isempty(t) || isempty(y)
          disp([mfilename,': this file has no data  ',filename])
          skip_save_mat=true;
          obj=[];
        else
          iter=0;
          while any(diff(t)==0)
            %loop inits
            n0=numel(t);
            iter=iter+1;
            %need to remove duplicate entries with different job IDs
            mask=true(size(t));
            for i=2:numel(t)
              %get rid of those entries with zero or negative time stamp delta and lower ID
              if t(i)<=t(i-1) && mask(i)
                if y(i,2) > y(i-1,2)
                  mask(i-1)=false;
                else
                  mask(i)=false;
                end
              end
            end
            t=t(mask);
            y=y(mask,:);
            disp(['At iter ',num2str(iter),', removed ',num2str(n0-numel(t),'%04d'),' duplicate time entries (',filename,').'])
          end
          %need to monotonize the data (sometimes the entries are ordered according to arc number and not chronologically)
          if any(diff(t)<0)
            [t,i]=sort(t);
            y=y(i,:);
            disp(['Sorted ',num2str(sum(i~=transpose(1:numel(i))),'%04d'),' time entries (',filename,').'])
          end
          %building object
          obj=simpletimeseries(t,y,...
            'format','datetime',...
            'labels',labels,...
            'units',units,...
            'timesystem',timesystem,...
            'descriptor',filename...
          );
        end
      case 'ACC1B'
        %load data
        [raw,header]=file.textscan(filename,'%f %s %f %f %f %f %f %f %f %f %f %f');
        %retrieve GPS time epoch
        header_line='TIME EPOCH (GPS TIME)         : ';
        header=strsplit(header,'\n');
        for i=1:numel(header)
          if strfind(header{i},header_line)
            gps_time_epoch=strrep(header{i},header_line,'');
            break
          end
        end
        %building time domain
        t=time.gpssec2datetime(raw(:,1),gps_time_epoch);
        %gather data domain
        y=raw(:,2:4);
        %skip empty data files
        if isempty(t) || isempty(y)
          disp([mfilename,': this file has no data  ',filename])
          skip_save_mat=true;
          obj=[];
        else
          %building object
          obj=simpletimeseries(t,y,...
            'format','datetime',...
            'y_units',{'m/s^2','m/s^2','m/s^2'},...
            'labels', {'x','y','z'},...
            'timesystem','gps',...
            'descriptor',strjoin(header,'\n')...
           ).fill;
        end
      case 'grc[AB]_gps_orb_.*\.acc'
        %load data
        [raw,header]=file.textscan(filename,'%f %f %f %f %f %f %f %f',[],'%');
        %retrieve GPS time epoch
        header_line='+unitfacor ';
        header=strsplit(header,'\n');
        for i=1:numel(header)
          if strfind(header{i},header_line)
            unitfactor=str2double(strrep(header{i},header_line,''));
            break
          end
        end
        %building time domain
        t=time.ToDateTime(raw(:,1:3),'yeardoysec');
        %gather data domain
        y=raw(:,5:7)./unitfactor;
        %skip empty data files
        if isempty(t) || isempty(y)
          disp([mfilename,': this file has no data  ',filename{i}])
          skip_save_mat=true;
          obj=[];
        else
          %building object
          obj=simpletimeseries(t,y,...
            'format','datetime',...
            'y_units',{'m/s^2','m/s^2','m/s^2'},...
            'labels', {'x','y','z'},...
            'timesystem','gps',...
            'descriptor',strjoin(header,'\n')...
           ).fill;
        end
      case 'SCA1B'
        error([mfilename,': implementation needed'])
      case 'KBR1B'
        error([mfilename,': implementation needed'])
      case 'GNV1B'
        error([mfilename,': implementation needed'])
      otherwise
        error([mfilename,': cannot handle files of type ''',e,'''.'])
      end
      %save mat file if requested
      if p.Results.save_mat && ~skip_save_mat
        save(datafile,'obj')
      end
      %delete uncompressed file if compressed file is there
      if p.Results.del_arch
        for i={'.z','.zip','.tgz','.gz','.tar','.gzip'}
          if ~isempty(dir(fullfile(d,[f,i{1}])))
            delete(filename)
            disp(['Deleted uncompressed file ''',in,'''.'])
          end
        end
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
    function [t,y]=monotonize(t,y)
      %get delta time
      dt=diff(t(:));
      %trivial call
      if all(dt>0)
        return
      end
      %sorting needed?
      if any(dt<0)
        error('implementation needed')
      end
      %get index of duplicates (last one is not duplicate by definition)
      idx0=[(dt==0); false];
      %get index+1 of duplicates
      idx1=circshift(idx0,1,1);
      %make sure the values are the same
      for i=1:size(y,2)
        if any(y(idx0,i)~=y(idx1,i))
          error([mfilename,': cannot monotonize the data because column ',num2str(i),' has different data in common epochs.'])
        end
      end
      %throw away duplicate epochs
      t=t(~idx0);
      y=y(~idx0,:);
    end
    function out=stats2(obj1,obj2,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('period', 30*max([obj1.step,obj2.step]), @(i) isduration(i));
      p.addParameter('overlap',seconds(0),                    @(i) isduration(i));
      % parse it
      p.parse(varargin{:});
      % call upstream method if period is infinite
      if ~isfinite(p.Results.period)
        out=simpledata.stats2(obj1,obj2,varargin{:});
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
        dat(i)=simpledata.stats2(...
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
  end
  methods
    %% constructor
    function obj=simpletimeseries(t,y,varargin)
      % parameter names
      pn=simpletimeseries.parameters;
      % input parsing
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y', @(i) simpledata.valid_y(i));
      %declare parameters
      for i=1:numel(pn)
        %declare parameters
        p.addParameter(pn{i},simpletimeseries.parameter_list.(pn{i}).default,simpletimeseries.parameter_list.(pn{i}).validation)
      end
      % parse it
      p.parse(t,y,varargin{:});
      % get datetime 
      [t,f]=time.ToDateTime(t,p.Results.format);
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
      for i=1:numel(pn)
        %parameter 'units' has already been handled when calling simpledata constructor, so skip it
        if strcmp(pn{i},'units')
          continue
        end
        if ~isscalar(p.Results.(pn{i}))
          %vectors are always lines (easier to handle strings)
          obj.(pn{i})=transpose(p.Results.(pn{i})(:));
        else
          obj.(pn{i})=p.Results.(pn{i});
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
    function out=metadata(obj)
      %call superclass
      out=metadata@simpledata(obj);
      %build cell array with parameters of this object
      parameters=simpletimeseries.parameters;
      for i=1:numel(parameters)
        if isprop(obj,parameters{i})
          out{end+1}=parameters{i};
          out{end+1}=obj.(parameters{i});
        end
      end
    end
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
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('period', seconds(inf), @(i) isduration(i));
      p.addParameter('overlap',seconds(0),   @(i) isduration(i));
      p.addParameter('mode',  'struct',      @(i) ischar(i));
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
        dat(i)=stats@simpledata(obj.trim(ts{i}(1),ts{i}(end)),varargin{:},'mode','struct'); %#ok<AGROW>
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
          case {'min','max','mean','std','rms','meanabs','stdabs','rmsabs'}; units=obj.y_units;
          case {'length','gaps'};                                            units=repmat({' '},1,obj.width);
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
    function out=span(obj)
      out=obj.stop-obj.start;
    end
    function out=t_domain(obj,step_now)
      if ~exist('step_now','var') || isempty(step_now)
        step_now=obj.step;
      end
      out=transpose(obj.start:step_now:obj.stop);
    end
    function out=ishomogeneous(obj)
      htd=obj.t_domain;
      out=(numel(htd)==numel(obj.t)) && all(obj.t(:)==htd(:));
    end
    function out=istavail(obj,t)
      if isscalar(t)
        out=any(simpletimeseries.ist('==',t,obj.t,obj.t_tol));
      else
        for i=1:numel(t)
          %scalar call
          out=obj.istavail(t(i));
          %no need to continue looping if found something
          if out;break;end
        end
      end
    end
    function out=isxavail(obj,t)
      %need to overload isxavail so that calls from simpledata come through
      %here and not through the function defined there
      out=obj.istavail(t);
    end
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
    function out=t_masked(obj,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      end
      out=obj.t(mask);
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
      if any(~simpletimeseries.ist('==',t_old,obj.t,obj.t_tol))
        error([mfilename,': changing epoch cause the time domain to also change.'])
      end
    end
    function out=get.epoch(obj)
      out=obj.epochi;
    end
    function obj=epoch_update(obj)
      obj.epoch=obj.t(1);
    end
    %% start/stop methods
    function out=get.start(obj)
      out=obj.t(1);
    end
    function out=get.stop(obj)
      out=obj.t(obj.length);
    end
    function obj=set.start(obj,start)
      if isempty(start) || simpletimeseries.ist('==',start,obj.start,obj.t_tol)
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
      if isempty(stop) || simpletimeseries.ist('==',stop,obj.stop,obj.t_tol)
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
        error([mfilename,': need a valid time system, i.e. one of ',strjoin(simpletimeseries.valid_timesystems,', '),'.'])
      end
      obj.t=simpletimeseries.([obj.timesystem,'2',in])(obj.t);
      obj.timesystem=in;
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
        if any(~simpletimeseries.ist('==',obj.t,t_now,obj.t_tol))
          error([mfilename,': the time domain is not consistent with input ''t_now''.'])
        end
      end
    end
    %% edit methods (overloaded with simpledata)
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
      %save current time domain and step
      t_now=obj.t;
%       step_now=obj.step;
      %call superclass
      obj=median@simpledata(obj,n);
      %resample (if needed, which is checked inside resample)
      obj=obj.interp(t_now,...
        'interp_over_gaps_narrower_than',0,...
        'interp1_args',{'linear'}...
      );
    end
    %% edit methods (specific to this class)
    function obj=extend(obj,nr_epochs)
      %sanity
      if ~obj.ishomogeneous
        error([mfilename,': cannot handle non-homogeneous time domains.'])
      end
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
        obj=extend(obj,floor((t_now-t_ref)/obj.step));
      otherwise
        error([mfilename,': cannot handle input ''nr_epochs'' of class ',class(nr_epochs),'.'])
      end
    end
    function obj=append_epochs(obj,t_now,y_now)
      %shortcuts
      if isscalar(y_now)
        y_now=y_now*ones(numel(t_now),obj.width);
      end
      %sanity
      assert(numel(t_now)==size(y_now,1) || size(y_now,2)~=obj.width,[mfilename,...
        'inputs ''t_now'' and/or ''y_now'' have sizes inconsistent with this obj.'])
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
      %        but it differs since the added time entries are set as explicit gaps.
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
        error([mfilename,': complete time domain has less entries than current time domain. Debug needed!'])
      end
      %find out where there are gaps larger than the step size
      gap_idx=find(diff(obj.t)>obj.step);
      %if there are no gaps and the time series is not homogeneous, we have a problem that needs fixing
      if isempty(gap_idx)
        error([mfilename,': implementation needed!'])
      end
      disp(['Need to fill in missing epochs: ',num2str(numel(t_new)-obj.length),' ('...
        num2str((numel(t_new)-obj.length)/numel(t_new)*1e2),'%).'])
      %loop over all implicit gaps (i.e. missing epochs)
      s.msg=[mfilename,': populating missing epochs (',datestr(obj.start),' to ',datestr(obj.stop),')',...
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
    function [obj_clean,obj_outlier]=despike(obj,n,nSigma)
      %get medianed timeseries
      obj_median=obj.median(n);
      %compute residual to median
      obj_res=obj-obj_median;
      %remove outliers from residual
      [obj_res_clean,obj_res_outlier]=obj_res.outlier(nSigma);
      %restore median
      obj_clean=obj_median+obj_res_clean;
      obj_outlier=obj_median+obj_res_outlier;
    end
    %% multiple object manipulation
    function out=isteq(obj1,obj2)
      out=~any(~simpletimeseries.ist('==',obj1.t,obj2.t,min([obj1.t_tol,obj2.t_tol])));
    end
    function compatible(obj1,obj2,varargin)
      %call mother routine
      compatible@simpledata(obj1,obj2,varargin{:});
      %shorter names
      par=simpletimeseries.compatible_parameter_list;
      for i=1:numel(par)
        % if a parameter is empty, no need to check it
        if ( iscell(obj1.(par{i})) && isempty([obj1.(par{i}){:}]) ) || ...
           ( ischar(obj1.(par{i})) && isempty( obj1.(par{i})    ) ) || ...
           ( iscell(obj2.(par{i})) && isempty([obj2.(par{i}){:}]) ) || ...
           ( ischar(obj2.(par{i})) && isempty( obj2.(par{i})    ) )
          continue
        end
        if ~isequal(obj1.(par{i}),obj2.(par{i}))
          error([mfilename,': discrepancy in parameter ',par{i},'.'])
        end 
      end
    end
    function [obj1,obj2,idx1,idx2]=merge(obj1,obj2)
      %add as gaps the t in obj1 that are in obj2 but not in obj1 (and vice-versa)
      %NOTICE:
      % - idx1 contains the index of the x in obj1 that were added from obj2
      % - idx2 contains the index of the x in obj2 that were added from obj1
      % - no data is propagated between objects, only the time domain is changed!
      if isprop(obj1,'epoch') && isprop(obj2,'epoch')
        [obj1,obj2]=matchepoch(obj1,obj2);
      end
      %call upstream method
      [obj1,obj2,idx1,idx2]=merge@simpledata(obj1,obj2);
      %sanity
      if ~isteq(obj1,obj2)
        error([mfilename,':BUG TRAP: failed to merge time domains. Debug needed!'])
      end
    end
    function [obj1,obj2]=interp2(obj1,obj2,varargin)
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
      %need to match the epoch
      if isa(obj1,'simpletimeseries') && isa(obj2,'simpletimeseries')
        [obj1,obj2]=matchepoch(obj1,obj2);
      end
      %call upstream method
      [obj1,obj2]=interp2@simpledata(obj1,obj2,varargin{:});
      %sanity
      if ~isteq(obj1,obj2)
        error([mfilename,':BUG TRAP: failed to merge time domains. Debug needed!'])
      end
    end
    function [obj,idx1,idx2]=append(obj1,obj2)
      if isa(obj1,'simpletimeseries') && isa(obj2,'simpletimeseries')
        [obj1,obj2]=matchepoch(obj1,obj2);
      end
      %call upstream method
      [obj,idx1,idx2]=append@simpledata(obj1,obj2);
    end
    function obj1_out=augment(obj1,obj2,varargin)
      if isa(obj1,'simpletimeseries') && isa(obj2,'simpletimeseries')
        [obj1,obj2]=matchepoch(obj1,obj2);
      end
      %call upstream method
      obj1_out=augment@simpledata(obj1,obj2,varargin{:});
    end
    %NOTICE: this function used to be called consolidate
    function [obj1,obj2]=interp2_lcm(obj1,obj2)
      %extends the time domain of both objects to be in agreement
      %with the each other
      compatible(obj1,obj2)
      %trivial call
      if isteq(obj1,obj2)
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
    function [obj1,obj2]=matchstep(obj1,obj2)
      %sanity
      if ~obj1.ishomogeneous || ~obj2.ishomogeneous
        error([mfilename,': can only handle homogeneous time domains.'])
      end
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
    function [obj1,obj2]=matchepoch(obj1,obj2)
      %trivial call
      if obj1.epoch==obj2.epoch
        return
      end
      %match epochs
      obj2.epoch=obj1.epoch;
    end
    function [obj1,obj2]=matchtime(obj1,obj2)
      %match step and epoch (checks for trivial call done inside)
      [obj1,obj2]=matchstep(obj1,obj2);
      [obj1,obj2]=matchepoch(obj1,obj2);
    end
    function obj1=glue(obj1,obj2)
      %objects need to have the same epoch
      assert(obj1.epoch==obj2.epoch,...
        'Input objects do not share the same time domain.')
      %call mother routine
      obj1=glue@simpledata(obj1,obj2);
    end
    %% plot methots
    function out=plot(obj,varargin)
      %call superclass
      out=plot@simpledata(obj,varargin{:});
%       %using internal Matlab representation for dates
%       lines_now=get(gca,'children');
%       for i=1:numel(lines_now)
%         for i=1:numel(out.handle)
%           if out.handle{i}==lines_now(i)
%             lines_now(i).XData=datenum(obj.x2t(lines_now(i).XData));
%           end
%         end
%       end
%       if strcmpi(get(gca,'XTick'),'auto')
%         set(gca,'XTick',datenum(obj.t));
%         datetick('x',time.format(seconds(obj.span)))
%       end
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
'# Column 3:    Time system (',obj.timesystem,')',10,...
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
        mjd=obj.mjd(obj.mask);
        %build format string
        fmt='%s UTC %14.8f';
        for i=1:numel(p.Results.columns)
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