classdef orbit
  %static
  properties(Constant,GetAccess=private)
    %list of data fields
    %NOTICE: needs updated when adding a new data type
    data_type_list=struct(...
      'pos',struct(...
        'label','position',...
        'size',3,...
        'xyz',struct('units',{{'m',  'm',  'm'}},'names',{{'x',  'y',  'z'}}),...
        'sph',struct('units',{{'deg','deg','m'}},'names',{{'lon','lat','rad'}})),...
      'vel',struct(...
        'label','velocity',...
        'size',3,...
        'xyz',struct('units',{{'m/s',  'm/s',  'm/s'}},'names',{{'x',   'y',   'z'}}),...
        'sph',struct('units',{{'deg/s','deg/s','m/s'}},'names',{{'azim','elev','rad'}})),...
      'acc',struct(...
        'label','acceleration',...
        'size',3,...
        'xyz',struct('units',{{'m/s^2',  'm/s^2',  'm/s^2'}},'names',{{'x',   'y',   'z'}}),...
        'sph',struct('units',{{'deg/s^2','deg/s^2','m/s^2'}},'names',{{'azim','elev','rad'}})),...
      'pos_cor',struct(...
        'label','pos. correlation',...
        'size',6,...
        'xyz',struct('units',{{'m^2',  'm^2',  'm^2', 'm^2', 'm^2',  'm^2'}},  'names',{{'xx','yy','zz','xy','xz','yz'}}),...
        'sph',struct('units',{{'deg^2','deg^2','m^2','deg^2','deg*m','deg*m'}},'names',{{'xx','yy','zz','xy','xz','yz'}})),...
      'clk',struct(...
        'label','clock correction',...
        'size',1,...
        'xyz',struct('units',{{'s'}},  'names',{{'t'}}),...
        'sph',struct('units',{{'1/s'}},'names',{{'f'}})),...
      'clk_cor',struct(...
        'label','clk. corr. correlation',...
        'size',4,....
        'xyz',struct('units',{{'s^2','m.s',  'm.s',  'm.s'}},'names',{{'tt','xt','yt','zt'}}),...
        'sph',struct('units',{{'s^2','deg.s','deg.s','m.s'}},'names',{{'tt','xt','yt','zt'}}))...
    );
    %default value of parameters
    %NOTICE: needs updated when adding a new parameter
    parameter_list={...
      'satname',  'unknown',@(i) ischar(i);...
      'frame',    'crs',    @(i) orbit.isframe(i);...
      'geodatum' ,'grs80',  @(i) ischar(i);...
      'sp3id',    'Xnn',    @(i) ischar(i);...
      'data_dir', fullfile(getenv('HOME'),'data'),@(i) ischar(i);...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTICE: needs updated when adding a new data type (if relevant)
    compatible_parameter_list={'satname','frame'};

  end
  %read only
  %NOTICE: needs updated when adding a new data type
  %NOTICE: needs updated when adding a new parameter
  properties(SetAccess=private)
    %parameters
    satname
    frame
    geodatum
    sp3id
    data_dir
    %data types
    pos
    vel
    acc
    pos_cor
    clk
    clk_cor
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    localframei
  end
  %calculated only when asked for
  properties(Dependent)
    time
    localframe
  end
  methods(Static)
    %interface methods to object constants
    function out=data_types
      out=fieldnames(orbit.data_type_list);
    end
    function out=parameters(i,method)
      persistent v parameter_names
      if isempty(v)
        v=varargs(orbit.parameter_list);
        parameter_names=v.Parameters;
      end
      if ~exist('i','var') || isempty(i)
        if ~exist('method','var') || isempty(method)
          out=parameter_names(:);
        else
          switch method
          case 'obj'
            out=v;
          otherwise
            out=v.(method);
          end
        end
      else
        if ~exist('method','var') || isempty(method)
          method='name';
        end
        if strcmp(method,'name') && isnumeric(i)
          out=parameter_names{i};
        else
          switch method
          case varargs.template_fields
            out=v.get(i).(method);
          otherwise
            out=v.(method);
          end
        end
      end
    end    %translation methods
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
        otherwise
          out='';
        end
      else
        out='';
      end
    end
    function out=isframe(in)
      out=~isempty(orbit.translateframe(in));
    end
    function out=translatesat(in)
      %define orbit types
      orbit_types={'l1b','l2','kin'};
      %assume there's no orbit type
      orbit_type_now='';
      %search for an orbit type specification
      for i=1:numel(orbit_types)
        if strfind(in,['-',orbit_types{i}])
          %save orbit type
          orbit_type_now=orbit_types{i};
          %remove it from input
          in=strrep(in,['-',orbit_type_now],'');
          %stop searching for more orbit types
          break
        end
      end
      %search for satellite name
      switch lower(in)
        case {'champ','ch'}
          out='ch';
        case {'grace-a','gracea','grace a','ga'}
          out='ga';
        case {'grace-b','graceb','grace b','gb'}
          out='gb';
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
          error([mfilename,': cannot handle satellite ''',in,'''.'])
      end
      %append orbit type, if given
      if ~isempty(orbit_type_now)
        out=[out,'-',orbit_type_now];
      end
    end
    function out=translatesp3id(in)
      switch lower(in)
        case 'ch'
          out='TODO!';
        case 'ga'
          out='TODO!';
        case 'gb'
          out='TODO!';
        case 'sa'
          out='L47';
        case 'sb'
          out='L48';
        case 'sc'
          out='L49';
        case 'go'
          out='TODO!';            
        otherwise
          out='unknown';
      end
    end
    %data source definitions
    function out=nrtdm_product(in)
      switch orbit.translatesat(in)
        case 'ch'
          out='CH_Basic/Orbit_CH-OG-3-RSO';
        case 'ga'
          out='GA_Basic/Orbit_NAVSOL';
        case 'gb'
          out='GB_Basic/Orbit_NAVSOL';
        case 'sa-l1b'
          out='SA_Basic/Orbit_L1B';
        case 'sb-l1b'
          out='SB_Basic/Orbit_L1B';
        case 'sc-l1b'
          out='SC_Basic/Orbit_L1B';
        case 'sa-kin'
          out='SA_Basic/Orbit_KIN';
        case 'sb-kin'
          out='SB_Basic/Orbit_KIN';
        case 'sc-kin'
          out='SC_Basic/Orbit_KIN';
        case 'sa-l2'
          out='SA_Basic/Orbit_L2';
        case 'sb-l2'
          out='SB_Basic/Orbit_L2';
        case 'sc-l2'
          out='SC_Basic/Orbit_L2';
        otherwise
          error([mfilenane,': unknown NRTDM product for satellite ''',in,''', debug needed!'])
      end
    end
    %the function <format>_filename below define the filenames given sat, date and dir
    %for the data of different formats
    function [filename,dirname]=aiub_filename(satname,start,data_dir)
      if ~exist('data_dir','var') || isempty(data_dir)
        data_dir=orbit.parameters('data_dir','value');
      end
      switch orbit.translatesat(satname)
        case 'sa-kin'
          prefix='SWMA';
        case 'sb-kin'
          prefix='SWMB';
        case 'sc-kin'
          prefix='SWMC';
        otherwise
          error([mfilenane,': unknown AIUB orbit for satellite ''',in,''', debug needed!'])
      end
      doy=simpletimeseries.FromDateTime(start,'yeardoysec');
      filename=[prefix,num2str(doy(1)-round(doy(1)/1e3)*1e3),num2str(doy(2),'%03d'),'_S20.KIN.gz'];
      dirname=fullfile(data_dir,'gswarm','aiub','orbit',num2str(year(start)));
    end
    function [filename,dirname]=ifg_filename(satname,start,data_dir)
      if ~exist('data_dir','var') || isempty(data_dir)
        data_dir=orbit.parameters('data_dir','value');
      end
      switch orbit.translatesat(satname)
        case 'sa-kin'
          prefix='SwarmA-kinematicOrbit';
        case 'sb-kin'
          prefix='SwarmB-kinematicOrbit';
        case 'sc-kin'
          prefix='SwarmC-kinematicOrbit';
        otherwise
          error([mfilenane,': unknown IfG orbit for satellite ''',in,''', debug needed!'])
      end
      filename=[prefix,'-',num2str(year(start)),'-',num2str(month(start),'%02d'),'-',num2str(day(start),'%02d'),'.tar.gz'];
      dirname=fullfile(data_dir,'gswarm','ifg','orbit','ascii',num2str(year(start)));
    end
    function [filename,dirname]=tudelft_filename(satname,start,data_dir)
      if ~exist('data_dir','var') || isempty(data_dir)
        data_dir=orbit.parameters('data_dir','value');
      end
      switch orbit.translatesat(satname)
        case 'sa-kin'
          prefix='SWARMA';
        case 'sb-kin'
          prefix='SWARMB';
        case 'sc-kin'
          prefix='SWARMC';
        otherwise
          error([mfilenane,': unknown TU Delft orbit for satellite ''',in,''', debug needed!'])
      end
      doy=simpletimeseries.FromDateTime(start,'yeardoysec');
      filename=[prefix,'.',num2str(doy(1)-round(doy(1)/1e3)*1e3),'.',num2str(doy(2),'%03d'),'_KIPP.sigma.gz'];
      dirname=fullfile(data_dir,'gswarm','tudelft','orbit',num2str(year(start)));
    end
    function [filename,dirname]=sp3xcom_filename(satname,start,data_dir)
      if ~exist('data_dir','var') || isempty(data_dir)
        data_dir=orbit.parameters('data_dir','value');
      end
      switch orbit.translatesat(satname)
        case 'sa-l2'
          prefix='SW_OPER_SP3ACOM_2__';
        case 'sb-l2'
          prefix='SW_OPER_SP3BCOM_2__';
        case 'sc-l2'
          prefix='SW_OPER_SP3CCOM_2__';
        otherwise
          error([mfilenane,': unknown SP3xCOM orbit for satellite ''',in,''', debug needed!'])
      end
      s=num2str(60-sum(simpletimeseries.leap_seconds<start));
      filename=[prefix,...
        datestr(start-days(1),'yyyymmdd'),'T2359',s,'_',...
        datestr(start,        'yyyymmdd'),'T2359',s,'_'....
        '0101.DBL'];
      dirname=fullfile(data_dir,'swarm','dissemination');
    end
    %wrapper for the <format>_filename routines, transparent for different <format>s
    function out=filename(format,satname,start,varargin)
      p=inputParser;
      p.addRequired( 'satname',     @(i) ischar(i));
      p.addRequired( 'format',      @(i) ischar(i));
      p.addRequired( 'start',       @(i) isdatetime(i) && isscalar(i));
      p.addParameter('data_dir',...
        orbit.parameters('data_dir','value'),...
        orbit.parameters('data_dir','validation');
      p.parse(format,satname,start,varargin{:});
      %picking interface routine
      interface=str2func(['orbit.',format,'_filename']);
      %call interface routines
      [filename,dirname]=interface(satname,start,p.Results.data_dir);
      out=fullfile(dirname,filename);
    end
    %loads data from one single ASCII file, the format can be given or it
    %is discovered from the header. Also handles compressed files.
    %TODO: most of the zip-handling functionalityhas been implemented in simpletimeseries as is probably duplicate here.
    function obj=load_ascii(filename,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'filename',                          @(i) ischar(i));
      p.addParameter('asciiformat','',                    @(i) ischar(i));
      p.addParameter('tmpdir',     ['/tmp/orbit.load_ascii-',str.rand(8)], @(i) ischar(i));
      % parse it
      p.parse(filename,varargin{:});
      %trivial call
      if isempty(dir(filename))
        obj=[];
        disp(['WARNING: cannot find file ''',filename,''', skipping it.'])
        return
      end
      % unzip if needed
      [~,f,e]=fileparts(filename);
      %prepend tar if there
      if strcmp(f(end-3:end),'.tar')
        e=['.tar',e];
      end
      try
        switch lower(e)
          case {'.z','.zip'}
            unzipped_filename=unzip(filename,p.Results.tmpdir);
          case {'.tgz','.tar.gz','.tar'}
            unzipped_filename=untar(filename,p.Results.tmpdir);
            %get rid of PaxHeaders
            unzipped_filename(~cellfun(@isempty,strfind(unzipped_filename,'PaxHeaders')))=[];
          case {'.gz','.gzip'}
            unzipped_filename=gunzip(filename,p.Results.tmpdir);
          otherwise
            unzipped_filename={};
        end
      catch
        %if the zip file is corrupted, assume data file is missing
        disp(['WARNING: error extracting archive ''',filename,'''.'])
        obj=[];
        return
      end
      %handle zipped files
      if ~isempty(unzipped_filename)
        %some sanity
        if ~iscell(unzipped_filename)
          error([mfilename,': expecting variable ''unzipped_filename'' to be a cellstr, not a ''',...
            class(unzipped_filename),'''.'])
        end
        if numel(unzipped_filename)~=1
          error([mfilename,': expecting zip archive ''',filename,''' to contain one file only, not ',...
            num2str(numel(unzipped_filename)),':',10,strjoin(unzipped_filename,'\n')])
        end
        %recursive call
        obj=orbit.load_ascii(unzipped_filename{1},varargin{:});
        %clean up
        rmdir(p.Results.tmpdir,'s')
        %and we're done
        return
      end
      % retrieve the ascii format from the first line of the file, if not given in input arguments.
      if isempty(p.Results.asciiformat)
        %read first line of file
        [fid,errmsg] = fopen(filename);
        if fid<0
          error([mfilename,': ',errmsg])
        end
        hline = fgetl(fid);
        fclose(fid);
        %assign format ID
        if ~isempty(strfind(hline,'ITSG'))
          formatID='ifg';
        elseif ~isempty(strfind(hline,'LEOPOD')) || ~isempty(strfind(hline,'AIUB'))
          formatID='aiub';
        elseif ~isempty(strfind(hline,'#c'))
          formatID='sp3';
        else
          formatID='numeric';
        end
      else
        %propagate
        formatID=p.Results.asciiformat;
      end
      %branch on format
      switch lower(formatID)
      case 'ifg'
        [t,p,pc,header] = read_ifg(filename);
        if ~isempty(t)
          obj=orbit(simpletimeseries.ToDateTime(t,'modifiedjuliandate'),...
            'pos',      p,...
            'pos_cor',  pc,...
            'format',   header.timeformat,...
            'satname',  header.satname,...
            'sp3id',    orbit.translatesp3id(header.satname),...
            'frame',    header.frame,...
            'geodatum', header.geodatum,...
            'descriptor',lower(formatID)...
          );
        else
          obj=[];
        end
      case 'aiub'
        [t,p,pc,m,header] = read_aiub(filename);
        if ~isempty(t)
          obj=orbit(simpletimeseries.ToDateTime(t,'gpsweeksecond'),...
            'pos',      p,...
            'pos_cor',  pc,...
            'mask',     m,...
            'format',   header.timeformat,...
            'satname',  header.satname,...
            'sp3id',    header.sp3id,...
            'frame',    header.frame,...
            'geodatum', header.geodatum,...
            'descriptor',lower(formatID)...
          );
        else
          obj=[];
        end
      case {'numeric','tudelft'}
        [t,p,pc,c,cc,header] = read_numeric(filename);
        if ~isempty(t)
          obj=orbit(simpletimeseries.ToDateTime(t,'datevector'),...
            'pos',      p,...
            'pos_cor',  pc,...
            'clk',      c,...
            'clk_cor',  cc,...
            'format',   header.timeformat,...
            'satname',  header.satname,...
            'sp3id',    orbit.translatesp3id(header.satname),...
            'frame',    header.frame,...
            'geodatum', header.geodatum,...
            'descriptor',lower(formatID)...
          );
        else
          obj=[];
        end
      case {'sp3','sp3xcom'}
        [t,p,v,header] = read_sp3c(filename);
        if ~isempty(t)
          obj=orbit(simpletimeseries.ToDateTime(t,'datevector'),...
            'pos',      p,...
            'vel',      v,...
            'format',   header.timeformat,...
            'satname',  header.satname,...
            'sp3id',    header.sp3id,...
            'frame',    header.frame,...
            'geodatum', header.geodatum,...
            'descriptor',lower(formatID)...
          );
        else
          obj=[];
        end
      otherwise
        error([mfilename,': unknown format ''',format,'''.'])
      end
    end
    %reads data in any format, over any time period
    function obj=import(format,satname,start,stop,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'format',             @(i) ischar(i));
      p.addRequired( 'satname',            @(i) ischar(i));
      p.addRequired( 'start',              @(i) isdatetime(i) && isscalar(i));
      p.addRequired( 'stop',               @(i) isdatetime(i) && isscalar(i));
      p.addParameter('cut24h',  false,     @(i) islogical(i));
      p.addParameter('resample',seconds(0),@(i) isduration(i));
      p.addParameter('only_convert_to_mat', false,     @(i) islogical(i));
      % parse it
      p.parse(format,satname,start,stop,varargin{:});
      %clean
      varargin=cells.vararginclean(varargin,{'cut24h','resample','only_convert_to_mat'});
      % branch on format
      switch lower(format)
      case 'nrtdm'
        % retrieve product name
        product_name=orbit.nrtdm_product(satname);
        % retrieve data
        p=nrtdm(product_name,start,stop,varargin{:});
        % initialize data according to its dimension
        switch p.metadata.dimension
          case 3
            obj=orbit(p.ts.t,...
              'pos',p.ts.y,'pos_units',p.metadata.units,...
              'frame','ecf',...
              varargin{:}...
            );
          case 6
            obj=orbit(p.ts.t,...
              'pos',p.ts.y(:,1:3),'pos_units',p.metadata.units(1:3),...
              'vel',p.ts.y(:,4:6),'vel_units',p.metadata.units(4:6),...
              'frame','ecf',...
              varargin{:}...
            );
          case 9
            obj=orbit(p.ts.t,...
              'pos',p.ts.y(:,1:3),'pos_units',p.metadata.units(1:3),...
              'vel',p.ts.y(:,4:6),'vel_units',p.metadata.units(4:6),...
              'acc',p.ts.y(:,7:9),'acc_units',p.metadata.units(7:9),...
              'frame','ecf',...
              varargin{:}...
            );
          otherwise
            error([mfilename,': cannot handle data of size ',num2str(p.metadata.dimension).'.'])
        end
      otherwise
        % build required file list (rounding start/stop to the start of the day
        day_list=simpletimeseries.list(dateshift(start,'start','day'),dateshift(stop,'start','day'),days(1));
        file_list=cell(size(day_list));
        for i=1:numel(day_list)
          file_list{i}=orbit.filename(format,satname,day_list(i),varargin{:});
        end
        % load data
        first=true;
        for i=1:numel(file_list)
          %check if mat file is already available
          [d,f]=fileparts(file_list{i});
          mat_file=fullfile(d,[f,'.mat']);
          if isempty(dir(mat_file))
            % load ascii data
            obj_now=orbit.load_ascii(file_list{i},'asciiformat',format);
            if ~isempty(obj_now)
              if p.Results.resample>seconds(0)
                obj_now=obj_now.op('resample',p.Results.resample);
              end
              %save satname (those derived from the headers miss the '-kin' or '-l2' part)
              obj_now.satname=satname;
              %update sp3id if needed
              if strcmp(obj_now.sp3id,'unknown') || strcmp(obj_now.sp3id,orbit.parameters('sp3id','value'))
                obj_now.sp3id=orbit.translatesp3id(obj_now.satname);
              end
              % save it in mat format for next time
              save(mat_file,'obj_now');
            end
          else
            if ~p.Results.only_convert_to_mat
              % load mat data
              disp([mfilename,': loading file ',mat_file])
              S=load(mat_file);
              obj_now=S.obj_now;
            end
          end
          % if only converting orbits to mat, skip further operations
          if p.Results.only_convert_to_mat
            obj_now=[];
          end
          %skip if orbit file is missing
          if isempty(obj_now)
              continue
          end
          %cut to 24h if requested or if numerous files are being loaded
          %(otherwise appending doesn't work)
          if p.Results.cut24h || numel(file_list)>1
            obj_now=obj_now.op('trim',day_list(i),day_list(i)+hours(24)-seconds(1));
          end
          %skip if trimming removes all data
          if isempty(obj_now)
            continue
          end
          %create/append
          if first
            obj=obj_now;
            first=false;
          else
            % append remaining days
            obj=obj.op('append',obj_now);
          end
        end
        %nothing else to do if only converting orbits to mat
        if p.Results.only_convert_to_mat
          obj=[];
          return
        end
        %fill gaps
        obj=obj.op('trim',start,stop).op('resample');
      end
    end
    %general test for the current object
    function out=test_parameters(field,varargin)
      %basic parameters
      switch field
      case 'l'; out=100; return
      end
      %optional parameters
      switch numel(varargin)
      case 0
        l=orbit.test_parameters('l');
      case 1
        l=varargin{1};
      end
      %more parameters
      switch lower(field)
      case 'pos'
        out=randn(l,3);
      case 'vel'
        out=randn(l,3)+1;
      case 'acc'
        out=randn(l,3)+2;
      case 'step'
        out=10;
      case 'time'
        k=orbit.test_parameters('step');
        out=now+(1:k:k*l);
      case 'start'
        out=datetime(2015,2,1,23,0,0);
      case 'stop'
        if ~isduration(l)
          error([mfilename,': expecting input ''l'' to be of class ''duration'', not ''',class(l),'''.'])
        end
        out=orbit.test_parameters('start')+l;
      case 'duration'
        out=hours(3);  
      case 'satname'
        out='swarm-a-kin';
      case 'satname-l2'
        out='swarm-a-l2';
      case {'tudelft','aiub','ifg'}
        out=orbit.import(...
          field,...
          orbit.test_parameters('satname'),...
          orbit.test_parameters('start'),...
          orbit.test_parameters('stop',l));
      case 'sp3xcom'
        out=orbit.import(...
          'sp3xcom',...
          orbit.test_parameters('satname-l2'),...
          orbit.test_parameters('start'),...
          orbit.test_parameters('stop',l));
      otherwise
        error([mfilename,': unknown field ',field,'.'])
      end
    end
    function out=test(l)
    
      if ~exist('l','var') || isempty(l)
        l=1e4;
      end
      
      switch class(l)
      case 'cell'
        figure
        out=cell(size(l));
        for i=1:numel(l)
          out{i}=orbit.test(l{i}); hold on
        end
        plot_line_color
        legend(l)
      case 'char'
        switch lower(l)
        case 'formats'
          a=orbit.test({'tudelft','aiub','ifg','sp3xcom'});
        case 'rel'
          a=orbit.test_parameters('tudelft',orbit.test_parameters('duration'));
          b=orbit.test_parameters('ifg',    orbit.test_parameters('duration'));
          c=orbit.test_parameters('aiub',   orbit.test_parameters('duration'));
          d=orbit.test_parameters('sp3xcom',orbit.test_parameters('duration'));
          out={d.relative(a),...
               d.relative(b),...
               d.relative(c)...
          };
          out{1}=out{1}.op('descriptor','tudelft');
          out{2}=out{2}.op('descriptor','ifg');
          out{3}=out{3}.op('descriptor','aiub');
          if nargout==0
            figure
            for j=1:numel(out)
              for i=1:3
                subplot(3,numel(out),i+numel(out)*(j-1))
                out{j}.pos.plot('column',i)
              end
              out{j}.print
            end
          end
          
        case 'stats'
          a=orbit.test('rel');
          out=a.periodic_stats(orbit.test_parameters('duration')/10);
          if nargout==0
            figure
            for i=1:3
              subplot(3,1,i)
              out.mean.pos.plot('column',i)
            end
          end
        otherwise
          out=orbit.test_parameters(l,hours(8));
          if ~isa(out,'orbit')
            error([mfilename,': cannot handle test of type ''',l,'''.'])
          end
          out.pos.plot('columns',1,'line',{'o-'})
          out.print
        end
      case 'double'
        a=orbit(...
          orbit.test_parameters('time',l),...
          'pos',orbit.test_parameters('pos',l),...
          'vel',orbit.test_parameters('vel',l),...
          'acc',orbit.test_parameters('acc',l)...
        );
        figure
        subplot(3,1,1)
        a.pos.plot('columns',1)
        subplot(3,1,2)
        a.vel.plot('columns',1)
        subplot(3,1,3)
        a.acc.plot('columns',1)
      otherwise
        error([mfilename,': cannot handle input ''l'' of class ''',class(l),'''.'])
      end

    end 
  end
  methods
    %% constructor
    function obj=orbit(t,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('t',@(i) ~isscalar(i) || (isstruct(i) && isscalar(i))); %this can be vector char, double, datetime or scalar struct
      %parse the arguments with the names defined in orbit.data_type_list
      for j=1:numel(orbit.data_types)
        %shorter names
        dtn=orbit.data_types{j};
        dts=orbit.data_type_list.(dtn).size;
        %declare data types
        p.addParameter( dtn,           [],         (@(i) isnumeric(i) && size(i,2)==dts && size(i,1)>0));
        p.addParameter([dtn,'_units'], cell(1,dts),(@(i) iscellstr(i) && numel(i)==dts));
        p.addParameter([dtn,'_labels'],cell(1,dts),(@(i) iscellstr(i) && numel(i)==dts));
      end
      %create argument object v and declare parameters p
      [v,~,obj]=varargs.wrap('parser',p,'sinks',{obj},'sources',{orbit.parameters([],'obj')},'mandatory',{t},varargin{:});
      %clean varargin
      varargin=cells.vararginclean(varargin,p.Parameters);
      % retrieve each data type
      for j=1:numel(orbit.data_types)
        %shorter names
        data_type=orbit.data_types{j};
        %skip if this data type is empty
        if ~isempty(v.(data_type))
          %add new data type
          obj=obj.add_data_type(...
            t,...
            data_type,v.(data_type),...
            'units',  v.([data_type,'_units']),...
            'labels', v.([data_type,'_labels']),...
            varargin{:}...
          );
        end
      end
      %initialize internal records
      obj.localframei=[];
    end
    function obj=add_data_type(obj,t,data_type,data_value,varargin)
      %simplify things
      data_type=lower(data_type);
      %parse input
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('t'         ,@(i) ~isscalar(i)); %this can be char, double, datetime
      p.addRequired('data_type' ,@(i)    ischar(i)); 
      p.addRequired('data_value',@(i) ~isscalar(i)); 
      p.addParameter('units' ,cell(1,size(data_value,2)),@(i) iscell(i) && numel(i)==size(data_value,2))
      p.addParameter('labels',cell(1,size(data_value,2)),@(i) iscell(i) && numel(i)==size(data_value,2))
      % parse it
      p.parse(t,data_type,data_value,varargin{:});
      varargin=cells.vararginclean(varargin,p.Parameters);
      %sanity
      if ~isempty(obj.(data_type))
        error([mfilename,': data of type ''',data_type,''' has already been created. Use another method to append data.'])
      end
      %check if this is a valid field and retrieve default details
      valid_data_type=false;
      for i=1:numel(orbit.data_types)
        %shorter names
        dtn=orbit.data_types{i};
        if strcmp(data_type,dtn)
          valid_data_type=true;
          default_units=orbit.data_type_list.(dtn).xyz.units;
          default_names=orbit.data_type_list.(dtn).xyz.names;
          data_width=orbit.data_type_list.(dtn).size;
          default_label=orbit.data_type_list.(dtn).label;
          break
        end
      end
      if ~valid_data_type
        error([mfilename,': cannot handle data of type ''',data_type,'''.'])
      end
      %check data width
      if data_width ~= size(data_value,2)
        error([mfilename,': for data of type ''',data_type,...
          ''', expecting the nr of cols to be ',num2str(data_width),...
          ', not ',num2str(size(data_value,2)),'.'])
      end
      %propagate units (cannot write to parse objects)
      units=p.Results.units;
      %set default units
      for i=1:size(data_value,2)
        if isempty(units{i})
          units{i}=default_units{i};
        end
      end
      %propagate labels (cannot write to parse objects)
      labels=p.Results.labels;
      %set default lavels
      for i=1:size(data_value,2)
        if isempty(labels{i})
          labels{i}=[default_names{i},' ',default_label];
        end
      end
      %call superclass for this data type
      obj.(data_type)=simpletimeseries(...
        p.Results.t,p.Results.data_value,...
        'units',units,...
        'labels',labels,...
        varargin{:}...
      );
    end
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      pn=[orbit.parameters;more_parameters(:)];
      for i=1:numel(pn)
        if isprop(obj,pn{i}) && isprop(obj_in,pn{i})
          obj.(pn{i})=obj_in.(pn{i});
        end
      end
      %propagate parameters of all non-empty data types
      for j=1:numel(orbit.data_types)
        %shorter names
        data_type=orbit.data_types{j};
        %sanity
        if xor(isempty(obj.(data_type)),isempty(obj_in.(data_type)))
          error([mfilename,': error propagating metadata of type ',data_type,': it does not exist in both objects.'])
        end
        %skip if data type is empty
        if ~isempty(obj.(data_type))
          obj.(data_type)=obj.(data_type).copy_metadata(obj_in.(data_type));
        end
      end
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      warning off MATLAB:structOnObject
      out=varargs(...
        structs.filter(struct(obj),[orbit.parameters;more_parameters(:)])...
      ).varargin;
      warning on MATLAB:structOnObject
    end
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=20;
      end
      disp(' --- Parameters --- ')
      for i=1:numel(orbit.parameters)
        %shorter names
        p=orbit.parameters{i};
        disp([p,repmat(' ',1,tab-length(p)),' : ',str.show(obj.(p))])
      end
%       d_list=orbit.data_types;
      d_list={'pos'};
      for i=1:numel(d_list)
        %shorter names
        d=d_list{i};
        if ~isempty(obj.(d))
          disp([' --- ',d,' --- '])
          obj.(d).print
        end
      end
    end
    %% time property
    function t=get.time(obj)
      t=[];
      odt=orbit.data_types;
      for i=1:numel(odt)
        if ~isempty(obj.(odt{i}))
          t=obj.(odt{i}).t;
          break
        end
      end
      if isempty(t)
        error([mfilename,': all data types are empty.'])
      end
    end
    function obj=set.time(obj,t)
      odt=orbit.data_types;
      for i=1:numel(odt)
        if ~isempty(obj.(odt{i}))
          obj.(odt{i})=obj.(odt{i}).interp(t,'interp1_args',{'spline'});
        end
      end
      obj.check_time
    end
    function check_time(obj)
      t=[];
      f='';
      odt=orbit.data_types;
      for i=1:numel(odt)
        if ~isempty(obj.(odt{i}))
          if isempty(t)
            t=obj.(odt{i}).t;
            f=odt{i};
          else
            if any(t~=obj.(odt{i}).t)
              error([mfilename,': time domain discrepancy between fields ''',f,''' and ''',odt{i},'''.'])
            end
          end
        end
      end
    end
    %% satname property
    function out=get.satname(obj)
      out=obj.satname;
    end
    function obj=set.satname(obj,in)
      obj.satname=orbit.translatesat(in);
    end
    %% management
    function compatible(obj1,obj2)
      %This method checks if the objectives are referring to the same
      %type of data, i.e. the data length is not important.
      parameters=orbit.compatible_parameter_list;
      for i=1:numel(parameters)
        if ~isequal(obj1.(parameters{i}),obj2.(parameters{i}))
          error([mfilename,': discrepancy in parameter ',parameters{i},': ''',...
            obj1.(parameters{i}),''' ~= ''',obj2.(parameters{i}),'''.'])
        end 
      end
      %check that all data type as compatible as well
      odt=orbit.data_types;
      for i=1:numel(odt)
        if ~isempty(obj1.(odt{i})) && ~isempty(obj2.(odt{i}))
          obj1.(odt{i}).compatible(obj2.(odt{i}))
        end
      end
    end
    %object obj1 will have the time domain of obj2 (interpolated if needed)
    function [obj1,obj2]=consolidate(obj1,obj2)
      %compatibility check
      compatible(obj1,obj2)
      %consolidate all data types
      counter=0;
      odt=orbit.data_types;
      for i=1:numel(odt)
        if ~isempty(obj1.(odt{i})) && ~isempty(obj2.(odt{i}))
          [obj1.(odt{i}),obj2.(odt{i})]=obj1.(odt{i}).consolidate(obj2.(odt{i}));
          counter=counter+1;
        end
      end
      if counter==0
        error([mfilename,': there were no common fields in the input objects.'])
      end
    end
    %% operator
    % uses a method from a superclass over all non-empty data types
    function obj=op(obj,operation,varargin)
      %operation counter
      counter=0;
      odt=orbit.data_types;
      %handle operations between two orbits
      if numel(varargin)==1 && isa(varargin{1},'orbit')
        err_msg='there are no common fields in the input objects';
        %shorter names
        obj1=varargin{1};
        %No need for consolidation, that is done inside the operation
        for i=1:numel(odt)
          %operate over all non-empty data types
          if ~isempty(obj.(odt{i})) && ~isempty(obj1.(odt{i}))
            %operate
            obj.(odt{i})=obj.(odt{i}).(operation)(obj1.(odt{i}));
            counter=counter+1;
          end
        end
      else
        err_msg='there are no non-empty fields in the input object';
        for i=1:numel(odt)
          %operate over all non-empty data types
          if ~isempty(obj.(odt{i}))
            if isprop(obj.(odt{i}),operation)
              %sanity
              if numel(varargin)>1
                error([mfilename,': when propagating data to field ',operation,...
                  ' can only handle one input argument, not ',num2str(numel(varargin)),'.'])
              end
              %propagate
              obj.(odt{i}).(operation)=varargin{1};
            else
              %operate
              obj.(odt{i})=obj.(odt{i}).(operation)(varargin{:});
            end
            counter=counter+1;
          end
        end
      end
      if counter==0
        error([mfilename,': ',err_msg,'.'])
      end
    end
    %TODO: need to implement sum and subtraction that consider the
    %correlations
    %% reference frame (delayed constructor)
    function out=get.localframe(obj)
      obj=localframe_refresh_if_empty(obj);
      out=obj.localframei;
    end
    function obj=localframe_refresh_if_empty(obj)
      if isempty(obj.localframei)
        obj=localframe_refresh(obj);
      end
    end
    function obj=localframe_refresh(obj)
      % radial vector (parallel to the position vector)
      obj.localframei.ra = simpletimeseries(obj.time,...
        obj.pos.unit,...
        'units', {'','',''},...
        'labels',orbit.data_type_list.pos.xyz.names,...
        'descriptor','radial direction'...
      );
      % cross-track vector (perpendicular to the orbital plane)
      obj.localframei.ct = simpletimeseries(obj.time,...
        obj.pos.autocross.unit,...
        'units', {'','',''},...
        'labels',orbit.data_type_list.pos.xyz.names,...
        'descriptor','cross-track direction'...
      );
      % along-track vector (perpendicular to the other two)
      obj.localframei.at = simpletimeseries(obj.time,...
        obj.localframei.ra.cross(obj.localframei.ct).y,...
        'units', {'','',''},...
        'labels',orbit.data_type_list.pos.xyz.names,...
        'descriptor','along-track direction'...
      );
    end
    function plot_localframe(obj,ts)
      %consolidate
      [ts,obj.pos]=ts.consolidate(obj.pos);
      [ts,ra]=ts.consolidate(obj.localframe.ra);
      [ts,ct]=ts.consolidate(obj.localframe.ct);
      [ts,at]=ts.consolidate(obj.localframe.at);
      %plot
      idx=5:obj.pos.length-5;
      quiver3(obj.pos.y(idx,1),obj.pos.y(idx,2),obj.pos.y(idx,3),...
                   ra.y(idx,1),     ra.y(idx,2),     ra.y(idx,3)), hold on
      quiver3(obj.pos.y(idx,1),obj.pos.y(idx,2),obj.pos.y(idx,3),...
                   ct.y(idx,1),     ct.y(idx,2),     ct.y(idx,3)), hold on
      quiver3(obj.pos.y(idx,1),obj.pos.y(idx,2),obj.pos.y(idx,3),...
                   at.y(idx,1),     at.y(idx,2),     at.y(idx,3)), hold on
      quiver3(obj.pos.y(idx,1),obj.pos.y(idx,2),obj.pos.y(idx,3),...
                   ts.y(idx,1),     ts.y(idx,2),     ts.y(idx,3)), hold on
      legend('ra','ct','at','vector')
      axis square
    end
    %defined as: along-track, cross-track and radial
    function obj=tolocalframe(obj,obj_lf)
      [obj,obj_lf]=consolidate(obj,obj_lf);
      odt=orbit.data_types;
      for i=1:numel(odt)
        if ~isempty(obj.(odt{i})) && obj.(odt{i}).width==3
          a1=obj.(odt{i});
          obj.(odt{i})=obj.(odt{i}).project(...
            obj_lf.localframe.at,...
            obj_lf.localframe.ct,...
            obj_lf.localframe.ra...
          );
          obj.(odt{i}).labels={'along-track','cross-track','radial'};
          a2=obj.(odt{i});
          disp(['Difference RMS between magnitude of unprojected/projected ',odt{i},':',...
            num2str(rms( sum(a1.y(5:end-5,:).^2,2) - sum(a2.y(5:end-5,:).^2,2) ))])
        end
      end
    end
    %% relative motion (in the local frame)
    function rel=relative(obj1,obj2)
      %compute difference
      rel=obj1.op('minus',obj2);
    end
    function out=periodic_stats(obj,period,varargin)
      % separate time series into segments
      [ts,idx]=segmentedfreqseries.time_segmented(obj.time,period,seconds(0));
      %get stats for all available stat types
      odt=orbit.data_types;
      for j=1:numel(odt)
        if ~isempty(obj.(odt{j}))
          % initialize
          s.msg=[mfilename,': cutting into segments data of type ''',odt{j},'''.'];s.n=numel(ts);
          clear tmp
          % propagate segments
          for i=1:numel(ts)
            %compute statistics
            tmp(i)=simpledata(...
              ts{i},...
              obj.(odt{j}).y(idx{i}(1):idx{i}(2),:),...
              'mask',obj.(odt{j}).mask(idx{i}(1):idx{i}(2),:),...
              varargin{:}...
            ).stats('struct'); %#ok<*AGROW>
            %inform
            s=simpledata.progress(s,i);
          end
          %propagate
          s_list=fields(tmp);
          for i=1:numel(s_list)
            stats.(odt{j}).(s_list{i})=transpose(reshape(...
              [tmp.(s_list{i})],...
              numel(tmp(1).(s_list{i})),...
              numel(ts)...
            ));
          end
        end
      end
      %build time domain
      t=datetime([],[],[]);
      for i=1:numel(ts)
        t(i)=mean(ts{i});
      end
      %build argument list
      o_list=fields(stats);
      args=cell(1,2*numel(o_list));
      %propagate all fields in statistics
      s_list=fields(tmp);
      for i=1:numel(s_list)
        %build argument list for this statistic
        for j=1:numel(o_list)
          args{2*j-1}=o_list{j};
          if size(stats.(o_list{j}).(s_list{i}),2)==1;
            args{2*j  }=repmat(stats.(o_list{j}).(s_list{i}),1,orbit.data_type_list.(o_list{j}).size);
          elseif size(stats.(o_list{j}).(s_list{i}),2)==orbit.data_type_list.(o_list{j}).size
            args{2*j  }=stats.(o_list{j}).(s_list{i});
          else
            error([mfilename,': BUG TRAP: statistic with non-comformant number of columns. Debug needed!'])
          end
        end
        %build orbit object for this statistic
        out.(s_list{i})=orbit(t,args{:}).copy_metadata(obj);
      end

    end
    
  end
end

function [t,pos,pos_cor,header] = read_ifg(filename)
  %parameters
  formatSpec='%21.15f %15.4f %15.4f %15.4f %15.9f %15.9f %15.9f %15.9f %15.9f %15.9f';
  HeaderLines=2;
  %open the file
  disp([mfilename,': loading file ',filename])
  fid=fopen(filename);
  %read the data (robustly)
  try
    d=textscan(fid,formatSpec,'HeaderLines',HeaderLines,'Delimiter',' ','MultipleDelimsAsOne',true);
  catch
    error([mfilename,': could not understand the format of the ASCII orbit in file ''',filename,'''.'])
  end
  fclose(fid);
  %sanity
  if ~iscell(d)
    error([mfilename,': BUG TRAP: expecting variable ''d'' to be a cell array, not a ',class(d),'.'])
  end
  %propagate
  t=d{1};
  pos=[d{2:4}];
  pos_cor=[d{5:10}];
  %retrieve first line
  fid=fopen(filename);
  hline = fgetl(fid);
  fclose(fid);
  %find geodetic datum
  GDstr='Geodetic Datum: ';
  idx=strfind(hline,GDstr);
  if isempty(idx)
    disp(['WARNING: cannot find string ''',GDstr,''' in the header of file ',filename,'; setting default value.'])
    header.geodatum='unknown';
    header.satname=hline; %there is no dedicated field for the satellite name but it is (probably) mentioned in the first line of the header
  else
    header.geodatum=strtrim(hline(idx+length(GDstr):end));
    header.satname=strrep(strrep(strrep(hline,...
      [', ',GDstr,header.geodatum],''),...
      '# ITSG - ',''),...
      ' kinematic orbit','');
  end
  %set known parameters
  header.timeformat='gpstime';
  header.frame='trf';
end
function [t,pos,pos_cor,mask,header] = read_aiub(filename)
  formatSpec=' %4s %3s         %4f %8f  %14f %14f %14f %1s   %13f%13f%13f%13f%13f%13f';
  HeaderLines=6;
  %open the file
  disp([mfilename,': loading file ',filename])
  fid=fopen(filename);
  %read the data (robustly)
  try
    d=textscan(fid,formatSpec,'HeaderLines',HeaderLines,'Delimiter',' ','MultipleDelimsAsOne',true);
  catch
    error([mfilename,': could not determine the format of the ASCII orbit in file ''',filename,'''.'])
  end
  fclose(fid);
  %sanity
  if ~iscell(d)
    error([mfilename,': BUG TRAP: expecting variable ''d'' to be a cell array, not a ',class(d),'.'])
  end
  %make sure the data refers to the same satellite
  if ~all(strcmp(d{1},d{1}{1})) || ~all(strcmp(d{2},d{2}{1}))
    error([mfilename,': there are orbits of multiple satellites in file ''',filename,'''; this is not supported.'])
  end
  %propagate
  t=[d{3:4}];
  pos=[d{5:7}];
  pos_cor=[d{9:14}];
  mask=strcmp(d{8},'K');
  %retrieve third line
  fid=fopen(filename);
  for i=1:3
    hline = fgetl(fid);
  end
  fclose(fid);
  %find geodetic datum
  GDstr='GEODETIC DATUM: ';
  GDidx=strfind(hline,GDstr);
  EPstr='EPOCH: ';
  EPidx=strfind(hline,EPstr);
  if isempty(GDidx) || isempty(EPidx)
    disp(['WARNING: cannot find string ''',GDstr,''' or ''',EPstr,''' in the header of file ',filename,'; setting default value.'])
    header.geodatum='unknown';
  else
    header.geodatum=strtrim(hline(GDidx+length(GDstr):EPidx-1));
  end
  %header stuff
  header.satname=d{1}{1};
  header.sp3id=d{2}{1};
  %set known parameters
  header.timeformat='gpstime';
  header.frame='trf';
end
function [t,pos,pos_cor,clk,clk_cor,header] = read_numeric(filename)
  formatSpec='';
  HeaderLines=0;
  %open the file
  disp([mfilename,': loading file ',filename])
  fid=fopen(filename);
  %read the data (robustly)
  try
    d=textscan(fid,formatSpec,'HeaderLines',HeaderLines,'Delimiter',' ','MultipleDelimsAsOne',true);
  catch
    error([mfilename,': could not determine the format of the ASCII orbit in file ''',filename,'''.'])
  end
  fclose(fid);
  %sanity
  if ~iscell(d)
    error([mfilename,': BUG TRAP: expecting variable ''d'' to be a cell array, not a ',class(d),'.'])
  end
  switch numel(d)
  case 20
    %propagate
    t=[d{1:6}];
    pos=[d{7:9}]*1e3;
    pos_cor=[d{[11:13,15:16,17]}];
    clk=sqrt(d{10})/299792458;
    clk_cor=([sqrt(d{14}),d{[17,19,20]}]/299792458);
  otherwise
    error([mfilename,': cannot understand format of file ''',filename,'''.'])
  end
  %set known parameters
  header.timeformat='gpstime';
  header.frame='trf';
  %set also unknown but mandatory parameters
  header.satname='unknown';
  header.geodatum='unknown';
end
function [t,pos,vel,header] = read_sp3c(filename)

% function [t,pos,vel,header] = read_sp3c(filename)
%
% Reads SP3c orbit files
% 
% Parameters IN:
%   o  filename:   Filename of File to read
%
% Parameters OUT:
%   o  t:          nx6 datevec matrix (in the time system stated in the header of the file)
%   o  pos         nx3-Matrix with Position
%   o  vel         nx3-Matrix with Velocity (when available)
%   o header       structure with several details retrieved from the header
%
% ? 2010 IAPG (Delf Neubersch & Markus Heinze)
% Last changed: 2010-04-08
% Changes: 2010-04-08  MH  : Help-Text added
%          2010-04-08  MH  : Modified routine to handle float seconds correct
%
% Changes by Ales Bezdek:
% ab: 26/5/14 jezismarja, tady chybelo fclose
% ab: Pridal jsem tsek, protoze jinak to nesmyslne ztratilo aspon 2 platne cislice z case v sekundach
%     tsek jsou tedy sekundy s 8 platnymi cislicemi zkopirovane z SP3c souboru
% ab: a pridavam dale take rovnou vystup pro hodiny a minuty
% ab: spatne nacteni rychlosti:       Vel=str2num(line(5:46));
%
% Changes by Joao Encarnacao (2/3/2016):
% - removed input argument 'type' (detected automatically)
% - replaced output arguments 'n_hod' 'n_min' 'rok' 'mesic' 'den' ( hour min year month day) with datevec 't'
% - replaced output 'M' with outputs 't' 'pos' and 'vel'
% - added extensive sanity on the header section
% - added header structure
% - using fscanf to improve speed (about 2x faster)

  %open the file
  disp([mfilename,': loading file ',filename])
  fid=fopen(filename);

  %first reading of the data, to get header and count data
  n=0;i=0;
  while(~feof(fid))
    %read this line
    line=fgetl(fid);
    %increment line counter
    i=i+1;
    %retrieve header info
    if i<=22
      %quick format checking
      switch i
        case  1;    stop = ~strcmp(line(1:2),'#c') && ~strcmp(line(1:2),'#b');
        case  2;    stop = ~strcmp(line(1:2),'##');
        case  3:7;  stop = ~strcmp(line(1:2),'+ ');
        case  8:12; stop = ~strcmp(line(1:2),'++');
        case 13:14; stop = ~strcmp(line(1:2),'%c');
        case 15:16; stop = ~strcmp(line(1:2),'%f');
        case 17:18; stop = ~strcmp(line(1:2),'%i');
        case 19:22; stop = ~strcmp(line(1:2),'/*');
      end
      if stop
        error([mfilename,': error in header line nr ',num2str(n),'.'])
      end
      %retrieve detail header info
      switch i
        case 1
          header.version   = line(1:2);
          header.PV        = line(3);
          switch header.PV
            case 'P'
              header.vel_flag=false;
            case 'V'
              header.vel_flag=true;
            otherwise
              error([mfilename,': PV flag with unsupported value: ''',header.PV,'''.'])
          end
          header.year      = str2double(line( 4: 7));
          header.month     = str2double(line( 9:10));
          header.dom       = str2double(line(12:13));
          header.hour      = str2double(line(15:16));
          header.min       = str2double(line(18:19));
          header.sec       = str2double(line(21:31));
          header.nrepochs  = str2double(line(33:39));
          header.data_used = line(41:45);
          header.geodatum  = line(47:51);
          header.type      = line(53:55);
          header.agency    = line(57:60);
        case 2
          header.gpswk     = str2double(line( 4: 7));
          header.gpswksec  = str2double(line( 9:23));
          header.dtime     = str2double(line(25:38));
          header.mjd       = str2double(line(40:44));
          header.dayfrac   = str2double(line(46:60));
        case 3
          header.nrsats    = str2double(line( 5: 6));
          %TODO: handle multiple satellites
          if header.nrsats>1
            error([mfilename,': can only handle one satellite per SP3 file. Implementation needed!'])
          end
          header.sp3id     = line(10:12);
  %       case 4:7   %TODO: implement reading sp3ids when multiple satellites are given in one SP3 file
  %       case 8:12  %TODO: implement reading satellite accuracies (if ever needed)
        case 13
          header.timesystem = line(10:12);
          if strcmp(header.timesystem,'GPS')
            header.timeformat='gpstime';
          else
            error([mfilename,': cannot handle time system ''',header.timesystem,'''. Implementation needed!'])
          end
  %       case 15:16
  %       case 17:18
        case 19:22
          if ~isfield(header,'comment')
            header.comment=line(4:end);
          else
            header.comment=[header.comment,' ',line(4:end)];
          end
      end
    else
      %count the number of epochs
      if line(1)=='*'
        n=n+1;
      end
    end
  end
  %sanity
  if n==0
    error([mfilename,': could not find any valid epochs in file ''',filename,'''.'])
  end
%   %user feedback
%   disp([...
%     'File               : ',filename,10,...
%     'Version Symbol     : ',header.version,10,...
%     'Pos or Vel Flag    : ',header.PV,10,...
%     'Year Start         : ',num2str(header.year),10,...
%     'Month Start        : ',num2str(header.month),10,...
%     'Day of Month St    : ',num2str(header.dom),10,...
%     'Hour Start         : ',num2str(header.hour),10,...
%     'Minute Start       : ',num2str(header.min),10,...
%     'Second Start       : ',num2str(header.sec),10,...
%     'Data Used          : ',num2str(header.data_used),10,...
%     'Geodetic Datum     : ',header.geodatum,10,...
%     'Orbit type         : ',header.type,10,...
%     'Agency             : ',header.agency,10,...
%     'GPS Week           : ',num2str(header.gpswk),10,...
%     'Seconds of Week    : ',num2str(header.gpswksec),10,...
%     'Mod Jul Day St     : ',num2str(header.mjd),10,...
%     'Fractional Day     : ',num2str(header.dayfrac),10,...
%     'Nr epochs (header) : ',num2str(header.nrepochs),10,...
%     'Nr epochs (actual) : ',num2str(n)...
%   ]);
  %fixing nr of epochs
  if n~=header.nrepochs
    disp(['WARNING: the number of epochs stated in the header (',num2str(header.nrepochs),...
      ') does not correspond to the actual data (',num2str(n),').'])
    header.nrepochs=n;
  end

  %rewind
  frewind(fid)
  %skip header
  for i=1:22
    fgetl(fid);
  end
  %allocating data
  t=nan(n,6);
  pos=nan(n,3);
  if header.vel_flag
    vel=nan(n,3);
  end
  %read data
  for i=1:n
    %*  2014 02 01 00 00 00.00000000
    t(i,:)=fscanf(fid,'* %d %d %d %d %d %f\n',6);
    %PL47  2412.4847619 -2988.1083749 -5695.2157915 999999.999999
    pos(i,:)=fscanf(fid,'P%*3s %f %f %f %*f\n',5)*1e3;
    if header.vel_flag
      %VL47-37019.3690920 50917.4070091-42452.2025333      0.000000
      vel(i,:)=fscanf(fid,'V%*3s %14f %14f %14f %*f\n',5)*1e-1;
    end
  end
  fclose(fid);
  %add known/derived parameters
  header.satname   = header.sp3id;
  header.frame     = 'trf';
end
