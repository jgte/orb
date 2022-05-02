classdef grace
  properties(Constant)
    l1bdir_options={...
      '~/data/grace';...
      './data/grace';...
    };
    % data_name, version, ...
    l1b_data={...
      'KBR1B','03';...
      'SCA1B','03';...
      'AHK1B','02';...
      'GNV1B','02';...
      'MAS1B','02';...
      'THR1B','02';...
      'CLK1B','02';...
      'GPS1B','02';...
      'IHK1B','02';...
      'MAG1B','02';...
      'TIM1B','02';...
      'TNK1B','02';...
      'USO1B','02';...
      'VSL1B','02';...
     };
    %list of data fields
    %NOTICE: needs updated when adding a new data type
    data_type_list={'scaA','scaB'};
    %default value of parameters
    %NOTICE: needs updated when adding a new parameter
    parameter_list={...
      'verbose',   false,    @islogical;...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTICE: needs updated when adding a new data type (if relevant)
    compatible_parameter_list={};
  end
  %NOTICE: needs updated when adding a new parameter
  properties
    %parameters
    verbose
  end
  %NOTICE: needs updated when adding a new data type
  properties(SetAccess=private)
    %data types
    scaa,scab
    initialized
  end
  %calculated only when asked for
  properties(Dependent)
    time
  end
  methods(Static)
    %% directories
    function out=dir(type)
      switch type
        case 'l1b'
          for i=1:numel(grace.l1bdir_options)
            if file.exist(grace.l1bdir_options{i})
              out=grace.l1bdir_options{i};
              return
            end
          end
        %add more directories here
      end
    end
    %% interface methods to object constants
    function out=data_types
      out=grace.data_type_list;
    end
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(grace.parameter_list); end
      out=v.picker(varargin{:});
    end
    function out=datatype2product(data_type)
      out=[upper(data_type(1:3)),'1B'];
    end
    function out=datatype2satname(data_type)
      out=grace.translatesat(['grace-',data_type(4)]);
    end
    %% internally-consistent naming of satellites
    function out=translatesat(in)
      %search for satellite name
      switch lower(in)
        case {'grace-a','gracea','grace a','ga'}
          out='ga';
        case {'grace-b','graceb','grace b','gb'}
          out='gb';
        case {'grace-c','gracec','grace c','gc'}
          out='gc';
        case {'grace-d','graced','grace d','gd'}
          out='gd';
        otherwise
          error('grace:BadSat',['Cannot handle GRACE satname  ''',satname,'''.'])
      end
    end
    function out=translatesatname(in)
      %search for satellite name
      switch grace.translatesat(in)
        case 'ga'; out='GRACE-A';
        case 'gb'; out='GRACE-B';
        case 'gc'; out='GRACE-C';
        case 'gd'; out='GRACE-D';
        otherwise
          error('grace:BadSatName',['Cannot handle GRACE satname  ''',satname,'''.'])
      end
    end
    %% specific data name-handling methods
    function sat=grace_l1b_sat(satname)
      switch simpletimeseries.translatesat(satname)
        case 'ga'; sat='A';
        case 'gb'; sat='B';
        case 'gc'; sat='C';
        case 'gd'; sat='D';
        otherwise
          error([mfilenane,': cannot handle GRACE satname  ''',satname,'''.'])
      end
    end
    function satname=grace_l1b_satname(sat)
      switch grace.translatesat(sat)
        case 'A'; satname='ga';
        case 'B'; satname='gb';
        case 'C'; satname='gc';
        case 'D'; satname='gd';
        otherwise
          error([mfilenane,': cannot handle GRACE sat ''',sat,''', debug needed!'])
      end
    end
    function version=grace_l1b_version(product)
      %set hard-coded version
      idx=find(strcmp(grace.l1b_data(:,1),product),1);
      assert(~isempty(idx),['Cannot handle product ''',product,'''.']);
      version=grace.l1b_data{:,2};
      %check if the version is defined in the project.yaml file
      global PROJECT
      project_varname=[product,'_GRACE_RL'];
      if isfield(PROJECT,project_varname)
        assert(ischar(PROJECT.(project_varname)),...
          ['Project parameter ''',project_varname,''' must char, not ',...
          class(PROJECT.(project_varname))])
        version=PROJECT.(project_varname);
      end
    end
    %NOTICE: data_dir is the top-most data dir, without specifying the satellite, data, etc
    function filename=grace_l1b_filename(product,satname,start,varargin)
      v=varargs.wrap('sources',{{...
        'version',grace.grace_l1b_version(product), @ischar;...
      }},varargin{:});
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,{...
        'data_dir',grace.grace_l1b_dirname(start,v.version),@ischar;...
      }},varargin{:});
      sat=grace.grace_l1b_sat(satname);
      date=time.FromDateTime(start,'yyyy-MM-dd');
      filename=fullfile(v.data_dir,[product,'_',date,'_',sat,'_',v.version,'.dat']);
    end
    %passes data_dir to dirname unless it is non-existing, empty or '.' (in which case, it
    %builds the default directory structure of the GRACE data dir
    function dirname=grace_l1b_dirname(start,version,data_dir)
      if ~exist('data_dir','var') ...
          || isempty(data_dir) ...
          || strcmp(data_dir,'.') ...
          || strcmp(data_dir,simpletimeseries.parameters('value','data_dir'))
        year=time.FromDateTime(start,'yyyy');
        dirname=fullfile(grace.dir('l1b'),'L1B','JPL',['RL',version],year);
      else
        dirname=data_dir;
      end
    end
    %GRACE L1B file names must be: <product>_yyyy-mm-dd_<sat>_<version>.dat
    function [product,sat,date,version,dirname]=strings_from_grace_l1b_filename(filename)
      %split input
      [d,f]=fileparts(filename);
      %split name of file
      fp=strsplit(f,'_');
      product=fp{1};
      date=fp{2};
      sat=fp{3};
      version=fp{4};
      dirname=grace.grace_l1b_dirname(time.ToDateTime(date,'yyyy-MM-dd'),version,d);
    end
    function [product,satname,start,version,dirname]=details_from_grace_l1b_filename(filename)
      %get strings
      [product,sat,date,version,dirname]=strings_from_grace_l1b_filename(filename);
      %convert
      start=time.ToDateTime(date,'yyyy-MM-dd');
      satname=grace.grace_l1b_satname(sat);
    end
    %% loading data
    function obj=import_format(filename,format)
      %NOTICE: this function is called by simpletimeseries.import_format
      %make sure this is a filename
      assert(ischar(filename),['Input ''filename'' must be of class ''char'', not ''',...
        class(filename),'''; consider using the simpletimeseries.batch_import method.'])
      %branch on extension/format ID
      switch format
% --------------
% CSR formats
% --------------
      case '.resid'
        fid=file.open(filename);
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
        idx=cells.strfind(filename,coords);
        %sanity
        assert(sum(idx)==1,[mfilename,': the name for .resid files must include (one of) AC0X, AC0Y or AC0Z, not ''',...
          filename,'''.'])
        %building object
        obj=simpletimeseries(t,y,...
          'format','datetime',...
          'units',{'m^2','m^2','m^2'},...
          'labels', {coords{idx},[coords{idx},'D'],[coords{idx},'Q']},...
          'timesystem','gps',...
          'descriptor',['model response from file ',filename]...
         );
      case '.sigma'
        fid=file.open(filename);
        raw = textscan(fid,'%d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f','delimiter',' ','MultipleDelimsAsOne',1);
        fclose(fid);
        %building time domain
        t=datetime([double([raw{1:5}]),raw{6}]);
        %building data domain
        y=[raw{7:end}];
        %building object
        obj=simpletimeseries(t,y,...
          'format','datetime',...
          'units',{'m','m','m','s','m^2','m^2','m^2','s^2','m^2','m^2','ms','m^s','ms','ms'},...
          'labels', {'x','y','z','t','xx', 'yy', 'zz', 'tt', 'xy', 'xz', 'xt','yz', 'yt','zt'},...
          'timesystem','utc',...
          'descriptor',['kinematic orbit from file ',filename]...
         );
% --------------
% GRACE L1B formats
% --------------
      case 'ACC1B-asc'
        obj=load_ACC1B(filename);
      case 'AHK1B-asc'
        obj=load_AHK1B(filename);
      case 'SCA1B-asc'
        obj=load_SCA1B(filename);
      case {...
        'ACC1B','AHK1B','GNV1B','KBR1B','MAS1B','SCA1B','THR1B','CLK1B',...
        'GPS1B','IHK1B','MAG1B','TIM1B','TNK1B','USO1B','VSL1B'}
        [d,f]=fileparts(filename);
        %get particles from filename
        [product,sat,date,version,dirname]=...
          grace.strings_from_grace_l1b_filename(fullfile(d,f));
        %define output file
        o=fullfile(dirname,[f,'.asc']);
        %invoke L1B cat script
        com=['./cat-l1b.sh ',...
          strrep(date,'-',''),' ',product,' ',sat,' ',version,' JPL > ',o];
        disp(com)
        %make sure there's a directory for o
        if ~exist(dirname,'dir'); mkdir(dirname); end
        %issue com
        file.system(com,'stop_if_error',true,'cd',grace.dir('l1b'));
        %recursive call to retrieve the data
        obj=grace.import_format(o,[format,'-asc']);
        %NOTICE: the data_dir was changed above, so need to bail to avoid writing a duplicate mat file in the default data_dir
        return
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
          obj=[];
        else
          %building object
          obj=simpletimeseries(t,y,...
            'format','datetime',...
            'units',{'m/s^2','m/s^2','m/s^2'},...
            'labels', {'x','y','z'},...
            'timesystem','gps',...
            'descriptor',strjoin(header,'\n')...
           ).fill;
        end
        %flip model upside down if needed
        for j=1:size(simpletimeseries.csr_acc_mod_invert_periods,1)
          invert_idx=simpletimeseries.csr_acc_mod_invert_periods(j,1) <= t & ...
                     simpletimeseries.csr_acc_mod_invert_periods(j,2) >= t;
          if any(invert_idx)
            assert(all(invert_idx),['Expecting all data between ',datestr(t(1)),' and ',datestr(t(end)),' to be withing ',...
              'one single inverted period, as defined in ''simpletimeseries.csr_acc_mod_invert_periods''.'])
            obj=obj.scale(-1);
          end
        end
      case 'msodp-acc'
        error([mfilename,': implementation needed'])
      case 'mjd'
        %define data cols
        data_cols=2;
        %load the data
        raw=file.textscan(filename,'%f %f');
        %get the labels (from the header)
        labels={'altitude'};
        %build units
        units=cell(size(data_cols));
        units(:)={'m'};
        %sort data
        [t,idx]=sort(raw(:,1));
        %remove duplicate data
        dup=diff(t)==0;
        %conver to datetime
        t=time.ToDateTime(t(~dup),'modifiedjuliandate');
        %building object
        obj=simpletimeseries(t,raw(idx(~dup),data_cols),...
          'format','datetime',...
          'units',units,...
          'labels', labels,...
          'timesystem','gps',...
          'descriptor','GRACE altitude'...
        );
      otherwise
        error('grace:BadImportFormat',['Cannot handle files of type ''',format,'''.'])
      end
    end
    %% utils
    function obj=GRACEaltitude(varargin)
      p=machinery.inputParser;
      p.addParameter('datafile',file.resolve_home(fullfile('~','data','grace','altitude','GRACE.altitude.dat')));
      p.parse(varargin{:});
      obj=simpletimeseries.import(p.Results.datafile,...
        'format','mjd',...
        'cut24hrs',false...
      );
    end
    %% exporting data
    function export(obj,filename,filetype,varargin)
      %NOTICE: this function is a functional duplicate of simpletimeseries.export, but
      %        implements GRACE-specific formats
      %TODO: needs testing
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
        case 'ACC1B' %expecting sat_name to be 'GRACE A' or 'GRACE B'
          gps_zero_epoch=datetime('2000-01-01 12:00:00');
          %enforce requested header type/value
          switch lower(v.header)
          case 'default'
            dh=[...
'PRODUCER AGENCY               : UTexas',10,...
'PRODUCER INSTITUTION          : CSR',10,...
'FILE TYPE ipACC1BF            : 8',10,...
'FILE FORMAT 0=BINARY 1=ASCII  : 1',10,...
'NUMBER OF HEADER RECORDS      : 23',10,...
'SOFTWARE VERSION              : N/A',10,...
'SOFTWARE LINK TIME            : N/A',10,...
'REFERENCE DOCUMENTATION       : N/A',10,...
'SATELLITE NAME                : ',v.sat_name,10,...
'SENSOR NAME                   : ACC',10,...
'TIME EPOCH (GPS TIME)         : ',datestr(gps_zero_epoch,'yyyy-mm-dd HH:MM:SS.FFF'),10,...
'TIME FIRST OBS(SEC PAST EPOCH): ',num2str(time.datetime2gpssec(obj.start,gps_zero_epoch)),...
  ' (',datestr(obj.start,'yyyy-mm-dd HH:MM:SS.FFF'),')',10,...
'TIME LAST OBS(SEC PAST EPOCH) : ',num2str(time.datetime2gpssec(obj.stop,gps_zero_epoch)),...
  ' (',datestr(obj.stop,'yyyy-mm-dd HH:MM:SS.FFF'),')',10,...
'NUMBER OF DATA RECORDS        : ',num2str(obj.length),10,...
'PRODUCT CREATE START TIME(UTC): ',datestr(datetime('now')),' by ',getenv('USER'),10,...
'PRODUCT CREATE END TIME(UTC)  : N/A',10,...
'FILESIZE (BYTES)              : N/A',10,...
'FILENAME                      : ',filename,10,...
'PROCESS LEVEL (1A OR 1B)      : 1B',10,...
'INPUT FILE NAME               : N/A',10,...
'INPUT FILE TIME TAG (UTC)     : N/A',10,...
'INPUT FILE NAME               : N/A',10,...
'INPUT FILE TIME TAG (UTC)     : N/A',10,...
'END OF HEADER',10];
            fprintf(fid,'%s',dh);
          case 'none'
            %do nothing
          otherwise
            fprintf(fid,'%s',v.header);
          end
          %build time vectors
          time_str=time.datetime2gpssec(obj.t_masked,gps_zero_epoch);
          %build format string (there's a lot of zeros because most of the original data is now lost)
          fmt=['%d ',strrep(v.sat_name,'GRACE ',''),...
            repmat(' %21.15e',1,numel(v.columns)),...
            repmat(' 0.0000000000000000000',1,9-numel(v.columns)),'  00000000\n'];
          %build output data
          y=obj.y_masked([],v.columns);
          %put everything together
          o=[num2cell(transpose(time_str));num2cell(transpose(y))];
          %fprintf it
          fprintf(fid,fmt,o{:});
        case 'msodp' %expecting sat_name to be '1201 GRACEA' or '1202 GRACEB'
          %translate satellite name: ACCREAD.f is very picky with this stuff
          switch lower(str.rep(v.sat_name,'-','',' ','','_','','.',''))
          %                             I7X                 A20
          case 'gracea'; sat_name='1201    GRACEA';
          case 'graceb'; sat_name='1202    GRACEB';
          otherwise; error(['unrecognized sat_name value ''',v.sat_name,'''.'])
          end
          %need only valid data
          obj=obj.masked;
          unitfacor=0.1e4;
          %enforce requested header type/value
          switch lower(v.header)
          case 'default'
            dh=cell(13,1);i=0;
%FORMAT (A10,X,A12,X,A12,X,I4,X,I2,X,I2,X,I2,X,I2,X,A9,X,A16)
%                    A10X         A12X         A12X  I4XI2XI2XI2XI2X       A9X             A16
i=i+1;dh{i}= '%grace.acc version 1.1  revision 2.1 2016 02 18 09:36 CSR/UT    Rick Pastor     ';
%FORMAT (A10,X,I7,X,A20)
%                    A10X
i=i+1;dh{i}=['+satellite ',sat_name];
%FORMAT (A10,10(X,A3))
%                    A10X A3X A3X
i=i+1;dh{i}= '+data_____ tim acl';
i=i+1;dh{i}= '+reference gps ifx';
%FORMAT (A10,X,I4,X,I2,X,I2,X,I2,X,I2,X,F10.7)
%                    A10X                       I4XI2XI2XI2XI2    X                              F10.7
i=i+1;dh{i}=['+first____ ',datestr(obj.start,'yyyy mm dd HH MM'),' ',num2str(second(obj.start),'%10.7e')];
i=i+1;dh{i}=['+last_____ ',datestr(obj.stop, 'yyyy mm dd HH MM'),' ',num2str(second(obj.stop ),'%10.7e')];
%FORMAT (A10X,I6)
i=i+1;dh{i}=['+interval_ ',num2str(seconds(obj.step),'%6i')];
i=i+1;dh{i}=['+datarecor ',num2str(obj.nr_valid,'%6i')];
%FORMAT (A10,X,E7.1)
i=i+1;dh{i}=['+unitfacor ',num2str(unitfacor,'%7.1e')];
%FORMAT (A10,X,A121)
i=i+1;dh{i}= '+format___ (I4,1X,I3,1X,I5,1X,I7,3(1X,F18.15),1X,I8)';
%FORMAT (A10,X,E7.3)
i=i+1;dh{i}= '+scmass___ 500.000';
%FORMAT (A10,X,A40)
i=i+1;dh{i}=['+software_ http://github.com/jgte/orb, data exported on ',...
  datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' by Joao Encarnacao'];
i=i+1;dh{i}= '+eoh______';
            fprintf(fid,'%s',strjoin(dh,newline));
          case 'none'
            %do nothing
          otherwise
            fprintf(fid,'%s',v.header);
          end
          %build time vectors
          sod=time.sod(obj.t);
          sod_floor=floor(sod);
          sod_fraction=sod-sod_floor;
          time_str=[year(obj.t),floor(time.doy(obj.t)),sod_floor,sod_fraction];
          %build format string (translated from header)
          fmt='%4d %3d %5d %7d %18.15f %18.15f %18.15f        0\n';
          %build output data
          y=obj.y(:,v.columns)*unitfacor;
          %put everything together
          o=[num2cell(transpose(time_str));num2cell(transpose(y))];
          %fprintf it
          fprintf(fid,fmt,o{:});
          %eof
          fprintf(fid,'%s','%eof');
        case 'seconds'
          error('implementation missing')
        otherwise
          error(['Cannot handle exporting time series to files of type ''',filetype,'''.'])
        end
        fclose(fid);
      end
    end
    %% testing for the current object
    function out=test_parameters(field)
      switch lower(field)
      case 'start'
        out=datetime('2010-01-01');
      case 'stop'
        out=datetime('2010-01-03');
      otherwise
        error([mfilename,': unknown field ',field,'.'])
      end
    end
    function out=test(method)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      start=grace.test_parameters('start');
      stop= grace.test_parameters('stop' );
      test_list={'import'};
      switch(method)
        case 'all'
          for i=1:numel(test_list)
            out{i}=grace.test(test_list{i});
          end
        case 'import'
          out=grace(start,stop);
        otherwise
          error(['Cannot handle test method ''',method,'''.'])
      end
    end
  end
  methods
    %% constructor
    function obj=grace(start,stop,varargin)
      p=machinery.inputParser;
      p.addRequired('start'); %defines the 1st  day of the data (can be char, double or datetime)
      p.addRequired('stop' ); %defines the last day of the data (can be char, double or datetime)
      %declare parameters p
      [~,p]=varargs.wrap('parser',p,'sinks',{obj},'sources',{grace.parameters('obj')},...
        'mandatory',{start,stop},varargin{:});
      %get day list
      startlist=time.day_list(p.Results.start,p.Results.stop);
      %clean varargin
      varargin=cells.vararginclean(varargin,p.Parameters);
      % retrieve each data type
      for i=1:numel(grace.data_types)
        %shorter names
        data_type=grace.data_types{i};
        %derive product and satname
        product=grace.datatype2product(data_type);
        satname=grace.datatype2satname(data_type);
        %loop over days
        for j=1:numel(startlist)
          %load the data
          switch product
            case 'SCA1B'
              data_value=attitude.import('grace_l1b',satname,startlist(j));
            otherwise
              error(['Cannot handle product ''',product,'''.'])
          end
          %add/append data type
          obj=obj.add_data_type(data_type,data_value,varargin{:});
        end
      end
      %initialize internal records
      obj.initialized=true;
    end
    function obj=add_data_type(obj,data_type,data_value)
      %simplify things
      data_type=lower(data_type);
      %parse input
      p=machinery.inputParser;
      p.addRequired('data_type' ,@(i) ischar(i) && cells.isincluded(lower(grace.data_types),i));
      p.addRequired('data_value',@(i) machinery.isa(i,{'simpletimeseries','attitude','orbit'}));
      % parse it
      p.parse(data_type,data_value);
      %create or append
      if isempty(obj.(data_type))
        obj.(data_type)=data_value;
      else
        obj.(data_type)=obj.(data_type).append(data_value);
      end
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      warning off MATLAB:structOnObject
      out=structs.filter(struct(obj),[grace.parameters('list');more_parameters(:)]);
      warning on MATLAB:structOnObject
    end
    function out=varargin(obj,more_parameters)
      out=varargs(obj.metadata(more_parameters)).varargin;
    end
    %% info methods
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=20;
      end
      disp(' --- Parameters --- ')
      for i=1:numel(grace.parameters('list'))
        %shorter names
        p=grace.parameters('value',i);
        disp([p,repmat(' ',1,tab-length(p)),' : ',str.show(obj.(p))])
      end
      d_list=grace.data_types;
      for i=1:numel(d_list)
        %shorter names
        d=d_list{i};
        if ~isempty(obj.(d))
          disp([' --- ',d,' --- '])
          obj.(d).print
        end
      end
    end
    function msg(obj,in)
      if obj.verbose
        disp(in)
      end
    end
  end
end

function [out,fid]=parse_grace_L1B_headers(filename)
  %define header details to be retrieved
  header_anchors={...
    'TIME EPOCH (GPS TIME)         : ','gps_time';...
    'NUMBER OF DATA RECORDS        : ','nr_data';...
    'SATELLITE NAME                : ','sat_name';...
  };
  out=struct([]);
  %open the file
  fid = file.open(filename);
  %run through the header and retrieve relevant parameters
  l='';
  while ~contains(l ,'END OF HEADER')
    l=fgetl(fid);
    for i=1:size(header_anchors,1)
      if contains(l,header_anchors{i,1})
        out(1).(header_anchors{i,2})=strrep(l,header_anchors{i,1},'');
      end
    end
  end
  %if fid is in output, don't close it
  if nargout<=1
    fclose(fid);
  end
end
function obj=load_ACC1B(filename)
  %parse headers
  header=parse_grace_L1B_headers(filename);
  %load data
  raw=file.textscan(filename,'%f %s %f %f %f %f %f %f %f %f %f %f');
  %building time domain
  t=time.gpssec2datetime(raw(:,1),header.gps_time);
  %gather data domain
  y=raw(:,2:4);
  %skip empty data files
  if isempty(t) || isempty(y)
    disp([mfilename,': this file has no data  ',filename])
    obj=[];
  else
    %building object
    obj=simpletimeseries(t,y,...
      'format','datetime',...
      'units',{'m/s^2','m/s^2','m/s^2'},...
      'labels', {'x','y','z'},...
      'timesystem','gps',...
      'descriptor',[strtrim(header.sat_name),' accelerometer (ACC1B)']...
     ).fill;
  end
end
function obj=load_AHK1B(filename)
  %parse headers, file opened inside
  [header,fid]=parse_grace_L1B_headers(filename);
  %init data vector and counter
  raw=cell(ceil(str2double(header.nr_data)/7),22);
  %inits
  c=0;
  %loop until the end
  while ~feof(fid)
    l=fgetl(fid);
    %get only the lines with temperature data
    if contains(l,'00000111111111111111111111000000')
      c=c+1;
%        1       2  3  4         5                                  6                  7           8               9              10              11              12            13              14         15          16          17          18           19          20          21          22           23          24           25        26     27
%222177690  893506  G  A  00000100   00000111111111111111111111000000  10.15329807017275 4.917429447 9.723482464e-09 2.481763683e-09 1.039316011e-08 1.959635476e-08 5.4534528e-09 5.740073306e-09 33.2671051 55.12400055 26.62501335 14.92129993 -14.91238022 5.074500084 21.62319946 14.02499962 -14.03499985 50.17699814 -50.44900131  00000001  60282
      ll=strsplit(strtrim(l),' ');
      raw(c,:)=[{[ll{1},'.',ll{2}]},ll(7:end)];
    end
  end
  fclose(fid);
  %convert string array to double
  raw=str2double(raw);
  %building time domain
  t=time.gpssec2datetime(raw(:,1),header.gps_time);
%               1  7                 8           9               10              11              12              13            14              15         16          17          18           19          20          21          22          23           24          25          26          27
%               1  2                 3           4               5               6               7               8             9               10         11          12          13           14          15          16          17          18           19          20          21          22
%222177690.893506  10.15329807017275 4.917429447 9.723482464e-09 2.481763683e-09 1.039316011e-08 1.959635476e-08 5.4534528e-09 5.740073306e-09 33.2671051 55.12400055 26.62501335 14.92129993 -14.91238022 5.074500084 21.62319946 14.02499962 -14.03499985 50.17699814 -50.44900131  00000001  60282
  %NOTICE: from the original file format, the following columns are relevant
  % Column 1: GPS seconds past J2000
  % Column 2: fraction of seconds
  % Column 6: the lines with the relevant information are tagged with ‘00000111111111111111111111000000’
  % Column 15: ‘SU’ temperature
  % Column 16: ‘ICU’ temperature
  % Column 17: ‘core’ temperature
  % Column 21: ‘ADC’ temperature
  %gather data domain
  y=raw(:,[10,12,11,16]);
  %skip empty data files
  if isempty(t) || isempty(y)
    disp([mfilename,': this file has no data  ',filename])
    skip_save_mat=true; %#ok<NASGU>
    obj=[];
  else
    %building object
    obj=simpletimeseries(t,y,...
      'format','datetime',...
      'units',{'deg C','deg C','deg C','deg C'},...
      'labels', {'SU','core','ICU','ADC'},...
      'timesystem','gps',...
      'descriptor',[strtrim(header.sat_name),' temperature (AHK1B)']...
     ).resample;
  end
end
function [obj,err,flag]=load_SCA1B(filename)
  %parse headers
  header=parse_grace_L1B_headers(filename);
  %load data                  t  A  7  qs q1 q2 q3 r  00000000
  raw=file.textscan(filename,'%f %s %f %f %f %f %f %f %f');
  %building time domain
  t=time.gpssec2datetime(raw(:,1),header.gps_time);
  %gather data domain
  y=raw(:,3:6);
  err=raw(:,7);
  flag=raw(:,8);
  %skip empty data files
  if isempty(t) || isempty(y)
    disp([mfilename,': this file has no data  ',filename])
    obj=[];
  else
    %building object
    obj=simpletimeseries(t,y,...
      'format','datetime',...
      'units',{' ',' ',' ',' '},...
      'labels', {'qs','qi','qj','qk'},...
      'timesystem','gps',...
      'descriptor',[strtrim(header.sat_name),' attitude (SCA1B)']...
     ).fill;
  end
end

% NOTICE: this is Christian's code, received by email on 20/1/2022 13:48 (relevant to GRACE-FO but applicable to GRACE too)
% function [D] = convert_grace_fo_l1_data(files)
% % [D] = convert_grace_fo_l1_data(files)
% %
% % Convert GRACE-FO Level 1A files from text to Matlab format. Return the
% % data of the last file.
% 
% % for testing:
% % files = dir('/Volumes/Elements/Data/GRACE-FO/L1A/ACC/ACC1A_2019-01-01_C_04.txt');
% % files = dir('/Volumes/Elements/Data/GRACE-FO/L1A/AHK/AHK1A_2019-01-01_C_04.txt');
% % files = dir('/Volumes/Elements/Data/GRACE-FO/L1A/THR/THR1A_2019-01-01_C_04.txt');
% % files = dir('/Volumes/Elements/Data/GRACE-FO/L1B/CLK/CLK1B_2018-06-14_D_04.txt');
% % files = dir('/Volumes/Elements/Data/GRACE-FO/L1B/ACT/ACT1B_2018-06-09_C_04.txt');
% 
% N = length(files);
% for n = 1:N
%     fn = filename(files(n));
%     
%     type = files(n).name(1:3);
%     
%     fid = fopen(fn);
%     
%     % skip header
%     while ~feof(fid)
%         line = fgetl(fid);
%         if contains(line,'# End of YAML header')
%             break
%         end
%     end
%     
%     % reset structure
%     clear D
%     
%     % default buffer size
%     B = 1e6;
%     
%     % read data
%     
%     %=========================================
%     %  ACC (L1A)
%     %=========================================
%     
%     if strcmp(type,'ACC')% || strcmp(type,'ACT')
%         
%         D.t = zeros(B,2);
%         D.tref = zeros(B,1);
%         D.id = zeros(B,1);
%         D.qual_flag = zeros(B,8);
%         D.prod_flag = zeros(B,32);
%         
%         % the leading space in front of the format specifier is needed for
%         % concatenating the strings
%         variables = {'lin_accl_x', ' %e'
%                      'lin_accl_y', ' %e'
%                      'lin_accl_z', ' %e'
%                      'ang_accl_x', ' %e'
%                      'ang_accl_y', ' %e'
%                      'ang_accl_z', ' %e'
%                      'bias_vol',   ' %e'
%                      'vd',         ' %e'
%                      'x1_out',     ' %e'
%                      'x2_out',     ' %e'
%                      'x3_out',     ' %e'
%                      'y1_out',     ' %e'
%                      'y2_out',     ' %e'
%                      'z1_out',     ' %e'
%                      'tesu',       ' %e'
%                      'taicu',      ' %e'
%                      'tisu',       ' %e'
%                      'v15picu',    ' %e'
%                      'v15micu',    ' %e'
%                      'vr5picu',    ' %e'
%                      'tcicu',      ' %e'
%                      'v15psu',     ' %e'
%                      'v15msu',     ' %e'
%                      'v48psu',     ' %e'
%                      'v48msu',     ' %e'
%                      '',''
%                      'icu_blk_nr', ' %d'
%                      'Tenhz_count',' %d'
%                      'Mhz_count',  ' %d'
%                      '',''
%                      '',''
%                      '',''};
%         
%         % init buffers for variables
%         for v = 1:size(variables,1)
%             if ~isempty(variables{v,1})
%                 D.(variables{v,1}) = zeros(B,1);
%             end
%         end
%         
%         % epoch counter
%         e = 0;
%         
%         while ~feof(fid)
%             % if mod(e,1000)==0
%             %     disp(e)
%             % end
%             
%             line = fgetl(fid);
%             
%             e = e + 1;
%             
%             [data,~,~,ind_next] = sscanf(line,'%d %d %c %c %8c %32c',[1,6]);
%             
%             % "seconds past 12:00:00 noon of January 1, 2000 in OBC Time"
%             D.t(e,:) = [data(1) data(2)/1e6]; % integer and fractional receiver time (seconds and microseconds)
%             
%             D.tref(e) = data(3); % R = Receiver, OBC, or LRI time
%                                  % G = GPS time
%                             
%             D.id(e) = data(4); % C or D
%             
%             D.qual_flag(e,:) = flip(data(5:12) == char('1')); % false of "0", true if "1"
%                                                               % first character = last bit
%             %           - bit 0 = 0 - Receiver Time, 1 - Spacecraft elapsed time
%             %           - bit 1 = 0 - Timing pulse sync, 1 - No Timing pulse sync
%             %           - bit 2 = Not defined
%             %           - bit 3 = Not defined
%             %           - bit 4 = ACC Mode. 0 - Normal range mode, 1 - Large range mode
%             %           - bit 5 = IPU nav/timing packet received
%             %           - bit 6 = No OBC-to-Receiver time mapping
%             %           - bit 7 = No ICU block number available for GRACE-FO
%         
%             % this flag defines which data is present in the line:
%             prod_flag = flip(data(13:44) == char('1'));
%             D.prod_flag(e,:) = prod_flag; % false of "0", true if "1"
%                                           % first character = last bit
%             %           - bit 0 = lin_accl_x
%             %           - bit 1 = lin_accl_y
%             %           - bit 2 = lin_accl_z
%             %           - bit 3 = ang_accl_x
%             %           - bit 4 = ang_accl_y
%             %           - bit 5 = ang_accl_z
%             %           - bit 6 = bias_vol
%             %           - bit 7 = Vd
%             %           - bit 8 = x1_out
%             %           - bit 9 = x2_out
%             %           - bit 10 = x3_out
%             %           - bit 11 = y1_out
%             %           - bit 12 = y2_out
%             %           - bit 13 = z1_out
%             %           - bit 14 = tesu
%             %           - bit 15 = taicu
%             %           - bit 16 = tisu
%             %           - bit 17 = v15picu
%             %           - bit 18 = v15micu
%             %           - bit 19 = vr5picu
%             %           - bit 20 = tcicu
%             %           - bit 21 = v15psu
%             %           - bit 22 = v15msu
%             %           - bit 23 = v48psu
%             %           - bit 24 = v48msu
%             %           - bit 25 = Not defined
%             %           - bit 26 = ICU block number
%             %           - bit 27 = 10Hz clock count
%             %           - bit 28 = MHz clock count
%             %           - bits 29-31 = Not defined
%             
%             % create the format string for reading the data
%             fmt = strcat(variables{prod_flag,2});
%             
%             % read the data from the line
%             data = sscanf(line(ind_next:end),fmt,[1,sum(prod_flag)]);
%             
%             % initialize data counter
%             d = 0;
%                         
%             % copy data into structure
%             for v = 1:size(variables,1)
%                 if prod_flag(v)
%                     d = d + 1;
%                     D.(variables{v,1})(e) = data(d);
%                 end
%             end
%         end
%         
%         % truncate buffers for variables to actual size
%         D.t = D.t(1:e,:);
%         D.tref = D.tref(1:e,:);
%         D.id = D.id(1:e,:);
%         D.qual_flag = D.qual_flag(1:e,:);
%         D.prod_flag = D.prod_flag(1:e,:);
%         
%         for v = 1:size(variables,1)
%             if ~isempty(variables{v,1})
%                 D.(variables{v,1}) = D.(variables{v,1})(1:e);
%             end
%         end
%         
%         
%     %=========================================
%     %  AHK (L1A)
%     %=========================================
%     
%     elseif strcmp(type,'AHK')
%         
%         D.t = zeros(B,2);
%         D.tref = zeros(B,1);
%         D.id = zeros(B,1);
%         D.qual_flag = zeros(B,8);
%         D.prod_flag = zeros(B,32);
%         
%         % the leading space in front of the format specifier is needed for
%         % concatenating the strings
%         variables = {'TFEEU_IF',   ' %e'
%                      'TFEEU_REF',  ' %e'
%                      'TFEEU_X',    ' %e'
%                      'TFEEU_YZ',   ' %e'
%                      'analog_GND', ' %e'
%                      'Vp3p3',      ' %e'
%                      'Vp',         ' %e'
%                      'MES_Vd',     ' %e'
%                      'MES_DET_X1', ' %e'
%                      'MES_DET_X2', ' %e'
%                      'MES_DET_X3', ' %e'
%                      'MES_DET_Y1', ' %e'
%                      'MES_DET_Y2', ' %e'
%                      'MES_DET_Z1', ' %e'
%                      'TSU_Yp',     ' %e'
%                      'TICUN',      ' %e'
%                      'TSU_Ym',     ' %e'
%                      'TSU_Zp',     ' %e'
%                      'TSU_Zm',     ' %e'
%                      'Vp5',        ' %e'
%                      'TICUR',      ' %e'
%                      'Vp15',       ' %e'
%                      'Vm15',       ' %e'
%                      'Vp48',       ' %e'
%                      'Vm48',       ' %e'
%                      '',''
%                      'icu_blk_nr', ' %d'
%                      'PPS_source', ' %d'
%                      'Sync_Qual',  ' %d'
%                      'status_flag',' %32c'
%                      '',''
%                      '',''};
%         
%         % init buffers for variables
%         for v = 1:size(variables,1)
%             if ~isempty(variables{v,1})
%                 D.(variables{v,1}) = zeros(B,1);
%             end
%         end
%         
%         % this variable has 32 bits
%         D.status_flag = zeros(B,32);
%         
%         % epoch counter
%         e = 0;
%         
%         while ~feof(fid)
%             % if mod(e,1000)==0
%             %     disp(e)
%             % end
%             
%             line = fgetl(fid);
%             
%             e = e + 1;
%             
%             [data,~,~,ind_next] = sscanf(line,'%d %d %c %c %8c %32c',[1,6]);
%             
%             % "seconds past 12:00:00 noon of January 1, 2000 in OBC Time"
%             D.t(e,:) = [data(1) data(2)/1e6]; % integer and fractional receiver time (seconds and microseconds)
%             
%             D.tref(e) = data(3); % R = Receiver, OBC, or LRI time
%                                  % G = GPS time
%                             
%             D.id(e) = data(4); % C or D
%             
%             D.qual_flag(e,:) = flip(data(5:12) == char('1')); % false of "0", true if "1"
%                                                               % first character = last bit
%             %           - bit 0 = 0 - Receiver Time, 1 - Spacecraft elapsed time
%             %           - bit 1 = 0 - Timing pulse sync, 1 - No Timing pulse sync
%             %           - bit 2 = Not defined
%             %           - bit 3 = Not defined
%             %           - bit 4 = ACC Mode. 0 - Normal range mode, 1 - Large range mode
%             %           - bit 5 = IPU nav/timing packet received
%             %           - bit 6 = No OBC-to-Receiver time mapping
%             %           - bit 7 = No ICU block number available for GRACE-FO
%         
%             % this flag defines which data is present in the line:
%             prod_flag = flip(data(13:44) == char('1'));
%             D.prod_flag(e,:) = prod_flag; % false of "0", true if "1"
%                                           % first character = last bit
%             %           - bit 0 = TFEEU_IF
%             %           - bit 1 = TFEEU_REF
%             %           - bit 2 = TFEEU_X
%             %           - bit 3 = TFEEU_YZ
%             %           - bit 4 = analog_GND
%             %           - bit 5 = +3.3V
%             %           - bit 6 = Vp
%             %           - bit 7 = MES_Vd
%             %           - bit 8 = MES_DET_X1
%             %           - bit 9 = MES_DET_X2
%             %           - bit 10 = MES_DET_X3
%             %           - bit 11 = MES_DET_Y1
%             %           - bit 12 = MES_DET_Y2
%             %           - bit 13 = MES_DET_Z1
%             %           - bit 14 = TSU_Y+
%             %           - bit 15 = TICUN
%             %           - bit 16 = TSU_Y-
%             %           - bit 17 = TSU_Z+
%             %           - bit 18 = TSU_Z-
%             %           - bit 19 = +5V
%             %           - bit 20 = TICUR
%             %           - bit 21 = +15V
%             %           - bit 22 = -15V
%             %           - bit 23 = +48V
%             %           - bit 24 = -48V
%             %           - bit 25 = Not defined
%             %           - bit 26 = ICU block number
%             %           - bit 27 = PPS Source
%             %           - bit 28 = Sync Quality Index
%             %           - bit 29 = Status flag
%             %           - bits 30-31 = Not defined
%             
%             % PPS_source:
%             %           - 0 - IN SYNC WITH GPS
%             %           - 1 - OUT OF SYNC
%             %           - 2 - SYNC IN PROGRESS BAD PPS
%             %           - 3 - SYNC IN PROGRESS GOOD PPS
%             %           - 4 - SYNC TO GROUND UTC IN PROGRESS
%             %           - 5 - IN SYNC WITH GROUND UTC
%             % Sync_Qual:
%             %           - 0 - PPS RECEIVED TIME PACKET RECEIVED
%             %           - 1 - PPS NOT RECEIVED TIME PACKET RECEIVED
%             %           - 2 - PPS RECEIVED TIME PACKET NOT RECEIVED
%             %           - 3 - PPS NOT RECEIVED TIME PACKET NOT RECEIVED
%             
%             % statusflag:
%             %           - bit 0 = Test registers. 0 - Test ADC registers function is enabled, 1 - Test ADC registers function is disabled
%             %           - bit 1 = Status, Vp. 0 - Vp 10V, 1 - Vp 40V
%             %           - bit 2 = Status, Vd. 0 - Vd 5Vrms, 1 - Vp 1.25Vrms
%             %           - bit 3 = Flag, SCI_ACC_X. 0 - ADC is synchronised, 1 - Time-out is detected
%             %           - bit 4 = Flag, SCI_ACC_Y
%             %           - bit 5 = Flag, SCI_ACC_Z
%             %           - bit 6 = Flag, SCI_ACC_Phi
%             %           - bit 7 = Flag, SCI_ACC_Theta
%             %           - bit 8 = Flag, SCI_ACC_Psi
%             %           - bit 9 = Flag, SCI_VP
%             %           - bit 10 = Flag, SCI_MUX1
%             %           - bit 11 = Flag, SCI_MUX2
%             %           - bit 12 = Flag, SCI_MUX3
%             %           - bit 13 = Flag, test registers. Shows when the ADC registers check has proceeded. 0 - ADC bit registers are NOT tested, 1 - ADC registers are tested
%             %           - bit 14 = Flag, validity of registers. 0 - All the ADC registers are identical to the initial value, 1 - At least one of the ADC registers is corrupt
%             %           - bit 15 = GWT. 0 - The stops are connected, 1 - The stops are NOT connected
%             %           - bit 16 = Event Report Generation. 0 - Event Report Generation function is disabled, 1 - Event Report Generation function is enabled
%             %           - bit 17 = PPS Channel. 0 - PPS A, 1 - PPS B
%             %           - bit 18 = Uart Channel. 0 - Uart A, 1 - Uart B
%             %           - bits 19-31 = Not defined
%           
%             % create the format string for reading the data
%             fmt = strcat(variables{prod_flag,2});
%             
%             % read the data from the line
%             data = sscanf(line(ind_next:end),fmt,[1,sum(prod_flag)]);
%             
%             % initialize data counter
%             d = 0;
%                         
%             % copy data into structure
%             for v = 1:size(variables,1)
%                 if prod_flag(v) && strcmp(variables{v,1},'status_flag')
%                     d = d + 1;
%                     D.(variables{v,1})(e,:) = data(d:d+31);
%                     d = d + 31;
%                 elseif prod_flag(v) && ~strcmp(variables{v,1},'status_flag')
%                     d = d + 1;
%                     D.(variables{v,1})(e) = data(d);
%                 end
%             end
%         end
%         
%         % truncate buffers for variables to actual size
%         D.t = D.t(1:e,:);
%         D.tref = D.tref(1:e,:);
%         D.id = D.id(1:e,:);
%         D.qual_flag = D.qual_flag(1:e,:);
%         D.prod_flag = D.prod_flag(1:e,:);
%         
%         for v = 1:size(variables,1)
%             if ~isempty(variables{v,1})
%                 D.(variables{v,1}) = D.(variables{v,1})(1:e,:);
%             end
%         end
%         
%     %=========================================
%     %  THR (L1A)
%     %=========================================
%     
%     elseif strcmp(type,'THR')
%         
%         data = fscanf(fid,'%d %d %c %c %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %8c');
%         
%         data = reshape(data,54,length(data)/54)';
%         
%         % "Seconds past 12:00:00 noon of January 1, 2000 in GPS time"
%         D.t = [data(:,1) data(:,2)/1e6]; % integer and fractional receiver time (seconds and microseconds)
% 
%         D.tref = data(:,3); % R = Receiver, OBC, or LRI time
%                             % G = GPS time
% 
%         D.id = data(:,4); % C or D
%         
%         D.thrust_count_att_ctrl_1_1 = data(:,5);
%         D.thrust_count_att_ctrl_1_2 = data(:,6);
%         D.thrust_count_att_ctrl_1_3 = data(:,7);
%         D.thrust_count_att_ctrl_1_4 = data(:,8);
%         D.thrust_count_att_ctrl_1_5 = data(:,9);
%         D.thrust_count_att_ctrl_1_6 = data(:,10);
%         D.thrust_count_att_ctrl_2_1 = data(:,11);
%         D.thrust_count_att_ctrl_2_2 = data(:,12);
%         D.thrust_count_att_ctrl_2_3 = data(:,13);
%         D.thrust_count_att_ctrl_2_4 = data(:,14);
%         D.thrust_count_att_ctrl_2_5 = data(:,15);
%         D.thrust_count_att_ctrl_2_6 = data(:,16);
%         % D.thrust_count_undef_1 = data(:,17);
%         % D.thrust_count_undef_2 = data(:,18);
%         D.on_time_att_ctrl_1_1 = data(:,19);
%         D.on_time_att_ctrl_1_2 = data(:,20);
%         D.on_time_att_ctrl_1_3 = data(:,21);
%         D.on_time_att_ctrl_1_4 = data(:,22);
%         D.on_time_att_ctrl_1_5 = data(:,23);
%         D.on_time_att_ctrl_1_6 = data(:,24);
%         D.on_time_att_ctrl_2_1 = data(:,25);
%         D.on_time_att_ctrl_2_2 = data(:,26);
%         D.on_time_att_ctrl_2_3 = data(:,27);
%         D.on_time_att_ctrl_2_4 = data(:,28);
%         D.on_time_att_ctrl_2_5 = data(:,29);
%         D.on_time_att_ctrl_2_6 = data(:,30);
%         D.on_time_orb_ctrl_1 = data(:,31);
%         D.on_time_orb_ctrl_2 = data(:,32);
%         D.accum_dur_att_ctrl = data(:,33);
%         % D.accum_dur_undef_1 = data(:,34);
%         % D.accum_dur_undef_2 = data(:,35);
%         % D.accum_dur_undef_3 = data(:,36);
%         % D.accum_dur_undef_4 = data(:,37);
%         % D.accum_dur_undef_5 = data(:,38);
%         % D.accum_dur_undef_6 = data(:,39);
%         % D.accum_dur_undef_7 = data(:,40);
%         % D.accum_dur_undef_8 = data(:,41);
%         % D.accum_dur_undef_9 = data(:,42);
%         % D.accum_dur_undef_10 = data(:,43);
%         % D.accum_dur_undef_11 = data(:,44);
%         D.accum_dur_orb_ctrl = data(:,45);
%         % D.accum_dur_undef_12 = data(:,46);
%         
%         D.qualflg = flip(data(:,47:54) == char('1')); % false of "0", true if "1"
%                                                     % first character = last bit
%         %           - bit 0 = On-time not calculated
%         %           - bit 1 = Multiple unaccounted thrusts prior to current record
%         %           - bit 2 = 1, branch 1 is active
%         %           - bit 3 = 1, branch 2 is active
%         %           - bits 4-5 = Not defined
%         %           - bit 6 = No OBC-to-Receiver Time mapping
%         %           - bit 7 = No clock correction available
%         
%     %=========================================
%     %  CLK (L1B)
%     %=========================================
%     
%     elseif strcmp(type,'CLK')
%         
%         data = fscanf(fid,'%d %c %d %e %e %e %e %8c');
%         
%         data = reshape(data,15,length(data)/15)';
%         
%         % "Seconds past 12:00:00 noon of January 1, 2000 in receiver time"
%         D.t = data(:,1); % integer and fractional receiver time (seconds and microseconds)
% 
%         D.id = data(:,2); % C or D
%         
%         D.clock_id = data(:,3); % C or D
%         
%         % GPS Time = rcv_time + eps_time
%         D.eps_time = data(:,4);
%         D.eps_err = data(:,5);
%         D.eps_drift = data(:,6);
%         D.drift_err = data(:,7);
%         
%         D.qualflg = flip(data(:,8:15) == char('1')); % false of "0", true if "1"
%                                                     % first character = last bit
%         %           - bit 0 = Linear extrapolation of eps_time not valid AFTER rcv_time
%         %           - bit 1 = Linear extrapolation of eps_time not valid BEFORE rcv_time
%         %           - bit 2 = Smoothed data lies before center of smoothing for start midnight
%         %           - bit 3 = Smoothed data lies after center of smoothing for start midnight
%         %           - bit 4 = Smoothed data lies before center of smoothing for end midnight
%         %           - bit 5 = Smoothed data lies after center of smoothing for end midnight
%         %           - bits 6-7 = Not defined
%         
%     %=========================================
%     %  ACT (L1B)
%     %=========================================
%     
%     elseif strcmp(type,'ACT')
%         
%         data = fscanf(fid,'%d %c %e %e %e %e %e %e %e %e %e %8c');
%         
%         data = reshape(data,19,length(data)/19)';
%         
%         % "Seconds past 12:00:00 noon of January 1, 2000 in GPS time"
%         D.t = data(:,1); % integer and fractional receiver time (seconds and microseconds)
% 
%         D.id = data(:,2); % C or D
%                 
%         D.lin_accl_x = data(:,3);
%         D.lin_accl_y = data(:,4);
%         D.lin_accl_z = data(:,5);
%         D.ang_accl_x = data(:,6);
%         D.ang_accl_y = data(:,7);
%         D.ang_accl_z = data(:,8);
%         D.acl_x_res = data(:,9);
%         D.acl_y_res = data(:,10);
%         D.acl_z_res = data(:,11);
%         
%         D.qualflg = flip(data(:,12:19) == char('1')); % false of "0", true if "1"
%                                                       % first character = last bit
%         %           - bit 0 = Vp (proof mass voltage) out of nominal range
%         %           - bit 1 = Not defined
%         %           - bit 2 = At least one of acl_x_res, acl_y_res, and acl_z_res is > 10 microns/s^2
%         %           - bit 3 = Extrapolated CLK1B clock offset value used for data point > 5 s from center of CRN filter window
%         %           - bit 4 = Extrapolated CLK1B clock offset value used for data point <= 5 s from center of CRN filter window
%         %           - bit 5 = Interpolated data point (due to gap) exists > 15 s from center of CRN filter window
%         %           - bit 6 = Interpolated data point (due to gap) exists > 5 s but < 15 s from center of CRN filter window
%         %           - bit 7 = Interpolated data point (due to gap) exists < 5 s from center of CRN filter window
%           
%     else
%         error(['Unkown file type: ' files(n).name(1:3)])
%     end
%     
%     % rename before saving
%     eval([type ' = D;'])
%         
%     % save data in Matlab format
%     save([fn(1:end-3) 'mat'],type)
%     
%     % close file
%     fclose(fid);
% end


%this used to be part of grace.import_format; retired for now
%       case '.GraceAccCal'
%         fmt='';
%         if ~isempty(regexp(filename,'AC0[XYZ]\d?\.aak','once')) || ~isempty(regexp(filename,'AC0[XYZ]\d?\.accatt','once'))
%           % 2002 04 05 2002.4.4. 23.59.47.00000000 1498260002 0.2784215319157E-07
%           fmt='%d %d %d %s %s %d %f';
%           units={'m/s^2',''};
%           labels={str.clean(filename,{'file','grace','.'}),'Job ID','arc start'};
%           time_fh=@(raw) time.utc2gps(...
%             datetime(...
%               strcat(...
%                 strrep(cellfun(@(x) [x(1:end-1),' '],raw{4},'UniformOutput',false),'.','/'),...
%                 strrep(strrep(raw{5},'.00000000',''),'.',':')...
%               ),'InputFormat','yyyy/MM/dd HH:mm:ss'...
%             )...
%           );
%           data_fh=@(raw) [raw{7},double(raw{6})];
%           timesystem='gps';
%           sanity_check=@(raw) true;
%         end
%         if ~isempty(regexp(filename,'AC0[XYZ][QD]\d?\.aak','once')) || ~isempty(regexp(filename,'AC0[XYZ][QD]\d?\.accatt','once'))
%           % 2002 04 05 2002.4.4. 23.59.47.00000000 1498260002  0.1389481692269E-07 52368.99985
%           fmt='%d %d %d %s %s %d %f %f';
%           units={'m/s^2','','MJD days'};
%           labels={str.clean(filename,{'file','grace','.'}),'Job ID','t_0','arc start'};
%           time_fh=@(raw) time.utc2gps(...
%             datetime(...
%               strcat(...
%                 strrep(cellfun(@(x) [x(1:end-1),' '],raw{4},'UniformOutput',false),'.','/'),...
%                 strrep(strrep(raw{5},'.00000000',''),'.',':')...
%               ),'InputFormat','yyyy/MM/dd HH:mm:ss'...
%             )...
%           );
%           data_fh=@(raw) [raw{7},double(raw{6}),raw{8}];
%           timesystem='gps';
%           sanity_check=@(raw) true;
%         end
%         if ~isempty(regexp(filename,'AC0[XYZ]\d?\.estim','once')) || ~isempty(regexp(filename,'AC0[XYZ][DQ]\d?\.estim','once'))
%           % 1    2  3  4  5  6  7     8 9   10      11       12                    13                     14                    15           16
%           % 2002 04 05 04/05/02 52369 1 0.0 26400.0 1593715  3.774424464092000e-08 -3.585594302740665e-09 3.415865033817934e-08 2.82822E-09  71279987.
%           fmt='%d %d %d %d/%d/%d %f %d %f %f %d %f %f %f %f %f';
%           units={'m/s^2','','sec','sec','mjd','','m/s^2'};
%           labels={str.clean(filename,{'file','grace','.'}),'Job ID','arc duration','arc start','arc t0','arc nr','TBD'};
%           time_fh=@(raw) datetime(raw{7}+raw{9}/seconds(days(1)),...
%             'ConvertFrom','modifiedjuliandate'...
%           );
%           data_fh=@(raw) [raw{14},double(raw{11}),raw{10},raw{9},time.mjd(time.utc2gps(time.ToDateTime(double(raw{16}),'J2000sec'))),double(raw{8}),double(raw{15})];
%           timesystem='gps';
%           sanity_check=@(raw) all(all([ raw{1}-2000==raw{6},raw{2}==raw{4},raw{3}==raw{5}]));
%         end
%
%         if isempty(fmt)
%           error([mfilename,': cannot handle the GraceAccCal file ''',filename,'''.'])
%         end
%         %reading data
%         fid = file.open(filename);
%         raw = textscan(fid,fmt,'delimiter',' ','MultipleDelimsAsOne',1);
%         fclose(fid);
%         %keep some sanity
%         assert(sanity_check(raw),['Failed sanity check on ',filename]);
%         %building time domain
%         t=time_fh(raw);
%         %building data domain
%         y=data_fh(raw);
%         %sanity
%         if isempty(t) || isempty(y)
%           disp([mfilename,': this file has no data  ',filename])
%           obj=[];
%         else
%           iter=0;
%           while any(diff(t)==0)
%             %loop inits
%             n0=numel(t);
%             iter=iter+1;
%             %need to remove duplicate entries with different job IDs
%             mask=true(size(t));
%             for i=2:numel(t)
%               %get rid of those entries with zero or negative time stamp delta and lower ID
%               if t(i)<=t(i-1) && mask(i)
%                 if y(i,2) > y(i-1,2)
%                   mask(i-1)=false;
%                 else
%                   mask(i)=false;
%                 end
%               end
%             end
%             t=t(mask);
%             y=y(mask,:);
%             disp(['At iter ',num2str(iter),', removed ',num2str(n0-numel(t),'%04d'),' duplicate time entries (',filename,').'])
%           end
%           %need to monotonize the data (sometimes the entries are ordered according to arc number and not chronologically)
%           if any(diff(t)<0)
%             [t,i]=sort(t);
%             y=y(i,:);
%             disp(['Sorted ',num2str(sum(i~=transpose(1:numel(i))),'%04d'),' time entries (',filename,').'])
%           end
%           %building object
%           obj=simpletimeseries(t,y,...
%             'format','datetime',...
%             'labels',labels,...
%             'units',units,...
%             'timesystem',timesystem,...
%             'descriptor',filename,...
%             'monotonize','remove'...
%           );
%         end