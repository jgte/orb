classdef grace
  properties(Constant)
    l1bdir_options={...
      '~/data/grace';...
      './data/grace';...
    };
  end
  methods(Static)
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
    %% exporting data
    function export(obj,filename,filetype,varargin)
      %NOTICE: this function is a functional duplicate of simpletimeseries.export, but
      %        implements GRACE-specific formats
      %TODO: needs testing
      p=inputParser;
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