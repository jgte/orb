classdef slr < gravity
  properties(Constant)
    data_options={...
      '~/data/SLR';...
      './data/SLR';...
    };
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
    type
  end
  %calculated only when asked for
  properties(Dependent)
    time
  end
  methods(Static)
    %% directories
    function out=dir(type)
      switch type
        case 'data'
          for i=1:numel(slr.data_options)
            if file.exist(slr.data_options{i})
              out=slr.data_options{i};
              return
            end
          end
        %add more directories here
      end
    end
    %% interface methods to object constants
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(slr.parameter_list); end
      out=v.picker(varargin{:});
    end
    %% testing for the current object
    function out=test_parameters(field)
      switch lower(field)
      case 'start'
        out=datetime('2010-01-01');
      case 'stop'
        out=datetime('2020-12-31');
      otherwise
        error([mfilename,': unknown field ',field,'.'])
      end
    end
    function out=test(method)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      start=slr.test_parameters('start');
      stop= slr.test_parameters('stop' );
      test_list={'CSR2x2'};
      switch(method)
        case 'all'
          for i=1:numel(test_list)
            out{i}=slr(test_list{i},'start',start,'stop',stop); %#ok<AGROW>
          end
        case {'CSR5x5','CSR2x2'}
          out=slr(method);
          out.plot('method','timeseries','degrees',2,'orders',0);
        otherwise
          error(['Cannot handle test method ''',method,'''.'])
      end
    end
  end
  methods
    %% constructor
    function obj=slr(type,varargin)
      %NOTICE: to define a start/stop, pass it in varargin; gravity.common_ops will handle that
      %declare parameters p
      v=varargs.wrap('sources',{slr.parameters('obj')},'mandatory',{type},varargin{:});
      %branch on type of SLR data
      switch type
        case 'CSR2x2'
          [t,y]=import_CSR2x2(v.varargin{:});
        case 'CSR5x5'
          [t,y]=import_CSR5x5(v.varargin{:});
        otherwise
          error(['Cannot handle SLR data of type ''',type,'''.'])
      end
      obj=obj@gravity(t,y,varargin{:},'descriptor',['SLR ',type]);
      %apply model processing options
      obj=gravity.common_ops('all',obj,v.varargin{:},'product_name',obj.descriptor);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      warning off MATLAB:structOnObject
      out=varargs(...
        structs.filter(struct(obj),[slr.parameters('list');more_parameters(:)])...
      ).varargin;
      warning on MATLAB:structOnObject
    end
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=20;
      end
      disp(' --- Parameters --- ')
      for i=1:numel(slr.parameters('list'))
        %shorter names
        p=slr.parameters('value',i);
        disp([p,repmat(' ',1,tab-length(p)),' : ',str.show(obj.(p))])
      end
      d_list=slr.data_types;
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

function [t_out,y_out]=import_CSR2x2(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   '/tmp',@ischar;...
      'format',        'csr',@ischar;...
      'suffix',       'RL06',@ischar;...
      'prefixes',{'C20','C21_S21','C22_S22'},@iscellstr;
      'degrees', [    2,        2,       2] ,@isnumeric;
      'orders',  {    0,[  1, -1],[  2, -2]},@isnumeric;
    },...
  },varargin{:});
  %sanity
  assert(numel(v.prefixes)==numel(v.degrees)&& numel(v.prefixes)==numel(v.orders),...
    'Metadata ''prefixes'' ''degrees'' and ''orders'' must all have the same size.') 
  %init records
  t_out=cell(size(v.prefixes));
  y_out=[];
  %loop over all data
  for i=1:numel(v.prefixes)
    %define the data file name
    data_file=[v.prefixes{i},'_',v.suffix,'.txt'];
    local_data=fullfile(v.import_dir,data_file);
    if file.exist(local_data)
      disp(['NOTICE: not downloading updated SLR data because file exists: ',local_data])
    else
      %download the data
      file.system(...
        ['wget http://ftp.csr.utexas.edu/pub/slr/degree_2/',data_file],...
        'cd',v.import_dir,...
        'disp',true,...
        'stop_if_error',true...
      );
    end
    %load the header
    header=file.header(local_data,20);
    %branch on files with one or two coefficients
    if contains(header,'C21') || contains(header,'C22')
      %2002.0411  2.43934614E-06 -1.40026049E-06  0.4565  0.4247 -0.0056  0.1782   20020101.0000   20020201.0000
      file_fmt='%f %f %f %f %f %f %f %f %f';
      data_cols=[2 3];
      sigm_cols=[4 5];
      corr_cols=[6 7];
%       units={'',''};
%       if ~contains(header,'C21')
%         labels={'C2,1','C2,-1'};
%       else
%         labels={'C2,2','C2,-2'};
%       end
    else
      %2002.0411  -4.8416939379E-04  0.7852  0.3148  0.6149   20020101.0000   20020201.0000
      file_fmt='%f %f %f %f %f %f %f';
      data_cols=2;
      sigm_cols=4;
      corr_cols=5;
%       units={''};
%       if contains(header,'C20'); labels={'C2,0'}; end
%       if contains(header,'C40'); labels={'C4,0'}; end
    end
    raw=file.textscan(local_data,file_fmt);
    %build the time domain
    t=datetime([0 0 0 0 0 0])+years(raw(:,1));
    %building data domain
    switch v.format
    case 'csr'       %NOTICE: this includes AOD mean for the solution period (aka correction, or 'corr')
      y=raw(:,data_cols);
    case 'csr-grace' %NOTICE: this removes the AOD correction from SLR, making it more suitable to replace the GRACE C20
      y=raw(:,data_cols)-raw(:,corr_cols)*1e-10;
    case 'csr-corr'  %NOTICE: this is the AOD correction
      y=raw(:,corr_cols)*1e-10;
    case 'csr-sigma' %NOTICE: this is the solution sigma
      y=raw(:,sigm_cols)*1e-10;
    end
    %sanity
    if i==1
      t_out=t;
    else
      assert(~any(~simpletimeseries.ist('==',t,t_out)),'time domain inconsistency')
    end
    %building aggregated records
    for j=1:numel(v.orders{i})
      d=v.degrees(i);
      o=v.orders{i}(j);
      y_out(:,gravity.colidx(d,o,max(v.degrees)))=y(:,j); %#ok<AGROW>
    end
  end  
end


function [t_out,y_out,header]=import_CSR5x5(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   '/tmp',@ischar;...
    },...
  },varargin{:});
  %define the data file name
  data_file='CSR_Monthly_5x5_Gravity_Harmonics.txt';
  local_data=fullfile(v.import_dir,data_file);
  if file.exist(local_data)
    disp(['NOTICE: not downloading updated SLR data because file exists: ',local_data])
  else
    %download the data
    file.system(...
      ['wget http://ftp.csr.utexas.edu/pub/slr/degree_2/',data_file],...
      'cd',v.import_dir,...
      'disp',true,...
      'stop_if_error',true...
    );
  end
  %define known details
  modelname='UT/CSR monthly 5x5 gravity harmonics';
  %open the file
  fid=file.open(local_data);
  header=struct('GM','radius','Lmax','tide_system='};
  % Read header
  while true
     s=fgets(fid); 
     if keyword_search(s,'end of header')
       break
     end
     if keyword_search(s,'earth_gravity_constant')
        header.GM = str2double(strtrim(strrep(s,'earth_gravity_constant','')));
     end
     if (keyword_search(s, 'radius'))
        header.radius=str2double(strtrim(strrep(s,'radius','')));
     end
     if (keyword_search(s, 'tide_system'))
        header.tide_system=strtrim(strrep(s,'tide_system',''));
     end
  end
  % read data
  while true
    s=fgets(fid);
    if ~ischar(s)
      break
    end
    %split line into columns
    s=strsplit(s);
    %branch on the number of columns
    switch numel(s)
      case 7
      %do nothing
    case 10
      continuar aqui
    otherwise
      disp(['WARNING: ignoring line: ',strjoin(s,' ')])
    end
  end
  fclose(fid)
end

%% Aux functions
function out=keyword_search(line,keyword)
    out=strncmp(strtrim(line),       keyword,         length(keyword));
end

%NOTICE: this was used to load data in the following format:
% MJD      Year       c30       c31      c32      c33      s31      s32      s33   sc20   sc31  sc32  sc33  ss31  ss32  ss33
% 52275.0 2002.0000   1.7548   0.1293   0.7154   0.1866   0.6015  -0.3642   0.7701  1.61  2.38  4.04  6.09  2.27  3.95  5.75
% 52306.0 2002.0849   1.6075  -0.5264   0.0338  -2.0945   0.0514   0.7133   0.8110  1.55  2.20  3.55  5.21  2.14  3.65  5.64
% 52334.0 2002.1615   2.1569  -0.2098   0.8088  -1.4374  -0.1334   0.8748   1.4793  1.37  2.15  3.66  5.51  2.08  3.61  5.46
% 52365.0 2002.2464   1.9519   0.1812   1.1169  -0.1042   0.1096   1.0244   1.2449  1.66  2.19  3.82  5.73  2.24  3.80  5.52
% 52395.0 2002.3285   0.6036  -0.4515   0.3466   0.4565  -0.7234   0.9718   1.7062  1.60  2.21  3.93  5.74  2.28  3.93  5.62
% 52426.0 2002.4134   0.9127  -0.5447   1.3512  -0.2333  -1.2079  -0.3675   0.8941  1.53  2.09  3.72  5.48  2.13  3.65  5.23
% 52456.0 2002.4956   0.0842  -0.8563   1.0261  -0.7019  -0.9459  -0.5351   0.3056  1.57  2.26  4.11  5.83  2.35  4.26  5.73
% 52487.0 2002.5804  -0.3460  -0.3462   0.1635  -0.6590  -0.7002  -0.4621  -0.9735  1.49  2.18  3.68  5.22  2.12  3.78  5.36
% 52518.0 2002.6653  -1.1352   0.1460  -0.1374   0.0289  -0.7295  -0.1479  -0.7320  1.65  2.22  3.83  5.27  2.20  3.77  5.49
%NOTICE: these data were located at /corral-tacc/utexas/csr/byaa664/run17/grace/lisjob/geosol/csdeg3.5d561_192_stb
function obj=import_slr_Cheng(obj,product,varargin)
  obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
  % sanity
  assert(time.isfinite(obj.start) && time.isfinite(obj.stop),'Need valid obj.start and obj.stop.')
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',        '',@ischar;...
      'filename',          '',@ischar;...
      'format','slr-csr-Chen',@ischar;...
    },...
    product.args...
  },varargin{:});
  %sanity
  assert(numel(v.degrees)==numel(v.orders),...
    'Metadata ''degrees'' and ''orders'' must all have the same size.')
  %define data cols
  data_cols=3:9;
  %load the data
  [raw,header]=file.textscan(fullfile(v.import_dir,v.filename),...
    '%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
  %get the labels (from the header)
  labels=strsplit(header,' ');
  %build units
  units=cell(size(data_cols));
  units(:)={''};
  t=datetime([0 0 0 0 0 0])+years(raw(:,2));
  %building object
  obj=simpletimeseries(t,raw(:,data_cols)*1e-10,...
    'format','datetime',...
    'units',units,...
    'labels', labels(data_cols),...
    'timesystem','gps',...
    'descriptor',['SLR Stokes coeff. from ',filename]...
  );
  %init gravity object
  g=gravity.unit(max(v.degrees),'scale',0,'t',ts.t);
  %propagate the data
  g=g.setC(v.degrees,v.orders,ts.y,ts.t);
  %apply model processing options
  g=gravity.common_ops('all',v,product,g);
  %save model
  obj=obj.data_set(product.dataname,g);
  obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
end