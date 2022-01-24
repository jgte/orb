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
      test_list={'CSR2x2','CSR5x5'};
      switch(method)
        case 'all'
          for i=1:numel(test_list)
            out{i}=slr(test_list{i},'start',start,'stop',stop); %#ok<AGROW>
          end
        case {'CSR5x5','CSR2x2','GSFC5x5'}
          out=slr(method);
          out.plot('method','timeseries','degrees',[2,2,2,2,2],'orders',[-2,-1,-0,1,2],'zeromean',true);
        otherwise
          error(['Cannot handle test method ''',method,'''.'])
      end
    end
  end
  methods
    %% constructor
    function obj=slr(source,varargin)
      %NOTICE: to define a start/stop, pass it in varargin; gravity.common_ops will handle that
      % input parsing
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('source',@ischar)
      %declare parameters p
      v=varargs.wrap('parser',p,'sources',{slr.parameters('obj')},'mandatory',{source},varargin{:});
      %branch on type of SLR data
      switch source
        case 'CSR2x2'
          [t,y]=import_CSR2x2(v.varargin{:});
        case 'CSR5x5'
          [t,y]=import_CSR5x5(v.varargin{:});
        case 'GSFC5x5'
          [t,y]=import_GSFC5x5(v.varargin{:});
        otherwise
          error(['Cannot handle SLR data of type ''',source,'''.'])
      end
      obj=obj@gravity(t,y,varargin{:},'descriptor',['SLR ',source]);
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

%TODO: need to retreive y_out_error in this function
function [t_out,y_out]=import_CSR2x2(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   fullfile(slr.dir('data'),'csr','2x2'),@ischar;...
      'format',        'csr',@ischar;...
      'data_dir_url', 'http://ftp.csr.utexas.edu/pub/slr/degree_2',@ischar;...
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
    %download the data (done inside file.unwrap)
    local_data=file.unwrap(local_data,'remote_url',v.data_dir_url,'scalar_as_strings',true,varargin{:});
    %check the format of the data
    if file.isext(local_data,'.mat')
       %load the data
      [out,loaded_flag]=file.load_mat(local_data,'data_var','out');
      assert(loaded_flag,['BUG TRAP: could not load the data from ',local_data])
      %unpack the data
      t_out=out.t_out;
      y_out=out.y_out;
    else
      %need to make sure file.unwrap returned the txt file
      assert(file.isext(local_data,'.txt'),'BUG TRAP: expecting a txt file')
      %load the header (search for the end of the header up until line 30)
      header=file.header(local_data,30);
      %branch on files with one or two coefficients
      if contains(header,'C21') || contains(header,'C22')
        %2002.0411  2.43934614E-06 -1.40026049E-06  0.4565  0.4247 -0.0056  0.1782   20020101.0000   20020201.0000
        file_fmt='%f %f %f %f %f %f %f %f %f';
        data_cols=[2 3];
        sigm_cols=[4 5];
        corr_cols=[6 7];
%         units={'',''};
%         if ~contains(header,'C21')
%           labels={'C2,1','C2,-1'};
%         else
%           labels={'C2,2','C2,-2'};
%         end
      else
        %2002.0411  -4.8416939379E-04  0.7852  0.3148  0.6149   20020101.0000   20020201.0000
        file_fmt='%f %f %f %f %f %f %f';
        data_cols=2;
        sigm_cols=4;
        corr_cols=5;
%         units={''};
%         if contains(header,'C20'); labels={'C2,0'}; end
%         if contains(header,'C40'); labels={'C4,0'}; end
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
      %pack the data
      out=struct('t_out',t_out,'y_out',y_out);
      %save the data in mat format
      file.save_mat(out,local_data,'data_var','out')
    end
  end  
end
function [t_out,y_out,y_out_error,y_out_AOD,header]=import_CSR5x5(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   fullfile(slr.dir('data'),'csr','5x5'),@ischar;...
      'data_dir_url', 'http://ftp.csr.utexas.edu/pub/slr/degree_5',@ischar;...
      'data_file',    'CSR_Monthly_5x5_Gravity_Harmonics.txt',@ischar;...
      'data_labels',  {'n','m','Cnm','Snm','Cnm+AOD','Snm+AOD','C-sigma','S-sigma','Year_mid_point'},@iscellstr;
      'lmax',         5,   @(i) isscalar(i) && isnumeric(i);...
    },...
  },varargin{:});
  %define the local data file name
  local_data=fullfile(v.import_dir,v.data_file);
  %download the data (done inside file.unwrap)
  local_data=file.unwrap(local_data,'remote_url',v.data_dir_url,'scalar_as_strings',true,varargin{:});
  %check the format of the data
  if file.isext(local_data,'.mat')
     %load the data
    [out,loaded_flag]=file.load_mat(local_data,'data_var','out');
    assert(loaded_flag,['BUG TRAP: could not load the data from ',local_data])
    %unpack the data
    t_out      =out.t_out;
    y_out      =out.y_out;
    y_out_AOD  =out.y_out_AOD;
    y_out_error=out.y_out_error;
  else
    %need to make sure file.unwrap returned the txt file
    assert(file.isext(local_data,'.txt'),'BUG TRAP: expecting a txt file')
    %declare header structure
    header=struct('GM',0,'radius',0,'lmax',v.lmax,'tide_system','unknown','modelname','unknown',...
      'static',[],'labels',{{}},'idx',struct([]),'units',1);
    %define known details
    header.modelname='UT/CSR monthly 5x5 gravity harmonics';
    %open the file
    fid=file.open(local_data);
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
      if (keyword_search(s, 'Units'))
        header.units=str2double(strtrim(strrep(s,'Units','')));
      end
      if (keyword_search(s, 'Coefficients:'))
        header.labels=cells.rm_empty(strsplit(strtrim(str.rep(s,'Coefficients:','',',',''))));
        %build index records
        for i=1:numel(v.data_labels)
          %build fieldname (label may have ilegal characters)
          fn=str.clean(v.data_labels{i},'fieldname');
          %find where this data label is located in the header labels
          header.idx(1).(fn)=cells.strequal(header.labels,v.data_labels{i});
          %make sure it is found
          assert(~isempty(header.idx(1).(fn)),['Cannot find reference to data label ''',...
            v.data_labels{i},''' in the header of ',local_data])
        end
      end
      if (keyword_search(s, '===================='))
        %init loop variables
        counter=0;
        %read the static coefficients
        while true
          %get the next line
          s=fgets(fid);
          %split line into columns and remove empty entries
          s=cells.rm_empty(strsplit(s));
          %exist criteria
          if numel(s)~=6
            break
          end
          %increment line counter
          counter=counter+1;
          %save the data
          header.static(counter,:)=cells.c2m(cells.num(s));
        end
        %convert the header.static data to gravity-friendly format
        for j=1:size(header.static,1)
          d=header.static(j,1);
          %skip if this degree is above the requested lmax
          if d>header.lmax; continue; end
          %cosine coefficients are in the 3rd column (errors in the 5th)
          o=header.static(j,2);
          static_signal(gravity.colidx(d,o,header.lmax))=header.static(j,3); %#ok<AGROW>
          static_error( gravity.colidx(d,o,header.lmax))=header.static(j,5); %#ok<AGROW>
          if o==0, continue;end
          %sine coefficients are in the 4th column (errors in the 6th)
          o=-o;
          static_signal(gravity.colidx(d,o,header.lmax))=header.static(j,4); %#ok<AGROW>
          static_error( gravity.colidx(d,o,header.lmax))=header.static(j,6); %#ok<AGROW>
        end
      end
    end
    % some sanity
    assert(~isempty(header.static),'Could not retrieve the static gravity field coefficients')
    % read data
    while true
      s=fgets(fid);
      if ~ischar(s)
        break
      end
      %split line into columns and remove empty entries
      s=cells.c2m(cells.num(cells.rm_empty(strsplit(s))));
      %branch on the number of columns
      switch numel(s)
      case 7
        arc=s(1);
      case 10
        %get degree and order
        d=s(header.idx.n);
        o=s(header.idx.m);
        %skip if this degree is above the requested lmax
        if d>header.lmax; continue; end
        %get the epoch 
        t_out(arc)=datetime([0 0 0 0 0 0])+years(s(header.idx.Year_mid_point)); %#ok<AGROW>
        %save cosine coefficient, error and value with AOD
        y_out(      arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Cnm);    %#ok<AGROW>
        y_out_error(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Csigma); %#ok<AGROW>
        y_out_AOD(  arc,gravity.colidx(d,o,header.lmax))=s(header.idx.CnmAOD); %#ok<AGROW>
        if o==0, continue;end
        %save sine coefficient, error and value with AOD
        o=-o;
        y_out(      arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Snm);    %#ok<AGROW>
        y_out_error(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Ssigma); %#ok<AGROW>
        y_out_AOD(  arc,gravity.colidx(d,o,header.lmax))=s(header.idx.SnmAOD); %#ok<AGROW>
      otherwise
        disp(['WARNING: ignoring line: ',strjoin(s,' ')])
      end
    end
    fclose(fid);
    %add static signal and error
    units=header.units*ones(size(y_out,1),1);
    y_out      =     y_out          +  units*static_signal;
    y_out_AOD  =     y_out_AOD      +  units*static_signal;
    y_out_error=sqrt(y_out_error.^2 + (units*static_error).^2);
    %pack the data
    out=struct('t_out',t_out,'y_out',y_out,'y_out_AOD',y_out_AOD,'y_out_error',y_out_error);
    %save the data in mat format
    file.save_mat(out,local_data,'data_var','out')
  end
end
function [t_out,y_out,y_out_error,y_out_AOD,header]=import_GSFC5x5(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   fullfile(slr.dir('data'),'gsfc','5x5'),@ischar;...
      'data_dir_url', 'https://earth.gsfc.nasa.gov/sites/default/files/2021-12',@ischar;...
      'data_file',    'GSFC_SLR_5x5c61s61_200001_202111.txt',@ischar;...
      'data_labels',  {'n','m','C','S'},@iscellstr;
      'lmax',         5,   @(i) isscalar(i) && isnumeric(i);...
    },...
  },varargin{:});
  %define the local data file name
  local_data=fullfile(v.import_dir,v.data_file);
  %download the data (done inside file.unwrap)
  local_data=file.unwrap(local_data,'remote_url',v.data_dir_url,'scalar_as_strings',true,varargin{:});
  %check the format of the data
  if file.isext(local_data,'.mat')
     %load the data
    [out,loaded_flag]=file.load_mat(local_data,'data_var','out');
    assert(loaded_flag,['BUG TRAP: could not load the data from ',local_data])
    %unpack the data
    t_out=out.t_out;
    y_out=out.y_out;
  else
    %need to make sure file.unwrap returned the txt file
    assert(file.isext(local_data,'.txt'),'BUG TRAP: expecting a txt file')
    %declare header structure
    header=struct('GM',0,'radius',0,'lmax',v.lmax,'tide_system','unknown','modelname','unknown',...
      'labels',{{}},'idx',struct([]),'units',1);
    %define known details
    header.modelname='GSFC weekly 5x5 gravity harmonics';
    %open the file
    fid=file.open(local_data);
    % Read header
    while true
      s=fgets(fid); 
      if keyword_search(s,'end of header')
        break
      end
      if keyword_search(s,'GM:')
        l=strsplit(s);
        header.GM = str2double(l{3});
      end
      if (keyword_search(s, 'R:'))
        l=strsplit(s);
        header.radius=str2double(l{3});
      end
      if (keyword_search(s, 'C20 is'))
        header.tide_system=strtrim(strrep(s,'C20 is',''));
      end
      if (keyword_search(s, 'Coefficient lines:'))
        header.labels=cells.rm_empty(...
          strsplit(...
            str.clean(...
              s,{...
                '(I3,I3,1X,E20.13,1X,E20.13)','Coefficient','lines:',','...
              }...
            )...
          )...
        );
        %build index records
        for i=1:numel(v.data_labels)
          %build fieldname (label may have ilegal characters)
          fn=str.clean(v.data_labels{i},'fieldname');
          %find where this data label is located in the header labels
          %NOTICE: -2 is needed to offset the two words 'Coefficient lines:'
          header.idx(1).(fn)=cells.strequal(header.labels,v.data_labels{i});
          %make sure it is found
          assert(~isempty(header.idx(1).(fn)),['Cannot find reference to data label ''',...
            v.data_labels{i},''' in the header of ',local_data])
        end
      end
      if keyword_search(s, 'Product:')
        break      
      end
    end
    %init loop variables
    arc=0;
    % read data
    while true
      s=fgets(fid);
      if ~ischar(s)
        break
      end
      %split line into columns and remove empty entries
      s=cells.c2m(cells.num(cells.rm_empty(strsplit(s))));
      %branch on the number of columns
      switch numel(s)
      case 2
        %increment loop var
        arc=arc+1;
        %get the epoch 
        t_out(arc)=datetime([0 0 0 0 0 0])+years(s(2)); %#ok<AGROW>
      case 4
        %get degree and order
        d=s(header.idx.n);
        o=s(header.idx.m);
        %skip if this degree is above the requested lmax
        if d>header.lmax; continue; end
        %save cosine coefficient
        y_out(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.C);    %#ok<AGROW>
        if o==0, continue;end
        %save sine coefficient
        o=-o;
        y_out(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.S);    %#ok<AGROW>
      otherwise
        disp(['WARNING: ignoring line: ',strjoin(s,' ')])
      end
    end
    fclose(fid);
    %pack the data
    out=struct('t_out',t_out,'y_out',y_out);
    %save the data in mat format
    file.save_mat(out,local_data,'data_var','out')
  end
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