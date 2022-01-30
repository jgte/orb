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
      switch type %NEEDS UPDATE WHEN ADDING A NEW MODEL
        case 'data'
          for i=1:numel(slr.data_options)
            if file.exist(slr.data_options{i})
              out=slr.data_options{i};
              return
            end
          end
        case 'CSR2x2';         out=fullfile(slr.dir('data'),'csr' ,'2x2'  );
        case 'CSR5x5';         out=fullfile(slr.dir('data'),'csr' ,'5x5'  );
        case 'TN-07';          out=fullfile(slr.dir('data'),'csr' ,'TN-07');
        case 'TN-11';          out=fullfile(slr.dir('data'),'csr' ,'TN-11');
        case 'CSR-RL06';       out=fullfile(slr.dir('data'),'csr' ,'RL06' );
        case 'GSFC5x5';        out=fullfile(slr.dir('data'),'gsfc','5x5'  );
        case 'TN-14';          out=fullfile(slr.dir('data'),'gsfc','TN-14');
        case 'GSFC';           out=fullfile(slr.dir('data'),'gsfc','C20'  );
        case 'GSFC-7DAY';      out=fullfile(slr.dir('data'),'gsfc','7day' );
        %add more directories here
        otherwise
          error(['Cannot handle type ''',type,'''.'])
      end
    end
    %% interface methods to object constants
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(slr.parameter_list); end
      out=v.picker(varargin{:});
    end
    %% retrieves the Monthly estimates of C20 from 5 SLR satellites based on GRACE RL05/RL06 models
    function out=graceC20(varargin)
      %parse arguments that are required later
      v=varargs.wrap('sources',{...
        {...
          'time',        [], @(i) isdatetime(i) || isempty(i);...
          'mode',    'read', @ischar;...
          'descriptor',  '', @ischar;...
          'source', 'TN-14', @(i) ischar(i) || iscellstr(i);...
          'start',time.zero_date, @isdatetime;...
          'stop',  time.inf_date, @isdatetime;...
          'C20mean',-4.8416945732000E-04, @isscalar;...
        },...
      },varargin{:});
      %parse using a model or the original data
      if contains(v.source,'-model')
        new_mode=['model-',strrep(v.mode,'model-','')];
        new_version=strrep(v.source,'-model','');
        str.say('WARNING: over-writing input mode',str.quote(v.mode),'with',str.quote(new_mode),...
          'since input version is',str.quote(v.source),', now passed along as',str.quote(new_version),'.')
        out=gravity.graceC20(varargin{:},'mode',new_mode,'version',new_version);
        return
      end
      switch v.mode
      case 'model-poly'
        out=2;
      case 'model-periods-datfile'
        [p,n]=fileparts(GetGRACEC20('mode','data_file','version',v.source));
        out=fullfile(p,[n,'_periods.mat']);
      case {'model-compute','model-periods'}
        %get data file
        f=gravity.graceC20(varargin{:},'mode','model-periods-datfile');
        %check if periods were already computed
        if ~file.exist(f)
          %loading necessary data
          c20=gravity.graceC20(varargin{:},'mode','read');
           np=gravity.graceC20(varargin{:},'mode','model-poly');
          %compute periods
          [~,pd]=c20.parametric_decomposition_search('np',np,'T',[365.2426,182.6213],'timescale','days');
          out=pd.T;
          %save periods
          save(f,'out');
        else
          disp(['Loading parametric decomposition from ',f])
          %load periods
          load(f,'out');
        end
%         %build strings to save to file
%         h=strsplit(c20.descriptor,newline);
%         s=strjoin([...
%             cellfun(@(i) ['%',i],h(:),'UniformOutput',false);...
%             {'out=[...'};...
%             arrayfun(@(i) ['  ',num2str(i,'%.12e'),';...'],out(:),'UniformOutput',false);...
%             {'];'};...
%           ],newline);
%         f=gravity.graceC20(varargin{:},'mode','model-periods-datfile');
%         b=file.strsave(f,s);
%         str.say('Written',b,'bytes of data to file:',newline,f,newline,'related to the periods of:',newline,c20.descriptor);
      case 'model-datfile'
        [p,n]=fileparts(GetGRACEC20(varargin{:},'mode','data_file'));
        out=fullfile(p,[n,'_pd.mat']);
      case 'model-md5file'
        [p,n]=fileparts(GetGRACEC20(varargin{:},'mode','data_file'));
        out=fullfile(p,[n,'.md5']);
      case 'model-md5'
        out=file.md5(GetGRACEC20(varargin{:},'mode','data_file'));
      case 'model-md5set'
        out=file.strsave(...
          gravity.graceC20(varargin{:},'mode','model-md5file'),...
          gravity.graceC20(varargin{:},'mode','model-md5')...
        );
      case 'model-md5get'
        md5file=gravity.graceC20(varargin{:},'mode','model-md5file');
        if ~file.exist(md5file)
          gravity.graceC20(varargin{:},'mode','model-md5set');
        end
        out=file.strload(md5file);
      case 'model-md5check'
        out=strcmp(...
          gravity.graceC20(varargin{:},'mode','model-md5get'),...
          gravity.graceC20(varargin{:},'mode','model-md5')...
        );
      case {'model','model-get','model-set','model-read','model-reload'}
        %loading necessary data
        c20=gravity.graceC20(varargin{:},'mode','read');
         np=gravity.graceC20(varargin{:},'mode','model-poly');
          T=gravity.graceC20(varargin{:},'mode','model-periods');
        %check if pdset is already available
        f_pdset=gravity.graceC20(varargin{:},'mode','model-list-datfile');
        f_pd   =gravity.graceC20(varargin{:},'mode','model-datfile');
        if ~file.exist(f_pdset) || ~file.exist(f_pd) || ~...
          gravity.graceC20(varargin{:},'mode','model-md5check')  
          %get the coefficients; NOTICE: always use c20.t so that f_pdset is not dependent on inputs
          [~,pd_set]=c20.parametric_decomposition('np',np,'T',T,...
            'timescale','days','time',c20.t_domain(days(7))); 
          %save them
          save(f_pdset,'pd_set')
          %update md5 of data
          gravity.graceC20(varargin{:},'mode','model-md5set')  
        else
          load(f_pdset,'pd_set')
        end
        %handle different time domains
        if isempty(v.time);      v.time=c20.t; end
        %implement default start/stop for modelling
        if time.iszero(v.start); v.start=v.time(1  ); end
        if time.isinf( v.stop );  v.stop=v.time(end); end
        %enforce start/stop
        v.time=v.time(v.time>=v.start & v.time<=v.stop);
        %evaluate model at requested time domain
        out=pardecomp.join( pd_set,'time',v.time);
      case 'model-list-datfile'
        [p,n]=fileparts(GetGRACEC20(varargin{:},'mode','data_file'));
        out=fullfile(p,[n,'_pdset.mat']);
      case {'model-list','model-list-tex'}
        %check if pdset is already available
        f=gravity.graceC20(varargin{:},'mode','model-list-datfile');
        if ~file.exist(f)
          %compute the model (saving is done inside)
          gravity.graceC20(varargin{:},'mode','model');
        end
        %load it
        load(f,'pd_set')
        %output it
        switch v.mode
          case 'model-list';     out=pardecomp.table(pd_set,'tablify',true);
          case 'model-list-tex'; out=pardecomp.table(pd_set,'tablify',false,'latex_table',true);
        end       
      case 'model-plot'
        %retrieve the orignal data
        c20o=gravity.graceC20(varargin{:},'mode','read');
        %resample to a finer time domain
        c20r=c20o.interp(c20o.t_domain(days(7)),...
          'interp_over_gaps_narrower_than',days(45),...
          'interp1_args',{'pchip'}...
        );
        c20m=gravity.graceC20(varargin{:},'mode','model','time',c20o.t_domain(days(7)));
        c20e=c20r-c20m;
        %compute means
        c20om=c20o.stats('mode','mean');
        c20mm=c20m.stats('mode','mean');
        c20em=c20e.stats('mode','mean');
        %compute stds
        c20os=c20o.stats('mode','std');
        c20ms=c20m.stats('mode','std');
        c20es=c20e.stats('mode','std');
        %remove the mean
        c20o=c20o-c20om(1);
        c20m=c20m-c20mm(1);
        c20e=c20e-c20em(1);
        %plot it
        plotting.figure;
        c20o.addgaps(days(35)).plot('columns',1);
        c20m.plot('columns',1);
        c20e.plot('columns',1);
        plotting.enforce(varargin{:},...
          'plot_legend',...
          {...
            ['Original \mu=',num2str(c20om(1),'%.3e'),' \sigma=',num2str(c20os(1),'%.3e')],...
            ['Modelled \mu=',num2str(c20mm(1),'%.3e'),' \sigma=',num2str(c20ms(1),'%.3e')],...
            ['Residual \mu=',num2str(c20em(1),'%.3e'),' \sigma=',num2str(c20es(1),'%.3e')],...
          },...
          'plot_xlabel','none',...
          'plot_title',c20o.descriptor...
        );
        out=c20e;
      case 'interp'
        assert(~isempty(v.time),['If mode is ''',v.mode,''', need argument ''time''.'])
        out=gravity.graceC20(varargin{:},'mode','read');
        out=out.interp(v.time);
      case 'plot-all'
        if iscellstr(v.source)
          version_list=v.source;
        else
          version_list={'GSFC-7day','GSFC','CSR-RL06','TN-14','TN-11','TN-07'};
        end
        dat_list=cell(size(version_list));
        plotting.figure(varargin{:});
        for i=1:numel(version_list)
          %pick keywords
          kw=strsplit(version_list{i},'-');
          %branch on them (NOTICE: always last)
          switch kw{end}
          case 'model'
            if i>1
              %NOTICE: to get the model with the most complete time domain, out it last in the 'version' input
              dat_list{i}=gravity.graceC20('mode','model',...
                'version',strrep(version_list{i},'-model',''),...
                'time',time.union(cellfun(@(j) j.t,dat_list(1:i-1),'UniformOutput',false),days(7))...
              );
            else
              dat_list{i}=gravity.graceC20('mode','model',...
                'version',strrep(version_list{i},'-model','')...
              );
            end              
          otherwise
            dat_list{i}=gravity.graceC20('mode','read','version',version_list{i});
            dat_list{i}=dat_list{i}.interp(dat_list{i}.t_domain(days(7)),...
              'interp_over_gaps_narrower_than',days(45)...
            );
          end
          dat_list{i}=dat_list{i}-v.C20mean;
          dat_list{i}.plot('columns',1);
        end
        plotting.enforce(varargin{:},...
          'plot_legend',version_list,...
          'plot_xlabel','none');
        title(['mean C_{20}: ',num2str(v.C20mean,'%e')]);
        out=dat_list;
      otherwise
        %call mother routine
        try
          [t,s,e,d]=GetGRACEC20(varargin{:},'mode',v.mode);
        catch ME
          switch ME.identifier
            case 'MATLAB:unassignedOutputs'
              t=GetGRACEC20(varargin{:},'mode',v.mode);
              s=[];
            otherwise
              rethrow(ME)
          end
        end
        if ~isempty(s)
          %create time series
          out=simpletimeseries(t,[s,e],...
            'labels',{'C20','error C20'},...
            'units',{'',''},...
            'timesystem','gps',...
            'descriptor',d...
          );
        else
          out=t;
        end
      end
    end
    %% load predefined SLR data
    function obj=load(source,varargin)
     %handle producing the parametric model
      if contains(source,'-C20-model')
        %save data source
        data_source=strrep(source,'-C20-model','');
        %build parametric model data filename
        [p,n]=fileparts(fullfile(slr.dir(data_source),data_source));
        datafile=fullfile(p,[n,'_C20_pd.mat']);
        %check if data file is missing or it's old
        if ~file.exist(datafile) || file.age(datafile)>days(30)
          %need to clean vararing of 'start' and 'stop' so that the model can be built from all data
          varargin_now=cells.vararginclean(varargin,{'start','stop'});
          %load the data
          data=slr.load(data_source,varargin_now{:});
          %get the modelled coefficients; NOTICE: always use data.t so that f_pdset is not dependent on inputs
          [~,pd_set]=data.ts_C(2,0).parametric_decomposition_search(...
            'np',2,...
            'T',[1,1/2]*days(years(1)),...
            'timescale','days'...
          ); 
          %append important information
          pd_set.GM=data.GM;
          pd_set.R=data.R;
          pd_set.descriptor=data.descriptor;
          %save the data
          save(datafile,'pd_set')
        else
          %load the data
          load(datafile,'pd_set')
        end
        %evaluate model at the data source time domain
        model=pardecomp.join( pd_set,'time',pd_set.time);
        %propagate relevant information
        t=model.t;
        y=zeros(numel(t),gravity.y_length(2));
        y(:,gravity.colidx(2,0,2))=model.y;
        %TODO: test that the metadata gets all the way here
        header=struct(...
          'GM'  ,pd_set.GM,...
          'R'   ,pd_set.R,...
          'descriptor',['model of ',pd_set.descriptor]...
        );
      else
        %branch on type of SLR data
        switch source %NEEDS UPDATE WHEN ADDING A NEW MODEL
          case 'CSR2x2'
            [t,y,header]=import_CSR2x2(varargin{:});
          case 'CSR5x5'
            [t,y,header]=import_CSR5x5(varargin{:});
          case 'GSFC5x5'
            [t,y,header]=import_GSFC5x5(varargin{:});
          case {'TN-07','TN-11','CSR-RL06','TN-14','GSFC','GSFC-7DAY'}
            [t,y,header]=import_C20('source',source,varargin{:});
          otherwise
            error(['Cannot handle SLR data of type ''',source,'''.'])
        end
      end
      obj=slr(t,y,varargs(header).varargin{:});
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
      test_list={'CSR2x2','CSR5x5','GSFC5x5','TN-07','TN-11','CSR-RL06','TN-14','GSFC','GSFC-7DAY'};
      switch(method)
        case 'all'
          plotting.figure;
          cellfun(@(i) slr.load(...
            i,'start',start,'stop',stop...
          ).plot(...
            'method','timeseries','degrees',2,'orders',0 ...
          ),test_list);
          plotting.enforce(...
            'plot_line_color','spiral',...
            'plot_legend',test_list...
          );
        case 'C20-model'
          for i=1:numel(test_list)
            plotting.figure;
            slr.load( test_list{i},          'start',start,'stop',stop).plot('method','timeseries','degrees',2,'orders',0);
            slr.load([test_list{i},'-C20-model'],'start',start,'stop',stop).plot('method','timeseries');
            plotting.enforce(...
              'plot_legend',{test_list{i},'model'}...
            );
          end
        case {'CSR5x5','CSR2x2','GSFC5x5'} %MAY NEED UPDATE WHEN ADDING A NEW MODEL
          out=slr.load(method);
          out.plot('method','timeseries','degrees',[2,2,2,2,2],'orders',[-2,-1,-0,1,2],'zeromean',true);
        case {'TN-07','TN-11','CSR-RL06','TN-14','GSFC','GSFC-7DAY'} %MAY NEED UPDATE WHEN ADDING A NEW MODEL
          out=slr.load(method);
          out.plot('method','timeseries','degrees',2,'orders',0);
        otherwise
          error(['Cannot handle test method ''',method,'''.'])
      end
    end
  end
  methods
    %% constructor
    function obj=slr(t,y,varargin)
      %NOTICE: to define a start/stop, pass it in varargin; gravity.common_ops will handle that
      %input parsing
      p=machinery.inputParser;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y', @(i) simpledata.valid_y(i));
      % create argument object, declare and parse parameters, save them to obj
      [v,~]=varargs.wrap('parser',p,'sources',{slr.parameters('obj')},'mandatory',{t,y},varargin{:});
      %init the object
      %NOTICE: generally, the following options should be in varargin: 'GM', 'R' and 'descriptor'
      obj=obj@gravity(t,y,varargin{:});
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

%TODO: need to retreive y_out_error in all import_* function
function [t_out,y_out,header]=import_CSR2x2(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   slr.dir('CSR2x2'),@ischar;...
      'format',        'csr',@ischar;...
      'data_dir_url', 'http://ftp.csr.utexas.edu/pub/slr/degree_2',@ischar;...
      'suffix',       'RL06',@ischar;...
      'prefixes',{'C20','C21_S21','C22_S22'},@iscellstr;
      'degrees', [    2,        2,       2] ,@isnumeric; %NOTICE: this cannot be a cell because each entry must be scalar
      'orders',  {    0,[  1, -1],[  2, -2]},@(i) all(cellfun(@isnumeric,i));
    },...
  },varargin{:});
  %sanity
  assert(numel(v.prefixes)==numel(v.degrees)&& numel(v.prefixes)==numel(v.orders),...
    'Metadata ''prefixes'' ''degrees'' and ''orders'' must all have the same size.')
  %shortcuts
  lmax=max(v.degrees);
  %look for aggregated data file
  agg_data=fullfile(v.import_dir,'CSR2x2.mat');
  %check the format of the data
  if file.exist(agg_data)
     %load the data
    [out,loaded_flag]=file.load_mat(agg_data,'data_var','out');
    assert(loaded_flag,['BUG TRAP: could not load the data from ',agg_data])
    %unpack the data
    t_out=out.t_out;
    y_out=out.y_out;
    header=out.header;
  else
    %init records
    t=cell(1,numel(v.prefixes)); y=cell(1,numel(v.prefixes)); h=cell(1,numel(v.prefixes));
    %loop over all data files
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
        t{i}=out.t;
        y{i}=out.y;
        h{i}=out.header;
      else
        %need to make sure file.unwrap returned the txt file
        assert(file.isext(local_data,'.txt'),'BUG TRAP: expecting a txt file')
        %build the header info (GM, R and tide_system is assumed to be the same as CSR5x5)
        h{i}=struct(...
          'origin',local_data,...
          'raw',file.header(local_data,30),... %load the text (search for the end of the header up until line 30)
          'd',v.degrees(i),... 
          'o',v.orders{i}...
        );
        %branch on files with one or two coefficients
        if contains(h{i}.raw,'C21') || contains(h{i}.raw,'C22')
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
        t{i}=datetime('0000-01-01 00:00:00')+years(raw(:,1));
        %building data domain
        switch v.format
        case 'csr'       %NOTICE: this includes AOD mean for the solution period (aka correction, or 'corr')
          y_raw=raw(:,data_cols);
        case 'csr-grace' %NOTICE: this removes the AOD correction from SLR, making it more suitable to replace the GRACE C20
          y_raw=raw(:,data_cols)-raw(:,corr_cols)*1e-10;
        case 'csr-corr'  %NOTICE: this is the AOD correction
          y_raw=raw(:,corr_cols)*1e-10;
        case 'csr-sigma' %NOTICE: this is the solution sigma
          y_raw=raw(:,sigm_cols)*1e-10;
        end
        %building aggregated records
        for j=1:numel(v.orders{i}) %NOTICE: this is why v.degrees{i} should be scalar
          d=v.degrees(i);
          o=v.orders{i}(j);
          y{i}(:,gravity.colidx(d,o,lmax))=y_raw(:,j);
        end
        %save the data in mat format
        file.save_mat(struct('t',t{i},'y',y{i},'header',h{i}),local_data,'data_var','out')
      end
    end  
    %aggregate the data from the separate files, define header
    header=struct(...
      'descriptor','UT/CSR monthly 2x2 RL-06 time series from SLR',...
      'GM',3.986004415E+14,...
      'R',6.378136300E+06,...
      'raw',{cell(1,numel(h))},... %NOTICE: this {} syntax is needed so that header is not a vector
      'origin','',...
      'lmax',lmax,...
      'tide_system','zero_tide'...
    );
    %init y_out
    y_out=zeros(numel(t{1}),gravity.y_length(lmax));
    for i=1:numel(v.prefixes)
      %retrieve/sanitize time
      if i==1
        t_out=t{i};
      else
        assert(~any(~simpletimeseries.ist('==',t{i},t_out)),'time domain inconsistency')
      end
      %aggregate coefficients
      for j=1:numel(v.orders{i}) %NOTICE: this is why v.degrees{i} should be scalar
        colidx=gravity.colidx(h{i}.d,h{i}.o(j),lmax);
        y_out(:,colidx)=y{i}(:,colidx); 
      end
      %save header info
      header.raw{i}=h{i}.raw;
      header.origin=strjoin({header.origin,h{i}.origin},'; ');
    end
    %save the data in mat format
    file.save_mat(struct('t_out',t_out,'y_out',y_out,'header',header),agg_data,'data_var','out')
  end
end
function [t_out,y_out,header,y_out_error,y_out_AOD]=import_CSR5x5(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   slr.dir('CSR5x5'),@ischar;...
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
    header     =out.header;
  else
    %need to make sure file.unwrap returned the txt file
    assert(file.isext(local_data,'.txt'),'BUG TRAP: expecting a txt file')
    %declare header structure
    header=struct(...
      'GM',0,...
      'R',0,...
      'lmax',v.lmax,...
      'tide_system','unknown',...
      'descriptor','unknown',...
      'origin',local_data,...
      'static',[],...
      'labels',{{}},...
      'idx',struct([]),...
      'scale',1 ...
    );
    %define known details
    header.descriptor='UT/CSR monthly 5x5 gravity harmonics';
    %open the file
    fid=file.open(local_data);
    % Read header
    while true
      s=fgets(fid); 
      if contains(s,'end of header')
        break
      end
      if contains(s,'earth_gravity_constant')
        header.GM = str2double(strtrim(strrep(s,'earth_gravity_constant','')));
      end
      if (contains(s, 'radius'))
        header.R=str2double(strtrim(strrep(s,'radius','')));
      end
      if (contains(s, 'tide_system'))
        header.tide_system=strtrim(strrep(s,'tide_system',''));
      end
      if (contains(s, 'Units'))
        header.scale=str2double(strtrim(strrep(s,'Units','')));
      end
      if (contains(s, 'Coefficients:'))
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
      if (contains(s, '===================='))
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
    %init output
    y_out=zeros(1,gravity.y_length(header.lmax));
    y_out_error=y_out;y_out_AOD=y_out;
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
        y_out(      arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Cnm);    
        y_out_error(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Csigma); 
        y_out_AOD(  arc,gravity.colidx(d,o,header.lmax))=s(header.idx.CnmAOD); 
        if o==0, continue;end
        %save sine coefficient, error and value with AOD
        o=-o;
        y_out(      arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Snm);    
        y_out_error(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.Ssigma); 
        y_out_AOD(  arc,gravity.colidx(d,o,header.lmax))=s(header.idx.SnmAOD); 
      otherwise
        disp(['WARNING: ignoring line: ',strjoin(s,' ')])
      end
    end
    fclose(fid);
    %add static signal and error
    scale=header.scale*ones(size(y_out,1),1);
    y_out      =      scale.*y_out           +  static_signal;
    y_out_AOD  =      scale.*y_out_AOD       +  static_signal;
    y_out_error=sqrt((scale.*y_out_error).^2 +  static_error.^2);
    %save the data in mat format
    file.save_mat(struct('t_out',t_out,'y_out',y_out,'y_out_AOD',y_out_AOD,'y_out_error',y_out_error,'header',header),local_data,'data_var','out')
  end
end
function [t_out,y_out,header]=import_GSFC5x5(varargin)
  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'import_dir',   slr.dir('GSFC5x5'),@ischar;...
      'data_dir_url', 'https://earth.gsfc.nasa.gov/sites/default/files/2021-12',@ischar;...
      'data_file',    'GSFC_SLR_5x5c61s61.txt',@ischar;...
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
    header=out.header;
  else
    %need to make sure file.unwrap returned the txt file
    assert(file.isext(local_data,'.txt'),'BUG TRAP: expecting a txt file')
    %declare header structure
    header=struct(...
      'GM',0,...
      'R',0,...
      'lmax',v.lmax,...
      'tide_system','unknown',...
      'descriptor','unknown',...
      'origin',local_data,...
      'labels',{{}},...
      'idx',struct([]),...
      'scale',1 ...
    );
    %define known details
    header.descriptor='NASA GSFC SLR 5x5+C61/S61 time variable gravity';
    %open the file
    fid=file.open(local_data);
    % Read header
    while true
      s=fgets(fid); 
      if contains(s,'end of header')
        break
      end
      if contains(s,'Title: ')
        header.descriptor = strtrim(strrep(s,'Title: ',''));
      end
      if contains(s,'GM:')
        l=strsplit(s);
        header.GM = str2double(l{3});
      end
      if (contains(s, 'R:'))
        l=strsplit(s);
        header.R=str2double(l{3});
      end
      if (contains(s, 'C20 is'))
        header.tide_system=strtrim(strrep(s,'C20 is',''));
      end
      if (contains(s, 'Coefficient lines:'))
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
      if contains(s, 'Product:')
        break      
      end
    end
    %init loop variables
    arc=0; y_out=zeros(1,gravity.y_length(header.lmax));
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
        t_out(arc)=time.ToDateTime(s(1),'modifiedjuliandate'); %#ok<AGROW>
      case 4
        %get degree and order
        d=s(header.idx.n);
        o=s(header.idx.m);
        %skip if this degree is above the requested lmax
        if d>header.lmax; continue; end
        %save cosine coefficient
        y_out(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.C);    
        if o==0, continue;end
        %save sine coefficient
        o=-o;
        y_out(arc,gravity.colidx(d,o,header.lmax))=s(header.idx.S);    
      otherwise
        disp(['WARNING: ignoring line: ',strjoin(s,' ')])
      end
    end
    fclose(fid);
    %save the data in mat format
    file.save_mat(struct('t_out',t_out,'y_out',y_out,'header',header),local_data,'data_var','out')
  end
end
function [t_out,y_out,header]=import_C20(varargin)

%TODO: there is a 1 month timeshift between:
% - GSFC-7DAY and GSFC5x5
% - CSR2x2 and CSR-RL06
% - TN-07 and TN-11 (half a month?)

  % add input arguments and metadata to collection of parameters 'v'
  v=varargs.wrap('sources',{...
    {...
      'source'       , 'TN-14',@ischar;...
      'import_dir'   , slr.dir('TN-14'),@ischar;...
      'data_dir_url' ,      '',@ischar;...
      'data_file'    ,      '',@ischar;...
      'time_column'  ,       1,@isnumeric;...
      'signal_column',       3,@isnumeric;...
      'error_column' ,       5,@isnumeric;...
      'comment_style',     '*',@iscahr;...
      'data_format'  , '%7.1f%10.4f%22.13f%8.4f%8.4f',@ischar;...
    },...
  },varargin{:});
  %upper-case version name 
  v.source=upper(v.source);
  %update local import dir
  v.import_dir=slr.dir(v.source);
  %parse dependent arguments (can be over-written)
  switch v.source
    case 'TN-07'
      v.data_file=[v.source,'_C20_SLR.txt'];
%       v.data_dir_url='ftp://podaac.jpl.nasa.gov/allData/grace/docs/'; %NOTICE: this no longer works
      v.data_dir_url='none';
    case 'TN-11'
      v.data_file=[v.source,'_C20_SLR.txt'];
%       v.data_dir_url=['https://podaac-w10n.jpl.nasa.gov/w10n/allData/grace/docs/',v.data_file];
      v.data_dir_url='https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/';
      v.data_format='%7.1f%10.4f%22.13f%8.4f%8.4f%8.1f%10.4f';
      v.time_column=[1 6];
    case 'CSR-RL06'
      v.data_file='C20_RL06.txt';
      v.data_dir_url='http://download.csr.utexas.edu/pub/slr/degree_2/';
      v.data_format='%10.4f%19.10f%8.4f%8.4f%8.4f%16.4f%16.4f';
      v.comment_style='#';
      v.signal_column=2;
      v.error_column=4;
    case 'TN-14'
      v.data_file=[v.source,'_C30_C20_GSFC_SLR.txt'];
      v.data_dir_url='https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/';
      v.data_format='%7.1f%10.4f%22.13f%8.4f%8.4f%22.13f%8.4f%8.4f%8.1f%10.4f';
    case 'GSFC'
      v.data_file='GSFC_SLR_C20_GSM_replacement.txt';
      v.data_dir_url='https://earth.gsfc.nasa.gov/sites/default/files/neptune/images_db/';
    case 'GSFC-7DAY'
      v.data_file='GSFC_SLR_C20_7day.txt';
      v.data_dir_url='none';
      v.data_format='%10.4f%23.13f';
      v.time_column=1;
      v.signal_column=2;
      v.error_column=0;
    otherwise
      warning([...
        'Loading ',v.source,' C20 timeseries with externally-defined details:',newline,...
        'import_dir   :',v.import_dir,newline,...
        'data_dir_url :',v.data_dir_url,newline,...
        'data_file    :',v.data_file,newline,...
        'time_column  :',v.time_column,newline, ...
        'signal_column:',v.signal_column,newline,...
        'error_column :',v.error_column,newline,...
        'comment_style:',v.comment_style,newline,...
        'data_format  :',v.data_format,newline,...
      ])
  end
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
    header=out.header;
  else
    %default header
    %NOTICE: most of this info is assumed, becase there's no info on the data files
    header=struct(...
      'descriptor',v.source,...
      'GM',3.986004415E+14,...
      'R',6.378136300E+06,...
      'lmax',2,...
      'origin',local_data,...
      'tide_system','zero_tide'... 
    );
    %open the file
    fid=file.open(local_data);
    %skip headers, except for some sources which textscan can handle directly
    while ~contains(v.source,{'GSFC-7DAY','CSR-RL06'})
      if contains(fgetl(fid),{'Product:','PRODUCT:'})
        break
      end
    end
    %read the data
    dat=textscan(fid,v.data_format,'CommentStyle',v.comment_style);
    %close the file
    fclose(fid);
    % outputs 
    switch v.source
    case {'GSFC-7DAY','CSR-RL06'}
      t_out=years(mean([dat{v.time_column}],2))+datetime('0000-01-01 00:00:00');
    otherwise
      t_out=datetime(mean([dat{v.time_column}],2),'ConvertFrom','modifiedjuliandate');
%       %plot difference between MJD and years columns
%       t_alt=years(mean([dat{v.time_column+1}],2))+datetime('0000-01-01 00:00:00');
%       plotting.figure;
%       plot(t_out-t_alt)
%       plotting.enforce('plot_title',v.source);
%       keyboard
    end
    %init outputs
    y_out=zeros(numel(t_out),gravity.y_length(header.lmax));
    %propagate data
    y_out(:,gravity.colidx(2,0,header.lmax))=dat{v.signal_column};
    %inform
    str.say('start/stop of',header.descriptor,':',t_out(1),t_out(end))
    %save the data in mat format
    file.save_mat(struct('t_out',t_out,'y_out',y_out,'header',header),local_data,'data_var','out')

%     %TODO: implement importing errors
%     if e_idx>0
%       e=dat{v.error_column}(idx)*1e-10;
%     else
%       e=zeros(size(s));
%     end
  
  end
end

%% Auxiliarly data
function [t,s,e,d]=GetGRACEC20(varargin)
  %parse arguments that are required later
  v=varargs.wrap('sources',{...
    {...
      'source',     'TN-14', @ischar;...
      'mode',         'read', @ischar;...
      'start',time.zero_date, @isdatetime;...
      'stop',  time.inf_date, @isdatetime;...
      'data_dir', file.orbdir('auxiliary'),  @(i) ischar(i) && exist(i,'dir')~=0;
    },...
  },varargin{:});
  %some default parameters
  t_idx=1;
  s_idx=3;
  e_idx=5;
  CommentStyle='*';
  datfmt='%7.1f%10.4f%22.13f%8.4f%8.4f';
  %upper-case version name 
  v.source=upper(v.source);
  %parse dependent arguments (can be over-written)
  %(NOTICE: upper is redundant but the preprocessor shows non-capitalized cases)
  switch upper(v.source)
    case 'TN-07'
      datfil=[v.source,'_C20_SLR.txt'];
      daturl=['ftp://podaac.jpl.nasa.gov/allData/grace/docs/',datfil]; %NOTICE: this no longer works
    case 'TN-11'
      datfil=[v.source,'_C20_SLR.txt'];
%       daturl=['https://podaac-w10n.jpl.nasa.gov/w10n/allData/grace/docs/',datfil];
      daturl=['https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/',datfil];
      datfmt='%7.1f%10.4f%22.13f%8.4f%8.4f%8.1f%10.4f';
      t_idx=[1 6];
    case 'CSR-RL06'
      datfil='C20_RL06.txt';
      daturl=['http://download.csr.utexas.edu/pub/slr/degree_2/',datfil];
      datfmt='%10.4f%19.10f%8.4f%8.4f%8.4f%16.4f%16.4f';
      CommentStyle='#';
      s_idx=2;
      e_idx=4;
    case 'TN-14'
      datfil=[v.source,'_C30_C20_GSFC_SLR.txt'];
      daturl=['https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/',datfil];
      datfmt='%7.1f%10.4f%22.13f%8.4f%8.4f%22.13f%8.4f%8.4f%8.1f%10.4f';
    case 'GSFC'
      datfil='GSFC_SLR_C20_GSM_replacement.txt';
      daturl=['https://earth.gsfc.nasa.gov/sites/default/files/neptune/images_db/',datfil];
    case 'GSFC-7DAY'
      datfil='GSFC_SLR_C20_7day.txt';
      daturl='personal communication: email from bryant.d.loomis@nasa.gov';
      datfmt='%10.4f%23.13f';
      t_idx=1;
      s_idx=2;
      e_idx=0;
    otherwise
      error(['Cannot handle version ''',v.source,'''.'])
  end
  %append v.file and v.url
  v=varargs.wrap('sources',{v,{...
        'file',fullfile(v.data_dir,datfil),  @ischar;...
        'url',daturl,                        @ischar;...
      }},varargin{:});
  switch lower(v.mode)
  case 'data_file'
    t=v.file;
  case 'read' %'set' if not already, then 'get'
    if ~file.exist(v.file)
      GetGRACEC20('mode','set','version',v.source,'data_dir',v.data_dir);
    end
    [t,s,e,d]=GetGRACEC20('mode','get','version',v.source,'file',v.file,'start',v.start,'stop',v.stop);
  case 'reload' %'set' then 'get'
    GetGRACEC20('mode','set','version',v.source,'data_dir',v.data_dir);
    [t,s,e,d]=GetGRACEC20('mode','get','version',v.source,'file',v.file,'start',v.start,'stop',v.stop);
  case 'plot'
    %NOTICE: this is a low-level plot, without much features
    [t,s,e,d]=GetGRACEC20(varargin{:},'mode','read');
    plotting.figure;
    plot(t ,s ,'o-','MarkerSize',4), hold on
    plotting.enforce('plot_title',d,'plot_legend_location','none');
  case 'get' %read the downloaded data
    %open the file
    fid=file.open(v.file);
    %sanity
    if fid<=0
      error([mfilename,': cannot open the data file ''',v.file,'''.'])
    end
    %default descriptor
    d=['C20 time series, version ',v.source];
    %get header info (NOTICE: upper is redundant but the preprocessor shows non-capitalized cases)
    switch upper(v.source)
    case {'GSFC-7DAY','CSR-RL06'}
      found=true;
    otherwise
      found=false;
    end
    c=0;
    while ~found
      c=c+1;
      line=fgetl(fid);
      if line<0
        error(['Could not find keyword PRODUCT in data file ',v.file])
      end
      found=any(str.contains(line,{'PRODUCT:','Product:'}));
      %these strings should come in the same order as here
      if any(str.contains(line,{'TITLE:','Title:'}))
        d=strtrim(str.clean(line,{'TITLE:','Title:'}));
      end
      if str.contains(line,'UPDATE HISTORY:')
        d=strjoin({d,strtrim(strrep(line,'UPDATE HISTORY:',''))},newline);
      end
      if str.contains(line,'Data span:')
        d=strjoin({d,strtrim(strrep(line,'Data span:',''))},newline);
      end
    end
%     str.say('read through ',c,'header lines')
    %read the data
    dat=textscan(fid,datfmt,'CommentStyle',CommentStyle);
    %close the file
    fclose(fid);
    % outputs (NOTICE: upper is redundant but the preprocessor shows non-capitalized cases)
    switch upper(v.source)
    case {'GSFC-7DAY','CSR-RL06'}
      t=years(mean([dat{t_idx}],2))+datetime('0000-01-01 00:00:00');
    otherwise
      t=datetime(mean([dat{t_idx}],2),'ConvertFrom','modifiedjuliandate');
    end
    % enforce start/stop
    idx=t>=v.start & t<=v.stop;
    t=t(idx);
    %inform
    str.say('start/stop of',d,':',t(1),t(end))
    s=dat{s_idx}(idx);
    if e_idx>0
      e=dat{e_idx}(idx)*1e-10;
    else
      e=zeros(size(s));
    end
  case 'set' %download the data
    if url.is(v.url)
      websave(v.file,v.url);
      t=file.strload(v.file);
      disp(t)
    else
      t=[];
    end
    s=[];e=[];d=[];
  otherwise
    error([mfilename,': unknown mode ''',v.mode,'''.'])
  end
end

%% Outdate functions
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
  t=datetime('0000-01-01 00:00:00')+years(raw(:,2));
  %building object
  ts=simpletimeseries(t,raw(:,data_cols)*1e-10,...
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