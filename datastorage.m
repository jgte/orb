classdef datastorage
  %static
  properties(Constant)
    
%         'fields',struct(...
%           'AC0X'  ,1e-7,...
%           'AC0Z'  ,1e-6,...
%           'AC0Y1' ,1e-5,...
%           'AC0Y2' ,1e-5,...
%           'AC0Y3' ,1e-5,...
%           'AC0Y4' ,2e-5,...
%           'AC0Y5' ,2e-5,...
%           'AC0Y6' ,5e-5,...
%           'AC0Y7', 5e-5,...
%           'AC0Y8', 5e-5,...
%           'AC0XD', 5e-7,...
%           'AC0ZD', 5e-6,...
%           'AC0YD1',5e-5,...
%           'AC0YD2',5e-5,...
%           'AC0YD3',1e-4,...
%           'AC0YD4',1e-4,...
%           'AC0YD5',1e-4,...
%           'AC0YD6',5e-3,...
%           'AC0YD7',5e-3,...
%           'AC0YD8',5e-4,...
%           'AC0XQ' ,3e-7,...
%           'AC0ZQ' ,5e-6,...
%           'AC0YQ1',5e-3,...
%           'AC0YQ2',5e-3,...
%           'AC0YQ3',5e-3,...
%           'AC0YQ4',5e-3,...
%           'AC0YQ5',5e-3,...
%           'AC0YQ6',5e-3,...
%           'AC0YQ7',5e-3,...
%           'AC0YQ8',5e-4 ...
%         )...
        
    parts={...
      'type',...
      'level',...
      'field',...
      'sat'...
    };

    %default value of some internal parameters
    default_list=struct(...
      'par_calpar_csr',struct(...
        'longtermbias',struct(...
          'A',fullfile('input','bsA2003'),...
          'B',fullfile('input','bsB2003')...
        ),...
        'coord',{{'X','Y','Z'}},...
        'data_col_idx',1,...
        'levels',struct(...
          'aak',   1,...
          'accatt',1,...
          'estim', 1 ...
        ),...
        'fields',struct(...
          'AC0X'  ,1e-7,...
          'AC0XD', 5e-7,...
          'AC0XQ' ,3e-7,...
          'AC0Z'  ,1e-6,...
          'AC0ZD', 5e-6,...
          'AC0ZQ' ,5e-6...
        )...
      ),...
      'par_acc',struct(...
        'data_col_name',{{'X','Y','Z'}},...
        'data_col_idx',[1 3]...
      ),...
      'par_data',struct(...
        'dir','data'...
      ),...
      'par_plot',struct(...
        'dir','plots',...
        'prefix','',...
        'suffix','',...
        'xlimits',[-inf,inf],...
        'ylimits',[-inf,inf],...
        'size',200+[0,0,21,9]*50,...
        'units','points',...
        'visible','on',...
        'fontsize',struct(...
          'axis', 24,...
          'title',32,...
          'label',28),...
        'line',struct(...
          'width',2)...
      )...
    );
%     sats={'A','B'};
    parameter_list=struct(...
      'par_modes',      struct('default',grace.default_list.par_modes,     'validation',@(i) isstruct(i)),...
      'par_calpar_csr', struct('default',grace.default_list.par_calpar_csr,'validation',@(i) isstruct(i)),...
      'par_acc',        struct('default',grace.default_list.par_acc,       'validation',@(i) isstruct(i)),...
      'par_data',       struct('default',grace.default_list.par_data,      'validation',@(i) isstruct(i)),...
      'par_plot',       struct('default',grace.default_list.par_plot,      'validation',@(i) isstruct(i)),...
      'start',          struct('default',datetime([0 0 31]),               'validation',@(i) isdatetime(i)),...
      'stop',           struct('default',datetime([0 0 31]),               'validation',@(i) isdatetime(i))...
    );
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={};
  end
  %read only
  properties(SetAccess=private)
    category
    par
    data
  end
  %private (visible only to this object)
  properties(GetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
    start
    stop
  end
  methods(Static)
    function out=parameters
      out=fieldnames(grace.parameter_list);
    end
  end
  methods
    %% constructor
    function obj=datastorage(varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      %declare parameters
      for j=1:numel(grace.parameters)
        %shorter names
        pn=grace.parameters{j};
        %declare parameters
        p.addParameter(pn,grace.parameter_list.(pn).default,grace.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(varargin{:});
      % reset data type list
      obj=obj.datatype_init;
      % save parameters with defaults first, they may be needed below
      for i=1:numel(grace.parameters)
        %shorter names
        pn=grace.parameters{i};
        %parameters have a dedicated member
        if strcmp(pn(1:4),'par_')
          obj.par.(pn(5:end))=collapse(p.Results.(pn));
        else
          obj.(pn)=collapse(p.Results.(pn));
        end
      end
      % data is added to this initialized object with the 'init' method (see below)
    end
    %% dataname handling
    function out=fix(obj,dataname)
      % The object dataname includes the 'category' field, which is common to each
      % datastorage instance. When calling the datatype_*/level_*/field_*/sat_* members
      % (either directly or through data_*/vector_* members), the output of 
      % dataname.cells_clean needs to be stripped of the 'category'. 
      if ~strcmp(dataname.category,obj.category)
        error([mfilename,...
          ': requesting a product of category ''',dataname.category,''', ',...
          'while this object is of category ''',obj.category,'''.'...
        ]);
      end
      out=dataname.cells;
      out=out(2:end);
    end
    %% datatype operations
    function obj=datatype_init(obj)
      obj.data=struct([]);
    end
    function obj=datatype_set(obj,datatype_name,datatype_value)
      obj.data(1).(datatype_name)=datatype_value;
    end
    function out=datatype_get(obj,datatype_name)
      if isfield(obj.data,datatype_name)
        out=obj.data.(datatype_name);
      else
        out=[];
      end
    end
    function out=datatype_list(obj)
      out=fieldnames(obj.data);
    end
    %% level operations
    function obj=level_init(obj,datatype)
      obj=obj.datatype_set(datatype,struct([]));
    end
    function obj=level_set(obj,datatype,level_name,level_value)
      datatype_value=obj.datatype_get(datatype);
      datatype_value.(level_name)=level_value;
      obj=obj.datatype_set(datatype,datatype_value);
    end
    function out=level_get(obj,datatype,level_name)
      datatype_value=obj.datatype_get(datatype);
      if isfield(datatype_value,level_name)
        out=datatype_value.(level_name);
      else
        out=[];
      end
    end
    function out=level_list(obj,datatype)
      out=obj.datatype_get(datatype);
      if ~isempty(out)
        out=fieldnames(out);
      end
    end
    %% field operations
    function obj=field_init(obj,datatype,level)
      obj=obj.level_set(datatype,level,struct([]));
    end
    function obj=field_set(obj,datatype,level,field_name,field_value)
      level_value=obj.level_get(datatype,level);
      level_value.(field_name)=field_value;
      obj=obj.level_set(datatype,level,level_value);
    end
    function out=field_get(obj,datatype,level,field_name)
      level_value=obj.level_get(datatype,level);
      if isfield(level_value,field_name)
        out=level_value.(field_name);
      else
        out=[];
      end
    end
    function out=field_list(obj,datatype,level)
      out=obj.level_get(datatype,level);
      if ~isempty(out)
        out=fieldnames(out);
      end
    end
    %% sat operations
    function obj=sat_init(obj,datatype,level,field)
      obj=obj.field_set(datatype,level,field,struct([]));
    end
    function obj=sat_set(obj,datatype,level,field,sat_name,sat_value)
      field_value=obj.field_get(datatype,level,field);
      field_value.(sat_name)=sat_value;
      obj=obj.field_set(datatype,level,field,field_value);
    end
    function out=sat_get(obj,datatype,level,field,sat_name)
      field_value=obj.field_get(datatype,level,field);
      if isfield(field_value,sat_name)
        out=field_value.(sat_name);
      else
        out=[];
      end
    end
    function out=sat_list(obj,datatype,level,field)
      out=obj.field_get(datatype,level,field);
      if ~isempty(out)
        out=fieldnames(out);
      end
    end
    %% generalized datanames operations
    function out=data_function(~,function_name,dataname)
      if ~isa(dataname,'datanames')
        error([mfilename,': input ''names'' must be of class ''datanames'', not a ',class(dataname),'.'])
      end
      if strcmp(function_name,'list')
        out=[dataname.leaf_type(1),'_',function_name];
      else
        out=[dataname.leaf_type,'_',function_name];
      end
    end  
    function obj=data_init(obj,dataname)
      dnf=obj.fix(dataname);
      obj=obj.(obj.data_function('init',dataname))(dnf{:});
    end
    function obj=data_set(obj,dataname,value)
      dnf=obj.fix(dataname);
      obj=obj.(obj.data_function('set',dataname))(dnf{:},value);
    end
    function out=data_get(obj,dataname)
      dnf=obj.fix(dataname);
      out=obj.(obj.data_function('get',dataname))(dnf{:});
    end
    function out=data_list(obj,dataname)
      dnf=obj.fix(dataname);
      out=obj.(obj.data_function('list',dataname))(dnf{:});
    end
    %% cell array operations
    function out=vector_names_sanity(~,in,name,default)
      if isempty(in)
        %assign default values
        out=default;
      else
        if ~iscellstr(in)
          if ischar(in)
            %if input is a string, convert it to cell
            out={in};
          else
            %cannot handle anything other than cell strings and strings
            error([mfilename,': input ''',name,''' must be of class ''cell'', not ''',class(in),'''.'])
          end
        else
          %we want cell strings
          out=in;
        end
      end
    end
    function out=vector_names_isvalid(obj,datatype,level,field,sat)
      try
        out=ischar(class(obj.data.(datatype).(level).(field).(sat)));
      catch
        out=false;
      end
    end
    function dataname_list=vector_names(obj,datatype,level,field,sat)
      %input arguments datatype, level, field and sat are either:
      % - a cell array of strings;
      % - a string (converted to a scalar cell array);
      % - empty (converted to all possible values of that type of argument).
      % This routine builds a list with all combinations of the given 
      % datatype, level, field and sat (after converting as mentioned above).
      if ~exist('sat',     'var'),      sat='';end
      if ~exist('field',   'var'),    field='';end
      if ~exist('level',   'var'),    level='';end
      if ~exist('datatype','var'), datatype='';end
      %if there's only one input, it may of datanames class, which is a compact way
      %to represent the inputs datatype,level,field,sat.
      if isa(datatype,'datanames')
        %recursive call
        dnf=obj.fix(datatype);
        dataname_list=obj.vector_names(dnf{:});
      else
        % datatype is root, so it does not depend on any other 
        datatype=obj.vector_names_sanity(datatype,'datatype',obj.datatype_list);
        %make room for outputs (this is just a guess, because level and field can be empty)
        dataname_list=cell(numel(datatype)*numel(level)*numel(field)*numel(sat),1);
        c=0;
        %loop over all datatypes, levels, fields and sats
        for i=1:numel(datatype)
          level_now=obj.vector_names_sanity(level,'level',obj.level_list(datatype{i}));
          for j=1:numel(level_now)
            field_now=obj.vector_names_sanity(field,'field',obj.field_list(datatype{i},level_now{j}));
            for k=1:numel(field_now)
              sat_now=obj.vector_names_sanity(sat,'sat',obj.sat_list(datatype{i},level_now{j},field_now{k}));
              for l=1:numel(sat_now)
                if obj.vector_names_isvalid(datatype{i},level_now{j},field_now{k},sat_now{l})
                  c=c+1;
                  dataname_list{c}=datanames({obj.category,datatype{i},level_now{j},field_now{k},sat_now{l}});
                end
              end
            end
          end
        end
        %remove empty cells
        dataname_list=dataname_list(cellfun(@(i)(~isempty(i)),dataname_list));
      end
    end
    function values=vector_get(obj,dataname_list)
      %make room for outputs
      values=cell(size(dataname_list));
      %populate outputs
      for i=1:numel(values)
        values{i}=obj.data_get(dataname_list{i});
      end
    end
    function obj=vector_set(obj,dataname_list,values)
      %sanity
      if numel(dataname_list) ~= numel(values)
        error([mfilename,': number of elements of input ''values'' (',num2str(numel(values)),') ',...
          'must be the same as in ''dataname_list'' (',num2str(numel(dataname_list)),').'])
      end
      %propagate inputs
      for i=1:numel(values)
        obj=obj.data_set(dataname_list{i},values{i});
      end
    end
    % Applies the function f to the cell array given by vector_get(dataname_list) and
    % returns the result, so it transforms (tr) a list of objects into something else
    function out=vector_tr(obj,f,dataname_list,varargin)
      %collapse requested data into cell array and operate on the cell array
      out=f(obj.vector_get(dataname_list),varargin{:});
    end
    % Applies the function f to the cell array given by vector_get(dataname_list) and
    % saves the resulting objects back to the same set, using vector_set(dataname_list,...)
    function obj=vector_op(obj,f,dataname_list,varargin)
      %operate on the requested set
      values=obj.vector_tr(f,dataname_list,varargin{:});
      %sanity
      if ~iscell(values)
        error([mfilename,': function ',fun2str(f),' must return a cell array, not a ',class(values),'.'])
      end
      %propagate cell array back to object
      obj=obj.vector_set(dataname_list,values);
    end
    % Retrieves the values of particular field or zero-input method from the 
    % simpletimeseries (sts) stored at the leaves defined by the input names.
    function out=vector_sts(obj,sts_field,dataname_list)
      out=cell(size(dataname_list));
      obj_list=obj.vector_get(dataname_list);
      for i=1:numel(out)
        out{i}=obj_list{i}.(sts_field);
      end
    end
    %% start/stop operations
    function out=get.start(obj)
      out=obj.vector_sts('start',obj.vector_names);
      out=min([out{:}]);
    end
    function out=get.stop(obj)
      out=obj.vector_sts('stop',obj.vector_names);
      out=max([out{:}]);
    end
    function obj=set.start(obj,start)
      %retrieve cell names
      names=obj.vector_names;
      %get cell array with simpletimeseries objects
      values=obj.vector_get(names);
      %loop over complete set of objects
      for i=1:numel(values)
        values{i}.start=start;
      end
      %back-propagate modified set of objects
      obj=obj.vector_set(names,values);
    end
    function obj=set.stop(obj,stop)
      %retrieve cell names
      names=obj.vector_names;
      %get cell array with simpletimeseries objects
      values=obj.vector_get(names);
      %loop over complete set of objects
      for i=1:numel(values)
        values{i}.stop=stop;
      end
      %back-propagate modified set of objects
      obj=obj.vector_set(names,values);
    end
    %% metadata interface
    function obj=mdset(obj,dataname,varargin)
      obj.par.(dataname.name_clean)=dataproduct(dataname,varargin{:});
    end
    function md=mdget(obj,dataname)
      md=[];delta=0;
      while isempty(md) && delta<4
        dataname_now=dataname.trunk_name_clean(delta);
        if isfield(obj.par,dataname_now)
          md=obj.par.(dataname_now);
        else
          delta=delta+1;
        end
      end
      if isempty(md)
        error([mfilename,': cannot find the metadata of product ',dataname.name,' or any of its trunk elements.'])
      end
    end
    %% datatype initialization
    function obj=init(obj,dataname,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('dataname',@(i) ischar(i) || isa(i,'datanames'));      
      % parse it
      p.parse(dataname,varargin{:});
      % convert to object
      dataname=datanames(dataname);
      % save the cateogry, if not done already
      if isempty(obj.category)
        obj.category=dataname.category;
      else
        if ~strcmp(obj.category,dataname.category)
          error([...
            mfilename,': current object is already attributed to category ''',obj.category,...
            ''', cannot change it to category ''',dataname.category,'''.'...
          ])
        end
      end
      % load the metadata for this product
      obj=obj.mdset(dataname,varargin{:});
      % make sure all data sources are loaded
      obj=obj.init_sources(dataname,varargin{:});
      % invoke init method
      ih=str2func(obj.mdget(dataname).mdget('method'));
      obj=ih(obj,dataname,varargin{:});
    end
    function obj=init_sources(obj,dataname,varargin)
      %loop over all source data
      for i=1:obj.mdget(dataname).nr_sources
        %load this source
        obj=obj.init(obj.mdget(dataname).sources(i),varargin{:});
      end
    end
    function obj=init_nrtdm(obj,dataname,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('dataname',@(i) ischar(i) || isa(i,'datanames'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(dataname,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
      %retrieve product info
      product=obj.mdget(dataname);
      %retrieve relevant parameters
      sats         =product.mdget('sats'); 
      indir        =product.mdget('nrtdm_data_dir');
      nrtdm_sats   =product.mdget('nrtdm_sats'); 
      nrtdm_product=product.mdget('nrtdm_product');
      %sanity
      str.sizetrap(sats,nrtdm_sats)
      %loop over the satellites
      for s=1:numel(sats)
        %load the data
        tmp=nrtdm([nrtdm_sats{s},'_',nrtdm_product],p.Results.start,p.Results.stop,'data_dir',indir);
        %get the data
        obj=obj.sat_set(dataname.type,dataname.level,dataname.field,sats{s},...
          tmp.ts...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
    %% generalized plotting
    function h=justplot(obj,dn,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=obj.mdget(dn).plot_args(p,varargin{:});
      %sanity on non-optional parameters
      if ~isa(dn,'datanames') && ~isscalar(dn)
        error([mfilename,': can only handle input ''dn'' as scalars of class ''datanames'', not ''',class(in),'''.'])
      end
      %retrieve the requested data
      d=obj.data_get(dn);
      %checking data class
      switch class(d)
      case 'simpletimeseries'
        h=d.plot(...
          'columns', p.Results.plot_columns,...
          'outlier', p.Results.plot_outlier,...
          'zeromean',p.Results.plot_zeromean...
        );
        %save units
        h.y_units=d.y_units{p.Results.plot_columns{1}};
        for i=2:numel(p.Results.plot_columns)
          if ~strcmp(h.y_units,d.y_units(p.Results.plot_columns{i}))
            error([mfilename,':BUG TRAP: y-units are not consistent in all plotted lines.'])
          end
        end
      otherwise
        error([mfilename,': cannot plot data of class ',class(d),'; implementation needed!'])
      end
    end
    function h=plot(obj,dataname,varargin)
      %sanity on non-optional parameters
      if ~isa(dataname,'datanames') && ~ischar(dataname) && ~iscell(dataname)
        error([mfilename,': can not handle input ''dataname'' of class ''',class(dataname),'''.'])
      end
      %handle input dataname
      switch class(dataname)
      case 'cell'
        %propagate
        dataname_list=dataname;
      case 'datanames'        
        %propagate
        dataname_list={dataname};
      case 'char'
        % loop over all data in lower levels (in case this isn't already the lowest level)
        dataname_list=vector_names(obj,datanames(dataname));
        % sanity
        if numel(dataname_list)==0
          error([mfilename,': cannot find any data with dataname ',dataname,'.'])
        end
      otherwise
        error([mfilename,': can not handle input ''dataname'' of class ''',class(dataname),'''.'])
      end
      % recursive call multiple plots
      if numel(dataname_list)>1
        h=cell(size(dataname_list));
        for i=1:numel(dataname_list)
          h{i}=obj.plot(dataname_list{i},varargin{:});
        end
        return
      end
      %save the single entry in dataname
      dataname_now=dataname_list{1};
      %parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=obj.mdget(dataname).plot_args(p,varargin{:});
      %if columns are not to be plotted together, need to expand the calls to obj.plot to include each column
      if ~p.Results.plot_columns_together
        %retrieve plot columns indexes and names
        plot_columns=p.Results.plot_columns;
        col_names   =p.Results.plot_column_names;
        %make room for handles
        h=cell(1,numel(plot_columns));
        %loop over all data columns to plot
        for i=1:numel(plot_columns)
          %buid suffixes for the plots of this data column
          if ~isempty(p.Results.plot_file_suffix)
            file_suffix=[col_names(plot_columns{i}),{p.Results.plot_file_suffix}];
          else
            file_suffix=col_names(plot_columns{i});
          end
          if ~isempty(p.Results.plot_title_suffix)
            title_suffix=[col_names(plot_columns{i}),{p.Results.plot_title_suffix}];
          else
            title_suffix=col_names(plot_columns{i});
          end
          h{i}=obj.plot(dataname_now,...
            varargin{:},...
            'plot_file_suffix', strjoin(file_suffix, '.'),...
            'plot_title_suffix',strjoin(title_suffix,' '),...
            'plot_columns_together',true,...
            'plot_columns',plot_columns(i)...
          );
        end
      else
        %plot filename arguments
        filename_args=[obj.mdget(dataname_now).file_args('plot'),{...
          'start',obj.start,...
          'stop',obj.stop,...
          'timestamp',true,...
          'remove_part',p.Results.plot_together,...
          'prefix',p.Results.plot_file_prefix...
          'suffix',p.Results.plot_file_suffix...
        }];
        %plot filename
        filename=dataname_now.file(filename_args{:});
        % check if plot is already there
        if isempty(dir(filename))
          figure('visible',p.Results.plot_visible);
          %retrive data names to be plotted here
          %NOTICE: this will cause plots to appear multiple times if they are not saved
          dataname_list_to_plot=obj.vector_names(dataname_now.edit(p.Results.plot_together,''));
          h=cell(size(dataname_list_to_plot));
          %loop over all data to plot
          for i=1:numel(dataname_list_to_plot)
            h{i}=obj.justplot(dataname_list_to_plot{i},varargin{:});
          end
          % paranoid sanity
          str.sizetrap(h,dataname_list_to_plot)
          %enforce plot preferences
          obj.mdget(dataname_now).enforce_plot
          %add prefixes, if there
          for i=1:numel(datastorage.parts)
            %only add prefix to title if:
            % - 'plot_label_prefix_<part>' is defined (otherwise cannot know which prefix to add)
            if ~isempty(p.Results.(['plot_label_prefix_',datastorage.parts{i}]))
              %loop over all plotted datanames
              for j=1:numel(dataname_list_to_plot)
                old_part_value=dataname_list_to_plot{j}.(datastorage.parts{i});
                new_part_value=[...
                  p.Results.(['plot_label_prefix_',datastorage.parts{i}]),...
                  old_part_value...
                ];
                dataname_list_to_plot{j}=dataname_list_to_plot{j}.edit(...
                  datastorage.parts{i},...  %part name
                  new_part_value...         %part value
                );
              end
            end
          end
          %add legend if there are multiple lines
          if numel(h)>1
            if ~isempty(p.Results.plot_legend)
              %some sanity
              str.sizetrap(h,p.Results.plot_legend)
              %propagate
              legend_str=p.Results.plot_legend;
            else
              %get unique parts in datanames 
              prefixes=datanames.unique(dataname_list_to_plot);
              % paranoid sanity
              str.sizetrap(h,prefixes)
              %get how many lines have been plotted
              n=0;
              for j=1:numel(h)
                n=n+numel(h{j}.legend);
              end
              %make room for legend strings
              legend_str=cell(1,n);
              %loop over all legend entries
              c=0;
              for j=1:numel(h)
                for k=1:numel(h{j}.y_mean)
                  %define sufix
                  if p.Results.plot_zeromean
                    suffix=num2str(h{j}.y_mean{k});
                  end
                  %build legend stirngs
                  c=c+1;
                  legend_str{c}=str.clean(...
                    strjoin({prefixes{j},suffix},' '),...
                    {'title','succ_blanks'}...
                  );
                end
              end
            end
            set(legend(legend_str),'fontname','courier')
          end
          %get title parts
          title_str=datanames.common(dataname_list_to_plot);
          %suppress some parts, if requested
          title_str=setdiff(title_str,p.Results.plot_title_suppress,'stable');
          %clean up the tile and put it there
          title(str.clean(strjoin([title_str,{p.Results.plot_title_suffix}],' '),'_'))
          %fix y-axis label
          if ~isempty(p.Results.plot_ylabel)
            ylabel(p.Results.plot_ylabel)
          else
            for i=2:numel(h)
              if ~strcmp(h{1}.y_units,h{i}.y_units)
                error([mfilename,':BUG TRAP: y-units are not consistent in all plotted lines.'])
              end
            end
            ylabel(['[',h{1}.y_units,']'])
          end
          grid on
          %save this plot
          saveas(gcf,filename)
        else
          h=[];
        end
      end
    end
    function h=plot_mult(obj,dataname_list,plot_columns,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=obj.mdget(dn).plot_args(p,varargin{:});
      %sanity on non-optional parameters
      if ~iscell(dataname_list)
        error([mfilename,': can not handle input ''dataname'' of class ''',class(dataname_list),''', only cell.'])
      end
      if ~isnumeric(plot_columns) && ~isfinite(plot_columns)
        error([mfilename,': can not handle input ''data_indexes'' of class ''',class(plot_columns),''' or non-finite.'])
      end
      %make room for outputs
      h=cell(size(dataname_list));
      %loop over all datanames
      for i=1:numel(dataname_list)
        h{i}=obj.justplot(dataname_list{i},varargin{:},...
          'plot_columns_together',true,...
          'plot_columns',plot_columns(i)...
        );
      end
      keyboard
    end
    %% generalized operations
    function obj=stats(obj,dataname)
      %retrieve products info
      product=obj.mdget(dataname);
      sourcep=obj.mdget(product.sources(1));
      %paranoid sanity
      if product.nr_sources~=1
        error([mfilename,': number of sources in product ',dataname,...
          ' is expected to be 1, not ',num2str(product.nr_sources),'.'])
      end
      %check if data is already there
      if ~product.isfile('data')
        %get names of parameters and levels
        levels=sourcep.mdget('levels');
        fields=sourcep.mdget('fields');
        sats  =sourcep.mdget('sats');
        %loop over all data
        for i=1:numel(levels)
          for j=1:numel(fields)
            for s=1:numel(sats)
              %retrive data for this satellite
              d=obj.sat_get(sourcep.dataname.type,levels{i},fields{j},sats{s});
              %compute stats
              stats_data=d.stats(...
                'period',product.mdget('stats_period'),...
                'overlap',product.mdget('stats_overlap'),...
                'outlier',product.mdget('stats_outlier'),...
                'struct_fields',product.mdget('stats')...
              );
              %loop over all requested stats
              stats_fields=fieldnames(stats_data);
              for si=1:numel(stats_fields)
                %build implicit 'sat' names
                sat_name=[sats{s},'_',stats_fields{si}];
                %save stats in implicit 'sat' names
                obj=obj.sat_set(...
                  product.dataname.type,levels{i},fields{j},sat_name,...
                  stats_data.(stats_fields{si})...
                ); 
              end
            end
          end
        end
        %save data
        s=obj.datatype_get(product.dataname.type); %#ok<*NASGU>
        save(char(product.file('data')),'s');
        clear s
      else
        %load data
        load(char(product.file('data')),'s');
        levels=fieldnames(s); %#ok<NODEF>
        for i=1:numel(levels)
          obj=obj.level_set(product.dataname.type,levels{i},s.(levels{i}));
        end
      end
    end
    function obj=corr(obj,dataname)
      %retrieve products info
      product=obj.mdget(dataname);
      sourcep=obj.mdget(product.sources(1));
      %paranoid sanity
      if product.nr_sources~=1
        error([mfilename,': number of sources in product ',dataname,...
          ' is expected to be 1, not ',num2str(product.nr_sources),'.'])
      end
      %check if data is already there
      if ~product.isfile('data')
        %get names of parameters and levels
        levels=sourcep.mdget('levels');
        fields=sourcep.mdget('fields');
        sats  =sourcep.mdget('sats');
        %loop over all data
        for i=1:numel(levels)
          for j=1:numel(fields)
            %retrive data for both satellites
            sats=obj.field_get(sourcep.dataname.type,levels{i},fields{j});
            %get sat namees
            satnames=fieldnames(sats);
            %can only deal with two sats at the moment
            if numel(fieldnames(sats))~=2
              error([mfilename,': can only deal with two sats at the moment. Implementation needed!'])
            end
            %compute stats
            stats_data=simpletimeseries.stats2(...
              sats.(satnames{1}),...
              sats.(satnames{2}),...
              'period',product.mdget('period'),...
              'overlap',product.mdget('overlap'),...
              'outlier',product.mdget('outlier'),...
              'struct_fields',product.mdget('stats')...
            );
            %loop over all requested stats
            stats_fields=fieldnames(stats_data);
            for si=1:numel(stats_fields)
              %save stats in 'sat' names given by the requested stats name
              obj=obj.sat_set(...
                product.dataname.type,levels{i},fields{j},stats_fields{si},...
                stats_data.(stats_fields{si})...
              ); 
            end
          end
        end
        %save data
        s=obj.datatype_get(product.dataname.type); %#ok<*NASGU>
        save(char(product.file('data')),'s');
        clear s
      else
        %load data
        load(char(product.file('data')),'s');
        levels=fieldnames(s); %#ok<NODEF>
        for i=1:numel(levels)
          obj=obj.level_set(product.dataname.type,levels{i},s.(levels{i}));
        end
      end
    end
    %% costumized operations
    function out=data_op(obj,opname,datatype,level,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('outlier',  0, @(i) isfinite(i));
      % parse it
      p.parse(varargin{:});
      switch lower(datatype)

      case 'acc'
        switch lower(opname)
        case 'load-stats'
          %compute residuals relative to the available models
          models=obj.field_list('acc','mod');
          for m=1:numel(models)
            filenames=simpletimeseries.unwrap_datafiles(...
              obj.id('data','accres_csr',level,models{m},'DATE_PLACE_HOLDER','stats'),...
              'start',obj.start,...
              'stop' ,obj.stop,...
              'only_existing',false...
            );
            for i=1:numel(filenames)
              if isempty(dir(filenames{i}))
                for s=1:numel(grace.sats)
                  %get the l1b acc data
                  l1bcsr=obj.sat_get('acc','l1b','csr',grace.sats{s});
                  %get the calibration model
                  calmod=obj.sat_get('acc','calmod',level,grace.sats{s});
                  %get the modeled acc
                  accmod=obj.sat_get('acc','mod',models{m},grace.sats{s});
                  if any([...
                    ~isa(l1bcsr,'simpletimeseries'),...
                    ~isa(calmod,'simpletimeseries'),...
                    ~isa(accmod,'simpletimeseries')...
                  ])
                    %patch nan residuals
                    accres=simpletimeseries(...
                      [obj.start;obj.stop],...
                      nan(2,numel(obj.par.acc.data_col_name))...
                    );
                  else
                    %calibrate the accelerometer
                    acccal=l1bcsr+calmod;
                    %compute residual
                    accres=acccal-accmod;
                  end
                  %save descriptor
                  accres.descriptor=str.clean(filenames{i},{'file','.'});
                  %compute stats per period
                  s.(grace.sats{s})=accres.stats(...
                    'period',days(1),...
                    'overlap',seconds(0),...
                    'outlier',p.Results.outlier,...
                    'struct_fields',obj.par.modes.stats...
                  );
                end
                save(filenames{i},'s')
              else
                load(filenames{i})
              end
            end
            %save this data to output var
            out.(models{m})=s;
          end
        case 'load-stats2'
          error([mfilename,': implementation needed.'])
        otherwise
          error([mfilename,': cannot handle operation ''',opname,''' for datatype ''',datatype,'''.'])
        end
      otherwise
        error([mfilename,': cannot handle datatype ''',datatype,'''.'])
      end
    end
    %% costumized plotting routine
    function oldplot(obj,plotname,datatype,level,varargin)
      switch lower(datatype)
      case 'calpar_csr'
        switch lower(plotname)
        case ''
          % loop over CSR cal pars
          fields=obj.field_list(datatype,level);
          for j=1:numel(fields)
            filename=obj.id('plot',datatype,level,fields{j},'');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %retrive data
              sat=struct(...
                'A',obj.sat_get(datatype,level,fields{j},'A'),...
                'B',obj.sat_get(datatype,level,fields{j},'B')...
              );
              %plot data
              sat.A.plot('columns',obj.par.calpar_csr.data_col_idx)
              sat.B.plot('columns',obj.par.calpar_csr.data_col_idx)
              %legend
              legend({'GRACE-A','GRACE-B'});
              obj.enforce_plot_par('autoscale',true,'automean',true,'outlier',3)
              grid on
              ylabel(sat.A.y_units{1}) %could also be sat.B.y_units{1}
              title([fields{j},' ',level])
              saveas(gcf,filename)
            end
          end
        case 'stats'
          s=obj.data_op('load-stats',datatype,level,varargin{:});
          fields=fieldnames(obj.par.calpar_csr.fields);
          %plot number of data
          obj.par.plot.prefix='length.bars';
          for i=1:numel(fields)
            filename=obj.id('plot',datatype,level,fields{i},'');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              b=bar(...
                mean([...
                  datenum(s.(fields{i}).A.length.t),...
                  datenum(s.(fields{i}).B.length.t)...
                ],2),[...
                  s.(fields{i}).A.length.scale(0.5).y(:,1),...
                  s.(fields{i}).B.length.scale(0.5).y(:,1)...
                ]...
              );
              b(1).EdgeColor = 'red';
              b(1).FaceColor = 'red';
              b(2).EdgeColor = 'black';
              b(2).FaceColor = 'black';
              %plot tweaking
              obj.enforce_plot_par
              legend('GRACE-A','GRACE-B')
              datetick('x','yyyy-mm')
              xlabel('time')
              grid on
              title(['Monthly count of ',fields{i},' ',level])
              saveas(gcf,filename)
            end
          end
          %plot more stats
          for mode=obj.par.modes.stats;
            %plot the requested statistics of the data
            obj.par.plot.prefix=mode{1};
            for i=1:numel(fields)
              filename=obj.id('plot',datatype,level,fields{i},'');
              if isempty(dir(filename))
                switch mode{1}
                case 'length'
                  continue
                  mask=true(s.(fields{i}).A.(mode{1}).length,1); %#ok<UNRCH>
                  scl=0.5;
                otherwise
                  mask=s.(fields{i}).A.length.y(:,1)>10;
                  scl=1;
                end
                sat=struct(...
                  'A',s.(fields{i}).A.(mode{1}).mask_and(mask).mask_update,...
                  'B',s.(fields{i}).B.(mode{1}).mask_and(mask).mask_update...
                );
                figure('visible',obj.par.plot.visible);
                sat.A.scale(scl).plot('columns',obj.par.calpar_csr.data_col_idx);
                sat.B.scale(scl).plot('columns',obj.par.calpar_csr.data_col_idx);
                %plot tweaking
                obj.enforce_plot_par('autoscale',true,'automean',true,'outlier',1)
                legend('GRACE-A','GRACE-B')
                ylabel(['[',sat.A.y_units{obj.par.calpar_csr.data_col_idx},']'])
                grid on
                title(['monthly ',mode{1},' of ',fields{i},' ',level])
                saveas(gcf,filename)
              end
            end
          end
        case 'corr'
          s=obj.data_op('load-stats2',datatype,level,varargin{:});
          fields=fieldnames(obj.par.calpar_csr.fields);
          % loop over all required modes
          for mode=obj.par.modes.stats2;
            %plot the requested statistics of the data
            obj.par.plot.prefix=mode{1};
            % loop over CSR cal pars
            for i=1:numel(fields)
              filename=obj.id('plot',datatype,level,fields{i},'');
              if isempty(dir(filename))
                switch mode{1}
                case 'length'
                  continue
                  mask=true(s.(fields{i}).(mode{1}).length,1); %#ok<UNRCH>
                  ylimits=[-Inf,Inf];
                  scl=0.5;
                otherwise
                  mask=s.(fields{i}).length.y(:,1)>=10;
                  ylimits=[-1,1];
                  scl=1;
                end    
                figure('visible',obj.par.plot.visible);
                s.(fields{i}).(mode{1}).mask_and(mask).mask_update.scale(scl).plot(...
                  'columns',obj.par.calpar_csr.data_col_idx...
                );
                obj.enforce_plot_par('ylimits',ylimits)
                datetick('x','yyyy-mm')
                xlabel('time')
                ylabel('[ ]')
                grid on
                title(['monthly ',mode{1},' of ',fields{i},' ',level])
                saveas(gcf,filename)
              end
            end
          end
        case 'overview'
          fields={...
            'AC0X',...
            'AC0Y1',...
            'AC0Z',...
            'AC0XD',...
            'AC0YD1',...
            'AC0ZD',...
            'AC0XQ',...
            'AC0YQ1',...
            'AC0ZQ'...
          };
          %single-sat stats
          s=obj.data_op('load-stats',datatype,level,varargin{:});
          stats_modes=obj.par.modes.stats;
          for j=1:numel(stats_modes);
            %skip length, not so interesting in this mode
            if strcmp(stats_modes{j},'length')
              continue
            end
            obj.par.plot.prefix=stats_modes{j};
            filename=obj.id('plot',datatype,level,'','');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %gathering data
              y=nan(numel(fields),2);
              for i=1:numel(fields)
                if ~isfield(s,fields{i})
                  continue
                end
                y(i,1)=s.(fields{i}).A.(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.calpar_csr.data_col_idx);
                y(i,2)=s.(fields{i}).B.(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.calpar_csr.data_col_idx);
              end
              pos_idx=(y>0);
              %plot negative parameters downwards
              y_now=abs(y); y_now(pos_idx)=nan;
              b{1}=bar(log10(y_now)); hold on
              %plot positive parameters upwards
              y_now=abs(y); y_now(~pos_idx)=nan;
              b{2}=bar(-log10(y_now)); hold on
              %annotating
              obj.enforce_plot_par
              for k=1:2
                b{k}(1).EdgeColor = 'red';
                b{k}(1).FaceColor = 'red';
                b{k}(2).EdgeColor = 'black';
                b{k}(2).FaceColor = 'black';
              end
              if all(pos_idx(~isnan(y)))
                axis([0 numel(fields)+1 0 10])
              else
                axis([0 numel(fields)+1 -10 10])
              end
              set(gca,'XTick',1:numel(fields))
              set(gca,'XTickLabel',fields)
              set(gca,'XTickLabelRotation',45)
              yticks=get(gca,'YTickLabel');
              for k=1:numel(yticks)
                ytickval=str2double(yticks{k});
                if ytickval==0
                  yticks{k}='1';
                elseif ytickval>0
                  yticks{k}=['y.10^{-',yticks{k},'}'];
                elseif ytickval<0
                  yticks{k}=['-y.10^{',yticks{k},'}'];
                end
              end
              set(gca,'YTickLabel',yticks);
              legend('GRACE-A','GRACE-B')
              grid on
              title([...
                stats_modes{j},' of ',level,' for ',...
                datestr(obj.par.plot.xlimits(1)),' to ',...
                datestr(obj.par.plot.xlimits(2))...
              ]);
              saveas(gcf,filename)
            end
          end
          %all-sat stats
          s=obj.data_op('load-stats2',datatype,level,varargin{:});
          stats_modes=obj.par.modes.stats2;
          for j=1:numel(stats_modes);
            %skip length, not so interesting in this mode
            if strcmp(stats_modes{j},'length')
              continue
            end
            obj.par.plot.prefix=stats_modes{j};
            filename=obj.id('plot',datatype,level,'','');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %gathering data
              y=nan(numel(fields),1);
              for i=1:numel(fields)
                if ~isfield(s,fields{i})
                  continue
                end
                y(i,1)=s.(fields{i}).(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.(datatype).data_col_idx);
              end
              bar(y,'k')
              %annotating
              obj.enforce_plot_par
              set(gca,'XTick',1:numel(fields))
              set(gca,'XTickLabel',fields)
              set(gca,'XTickLabelRotation',45)
              axis([0 numel(fields)+1 -1 1])
              grid on
              title([...
                stats_modes{j},' of ',level,' for ',...
                datestr(obj.par.plot.xlimits(1)),' to ',...
                datestr(obj.par.plot.xlimits(2))...
              ]);
              saveas(gcf,filename)
            end
          end  
%           case 'calpar_csr-tables'
%             fields=fieldnames(obj.par.calpar_csr.fields);
%             tab=12;
%             %single-sat stats
%             s=obj.data_op('load-stats',datatype,level,varargin{:});
%             stats_modes=obj.par.modes.stats;
%             for s=1:numel(grace.sats)
%               disp(['GRACE-',upper(grace.sats{s})])
%               for i=0:numel(fields)
%                 out=cell(1,numel(stats_modes)+1);
%                 if i==0
%                   calpar_clean='cal par';
%                 else
%                   calpar_clean=str.clean(fields{i},'_');
%                 end
%                 out{1}=str.tabbed(calpar_clean,tab);
%                 for j=1:numel(stats_modes);
%                   if i==0
%                     out{j+1}=str.tabbed(stats_modes{j},tab,true);
%                   else
%                     d=s{i}.(grace.sats{s}).(stats_modes{j}).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
%                     out{j+1}=str.tabbed(num2str(d(obj.par.calpar_csr.data_col_idx),'% .3g'),tab,true);
%                   end
%                 end
%                 disp(strjoin(out,' '))
%               end
%             end
%             %all-sat stats
%             s=obj.data_op('load-stats2',datatype,level,varargin{:});
%             stats_modes=obj.par.modes.stats2;
%             for i=0:numel(fields)
%               out=cell(1,numel(stats_modes)+1);
%               if i==0
%                 calpar_clean='cal par';
%               else
%                 calpar_clean=str.clean(fields{i},'_');
%               end
%               out{1}=str.tabbed(calpar_clean,tab);
%               for j=1:numel(stats_modes);
%                 if i==0
%                   out{j+1}=str.tabbed(stats_modes{j},tab,true);
%                 else
%                   d=s{i}.(stats_modes{j})(:,1).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
%                   out{j+1}=str.tabbed(num2str(mean(d(obj.par.calpar_csr.data_col_idx)),'% .3f'),tab,true);
%                 end
%               end
%               disp(strjoin(out,' '))
%             end
        otherwise
          error([mfilename,': cannot handle plotname ''',plotname,''' for datatype ''',datatype,'''.'])
        end
      case 'acc'
        switch lower(plotname)
        case {'x','y','z'}
          for s=1:numel(grace.sats)
            data_col_idx=find(~cellfun(@isempty,strfind(obj.par.acc.data_col_name,upper(plotname))));
            data_col_name=obj.par.acc.data_col_name{data_col_idx}; 
            filename=obj.id('plot',datatype,level,...
              [data_col_name,'.',datestr(obj.start,'yyyymmddTHHMMSS'),'-',datestr(obj.stop,'yyyymmddTHHMMSS')],...
            grace.sats{s});
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              ts={...
                obj.data.(datatype).l1b.csr.(       grace.sats{s}),...
                obj.data.(datatype).calmod.(level).(grace.sats{s}),...
                obj.data.(datatype).mod.csr.(       grace.sats{s}),...
                obj.data.(datatype).mod.nrtdm.(     grace.sats{s})...
              };
              idx=0;
              legend_str=cell(0);
              idx=idx+1;
              if isa(ts{1},'simpletimeseries')
                h=ts{1}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['L1B         ',num2str(h.y_mean{1})];
                y_units=ts{1}.y_units{data_col_idx};
              end
              idx=idx+1;
              if isa(ts{3},'simpletimeseries')
                calib=ts{1}+ts{2};
                h=calib.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['Calibrated  ',num2str(h.y_mean{1})];
              end
              idx=idx+1;
              if isa(ts{2},'simpletimeseries')
                h=ts{2}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['Calib. Mod. ',num2str(h.y_mean{1})];
                y_units=ts{2}.y_units{data_col_idx};
              end
              idx=idx+1;
              if isa(ts{3},'simpletimeseries')
                h=ts{3}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['CSR Mod.    ',num2str(h.y_mean{1})];
                y_units=ts{3}.y_units{data_col_idx};
              end
              idx=idx+1;
              if isa(ts{4},'simpletimeseries')
                h=ts{4}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['NRLMSISE-00 ',num2str(h.y_mean{1})];
                y_units=ts{4}.y_units{data_col_idx};
              end
              if all(cellfun(@isempty,legend_str))
                disp(['No data to plot ',filename])
                return
              end
              %plot tweaking
              obj.enforce_plot_par('autoscale',true,'automean',true)
              %add legend
              legend_idx=cellfun(@isempty,legend_str);
              lh=legend(legend_str{~legend_idx});
              set(lh,'fontname','courier')
              ylabel(['[',y_units,']'])
              grid on
              title([level,' ACC ',lower(plotname),'-componend GRACE-',grace.sats{s}])
              saveas(gcf,filename)
            end            
          end
        case 'stats'
          stats=obj.data_op('load-stats',datatype,level,varargin{:});
          models=fieldnames(stats.(level));
          stat_modes=obj.par.modes.stats;
          for s=1:numel(grace.sats)
            for sm=1:numel(stat_modes)
              for i=1:numel(obj.par.(datatype).data_col_idx)
                data_col_idx=obj.par.(datatype).data_col_idx(i);
                data_col_name=obj.par.acc.data_col_names{data_col_idx};
                filename=obj.id('plot',datatype,level,[stat_modes{sm},'-',data_col_name],grace.sats{s});
                if isempty(dir(filename))
                  figure('visible',obj.par.plot.visible);
                  legend_str=cell(1,numel(models));
                  for m=1:numel(models)
                    %get current stats timeseries
                    ts=stats.(level).(models{m}).(grace.sats{s}).(stat_modes{sm});
                    %only show stats when there is enough data
                    switch stat_modes{1}
                    case 'length'
                      mask=true(ts.length,1);
                      scl=0.5;
                    otherwise
                      mask=stats.(level).(models{m}).(grace.sats{s}).length.y(:,1)>10;
                      scl=1;
                    end
                    %plot it
                    h=ts.mask_and(mask).mask_update.scale(scl).plot('columns',data_col_idx);
                    legend_str{m}=[models{m},' \mu=',h.y_mean{data_col_idx}];
                  end
                  %plot tweaking
                  obj.enforce_plot_par('autoscale',true,'automean',true)
                  legend(legend_str)
                  ylabel(['[',ts.y_units{data_col_idx},']'])
                  grid on
                  title(['daily ',stat_modes{1},' of ',level,'-calib. acc. residuals along ',data_col_name,...
                    'for GRACE-',grace.sats{s}])
                  saveas(gcf,filename)
                end
              end
            end
          end
        otherwise
          error([mfilename,': cannot handle plotname ''',plotname,''' for datatype ''',datatype,'''.'])
        end
      otherwise
        error([mfilename,': unknown data of type ''',datatype,'''.'])
      end
    end
  end
end



function out=collapse(in)
  if ~isscalar(in) && ~iscell(in)
    %vectors are always lines (easier to handle strings)
    out=transpose(in(:));
  else
    out=in;
  end
end
