classdef dataproduct
  %static
  properties(Constant)
    default_metadata=struct(...
      'plot_product',false...
    );
%       'plot_columns',{{-1}},...
%       'plot_column_names',{{}},...
%       'plot_column_together',true,...
%       'plot_together','',...
%       'plot_file_prefix','',...
%       'plot_file_suffix','',...
%       'plot_legend',{{}},...
%       'plot_legend_suppress',{{}},...
%       'plot_legend_location','best',...
%       'plot_ylabel','',...
%       'plot_xdate',true,...
%       'plot_xdateformat','HH:MM',...
%       'plot_xlimits',[-inf,inf],...
%       'plot_ylimits',[-inf,inf],...
%       'plot_size',200+[0,0,21,9]*50,...
%       'plot_units','points',...
%       'plot_visible',true,...
%       'plot_fontsize_axis', 24,...
%       'plot_fontsize_title',32,...
%       'plot_fontsize_label',28,...
%       'plot_title','',...
%       'plot_title_suppress',{{}},...
%       'plot_title_suffix','',...
%       'plot_title_prefix','',...
%       'plot_grid',true,...
%       'plot_line_width',2,...
%       'plot_autoscale',false,... %y-scale is derived from the data (in dataproduct.enforce_plot)
%       'plot_automean',false,...  %middle-point of y axis is derived from the data (in dataproduct.enforce_plot)
%       'plot_zeromean',false,...  %mean of data is removed before plotting (in simpledata.plot)
%       'plot_colormap','',...
%       'plot_outlier',0,...
%       'plot_method','timeseries',...
%       'plot_normalize',false...
%     );
  
  %TODO: plot_dir needs to be duplicated in the metadata (so that plotting routines can use it)
  
    parameter_list={...
      'metadata_dir',   fullfile(dataproduct.scriptdir,'metadata'), @(i) ischar(i);...
      'plot_dir',       fullfile(dataproduct.scriptdir,    'plot'), @(i) ischar(i);...
      'data_dir',       fullfile(dataproduct.scriptdir,    'data'), @(i) ischar(i);...
      'metadata',       dataproduct.default_metadata,               @(i) isstruct(i);...
      'debug_plot',     false,                                      @(i) islogical(i)...
    };
    %this is the maximum depth of the structure inside on dataname, i.e. obj.data.(dataname.name)
    %this parameter is arbitrary but keep it as low as possible to prevent time-consuming searches.
    maxdepth=4;
  end
  %read-write
  properties
    dataname
    metadata
    debug_plot
    metadata_dir
    data_dir
    plot_dir
  end
  %calculated only when asked for
  properties(Dependent)
    name
    codename
  end
  methods(Static)
    function out=parameters(i,method)
      persistent v parameter_names
      if isempty(v)
        v=varargs(dataproduct.parameter_list);
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
    end
    function out=scriptdir
      out=fileparts(which(mfilename));
      if isempty(out)
        out='.';
      end
    end
    function out=parse_commands(in)
      %sanity
      switch class(in)
      case 'char'
        %do nothing
      case 'cell'
        out=cellfun(@dataproduct.parse_commands,in,'UniformOutput',false);
        return
      case 'struct'
        out=structs.fun(@dataproduct.parse_commands,in);
        return
      otherwise
        %pass-through
        out=in;
        return
      end
      start_idx=strfind(in,'<');
       stop_idx=strfind(in,'>');
      %check if matlab code is not part of the field
      if isempty(start_idx) && isempty(stop_idx)
        out=strtrim(in);
        return
      end
      %parse matlab code if nr of < and > match
      if ~isempty(start_idx) && numel(start_idx) == numel(stop_idx)
        out=cell(2*numel(start_idx)+1,1);
        out{1}=in(1:start_idx-1);
        for i=1:numel(start_idx)
          com_str=in(start_idx(i)+1:stop_idx(i)-1);
          out{2*i-1}=eval(com_str);
          if i==numel(start_idx)
            out{2*i}=in(stop_idx(i)+1:end);
          else
            out{2*i}=in(stop_idx(i)+1:start_idx(i)-1);
          end
        end
        %cleanup empty entries
        out=out (~cellfun(@isempty,out));
        %do some tedious parsing
        if iscellstr(out)
          %assemble components
          out=strjoin(out,'');
        elseif numel(out)==1
          %unwrap scalar cell
          out=out{1};
        else
          %need to be handled externally
        end
      else
        error([mfilename,': could not parse string ''',in,''' because of unmatching ''<'' and ''>'' characters.'])
      end
    end
    function out=level_vals_str(level)
      out=['level',num2str(level),'_vals'];
    end
    function out=level_name_str(level)
      out=['level',num2str(level),'_name'];
    end
    function products_out=unwrap_product(products_in,level)
      if ~exist('level','var')
        products_out=products_in;
        for i=1:dataproduct.maxdepth
          products_out=dataproduct.unwrap_product(products_out,i);
        end
      else
        %init outputs and output counter
        products_out=cell(size(products_in));c=0;
        %loop over all input products
        for i=1:numel(products_in)
          %check if unwrapping of this part is needed
          if products_in{i}.islevel_wrapped(level)
            %retrieve values of the wrapped field
            level_vals=products_in{i}.level_vals(level);
            %retrieve name of the wrapped field
            level_name=products_in{i}.level_name(level);
            for p=1:numel(level_vals)
              %patch product
              product_patched=products_in{i};
              product_patched.metadata=rmfield(product_patched.metadata,...
                {dataproduct.level_vals_str(level),dataproduct.level_name_str(level)}...
              );
              product_patched.metadata.(level_name)=level_vals{p};
              product_patched.dataname=products_in{i}.dataname.append_field_leaf([level_name,'_',str.show(level_vals{p})]);
              % append to product list
              c=c+1; products_out{c}=product_patched;
            end
          else
            %propagate current wrap-less product
            c=c+1; products_out{c}=products_in{i};
          end
        end
      end
    end
    function out=array(in)
      assert(iscell(in),[mfilename,': cannot handle input ''in'' of class ',class(in),', expecting a cell array.'])
      if isempty(in)
        out={};
      else
        out=cellfun(@dataproduct,in,'UniformOutput',false);
      end
    end
    function out=mdfile_static(dn,metadata_dir)
      assert(exist(metadata_dir,'dir')~=0,[mfilename,': ',...
        'cannot find metadata dir ''',metadata_dir,'''.'])
      % build filename and add path
      out=fullfile(metadata_dir,[dn.name,'.metadata']);
    end
  end
  methods
    %% constructor
    function obj=dataproduct(in,varargin)
      % trivial call (all other classes are handled in datanames)
      if isa(in,'dataproduct')
        obj=in;
        return
      end
      % input parsing (relevant to the parameters defined in dataproduct.parameter_list)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('field_path',    '',@(i) ischar(i) || iscell(i));
      p.addParameter('metadata_from', '',@(i) isstruct(i));
      %create argument object, declare and parse parameters, save them to obj
      [~,p,obj]=varargs.wrap('sinks',{obj},'parser',p,'sources',{dataproduct.parameters([],'obj')},varargin{:});
      % call superclass constructor
      if ~isempty(p.Results.field_path)
        obj.dataname=datanames(in,p.Results.field_path);
      else
        obj.dataname=datanames(in);
      end
      % make sure metadata file exists
      if ~obj.ismdfile
        %if there's no metadatafile, let's assume there's field_path in 'in' and start digging
        in_split=strsplit(obj.name,datanames.separator);
        for i=numel(in_split):-1:2
          dn_now=datanames(strjoin(in_split(1:i-1),datanames.separator),in_split(i:end));
          if exist(dataproduct.mdfile_static(dn_now,obj.metadata_dir),'file')~=0
            %found it!
            obj.dataname=dn_now;
          end
        end
      end
      % load metadata
      obj=obj.mdload(p.Results.metadata_from);
      % input parsing (relevant to the parameters defined in obj.metadata)
      [~,~,obj.metadata]=varargs.wrap('sources',{obj.default_metadata,obj.metadata},'sinks',{obj.metadata},varargin{:});
    end
    function obj=rm_data(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('rm_args',[], @(i) ischar(i));
      p.addParameter('mode','data', @(i) ischar(i));
      % parse it
      p.parse(varargin{:});
      %get list of files
      files=obj.file(p.Results.mode,varargin{:},'discover',true);
      %loop over all files and delete them
      for i=1:numel(files)
        if ~isempty(dir(files{i}))
          disp(['deleting ',files{i},'...'])
          com=['rm -vr ',p.Results.rm_args,' ',files{i}];
          if system(com)==0
            disp('done!')
          else
            error([mfilename,': command "',com,'" failed.'])
          end
        else
          disp([files{i},' non-existing, skipping...'])
        end
      end
    end
    %% filename handlers
    function out=file_args(obj,mode)
      switch lower(mode)
      case 'data'
        out={...
          'use_storage_period',true,... % I think this can be removed
          'dir',obj.data_dir,...
          'sub_dirs','single',...
          'ext','mat'...
        };
      case 'plot'
        out={...
          'use_storage_period',false,... % I think this can be removed
          'dir',obj.plot_dir,...
          'sub_dirs','none',...
          'ext','png'...
        };
      otherwise
        error([mfilename,': cannot handle mode ''',mode,'''.'])
      end
    end
    function [out,startlist,stoplist]=file(obj,mode,varargin)
      %NOTICE: this procedure ALWAYS returns cell arrays!
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('start',     datetime('now'), @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',      datetime('now'), @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('start_timestamp_only', true, @(i) islogical(i))
      p.addParameter('use_storage_period',   true, @(i) islogical(i))
      p.addParameter('discover',            false, @(i) islogical(i))
      % parse it
      p.parse(varargin{:});
      % retrive arguments for this file type
      typeargs=obj.file_args(mode);
      % build filename
      filename=obj.dataname.file(...
        typeargs{:},varargin{:},'keeptsplaceholder',p.Results.use_storage_period...
      );
      %branch on requested behaviour
      if p.Results.discover
        %it doesn't make sense to ask for the start/stop list and request the list of existing files
        assert(nargout==1,'when ''discover'' is true, outputs ''startlist'' and ''stoplist'' are ilegal.')
        %discover existing files
        out=file.wildcard(strrep(filename,'<TIMESTAMP>','*'),'disp',false);
        %enforce the right kind of empty
        if isempty(out);out={};end
      %resolve time stamp
      elseif p.Results.use_storage_period
        switch lower(obj.mdget('storage_period'))
        case {'day','daily'}
          [startlist,stoplist]=time.day_list(p.Results.start,p.Results.stop);
          timestamp_fmt='yyyymmdd';
        case {'month','monthly'}
          [startlist,stoplist]=time.month_list(p.Results.start,p.Results.stop);
          timestamp_fmt='yyyymm';
        case {'year','yearly'}
          [startlist,stoplist]=time.year_list(p.Results.start,p.Results.stop);
          timestamp_fmt='yyyy';
        case {'infinite','global'}
          out={strrep(filename,'.<TIMESTAMP>','')};
%           if iscell(obj.metadata.plot_xlimits(1));xl=obj.metadata.plot_xlimits{1};
%           else                                    xl=obj.metadata.plot_xlimits(1);
%           end
%           if isfinite(xl);startlist=xl;else startlist=[];end
          startlist=p.Results.start;
%           if iscell(obj.metadata.plot_xlimits(2));xl=obj.metadata.plot_xlimits{2};
%           else                                    xl=obj.metadata.plot_xlimits(2);
%           end
%           if isfinite(xl);stoplist=xl;else stoplist=[];end
          stoplist=p.Results.stop;
          return
        case {'none'}
          out={};
          startlist=datetime(0,'ConvertFrom','datenum');
          stoplist =datetime(0,'ConvertFrom','datenum');
          return
        case {'direct'}
          startlist=p.Results.start;
           stoplist=p.Results.stop;
           timestamp_fmt='yyyymmdd';
        otherwise
          error([mfilename,': cannot handle metadata key ''storage_period'' with value ''',obj.mdget('storage_period'),'''.'])
        end
        %build list of files
        out=cell(size(startlist));
        for i=1:numel(startlist)
          if p.Results.start_timestamp_only
            out{i}=strrep(filename,'<TIMESTAMP>',...
              datestr(startlist(i),timestamp_fmt)...
            );
          else
            out{i}=strrep(filename,'<TIMESTAMP>',[...
              datestr(startlist(i),timestamp_fmt),'T',datestr(stoplist(i),timestamp_fmt)...
            ]);
          end
        end
      else
        out={filename};
      end
    end
    function out=isfile(obj,varargin)
      file_list=obj.file(varargin{:});
      if isempty(file_list)
        out=false;
      else
        out=cellfun(@(i)exist(i,'file')~=0,file_list);
      end
    end
    function filecheck(obj,varargin)
      if any(~obj.isfile,varargin{:})
        f=obj.file(varargin{:});
        for i=1:numel(f)
          if isempty(dir(f))
            msg=[': for product ',obj.name,', could not find the following file: ',f];
            if numel(f)>1
              error([msg,', ',str.th(i),' of ',num2str(numel(f)),').'])
            else
              error([msg,').'])
            end
          end
        end
      end
    end
    %% name methods (easy viewing)
    function out=get.name(obj)
      out=obj.dataname.name;
    end
    function out=get.codename(obj)
      out=obj.dataname.codename;
    end
    function out=str(obj)
      out=obj.dataname.str;
    end
    %% metadata file
    function out=mdfile(obj)
      assert(logical(exist(obj.metadata_dir,'dir')),[mfilename,': ',...
        'cannot find metadata dir ''',obj.metadata_dir,'''.'])
      % build filename and add path
      out=dataproduct.mdfile_static(obj.dataname,obj.metadata_dir);
    end
    function out=ismdfile(obj)
      out=exist(obj.mdfile,'file')~=0;
    end
    function mdfile_check(obj)
      assert(obj.ismdfile,[mfilename,': ',...
        'could not find the metadata for product ',obj.name,' (expecting ',obj.mdfile,').'])
    end
    function obj=mdload(obj,metadata)
      if ~exist('metadata','var') || isempty(metadata)
        %need to read YAML
        addpath(genpath(fullfile(dataproduct.scriptdir,'yamlmatlab')));
        %make sure the metadata files is there
        obj.mdfile_check
        %load metadata
        metadata=ReadYaml(obj.mdfile);
        %convert numeric cell arrays to normal arrays
        metadata=structs.fun(@cells.c2m,metadata);
      end
      %load and merge with default/existing metadata
      obj=obj.mdmerge(metadata);
      %load sub-metadata files
      c=0;
      while true
        c=c+1;
        fieldname=['submetadata',num2str(c)];
        if obj.ismdfield(fieldname)
          submetadatafile=fullfile(obj.metadata_dir,obj.mdget(fieldname));
          %sanity
          assert(~isempty(dir(submetadatafile)),['Cannot find sub-metadatafile ',submetadatafile,'.'])
          obj=obj.mdmerge(ReadYaml(submetadatafile));
        else
          break
        end
      end
      %parse commands (if any)
      obj.metadata=dataproduct.parse_commands(obj.metadata);
    end
    function obj=mdmerge(obj,new_metadata)
      obj.metadata=structs.copy(new_metadata,obj.metadata);
      %need to translate scalar sources
      if isfield(obj.metadata,'sources') && ischar(obj.metadata.sources)
        obj.metadata.sources={obj.metadata.sources};
      end
      %translate sources as strings to datanames (the 'sources' metadata field could have been updated)
      for i=1:obj.nr_sources
        obj.metadata.sources{i}=datanames(obj.metadata.sources{i});
      end
    end
    %% metadata fields
    function out=ismdfield(obj,metadatafieldname)
      out=isfield(obj.metadata,metadatafieldname);
    end
    function mdfield_check(obj,metadatafieldname)
      if ~obj.ismdfield(metadatafieldname)
        error([mfilename,': cannot find field ',metadatafieldname,' in the metadata of ',obj.name,'.'])
      end
    end
    %% metadata parsing
    function out=mdget(obj,metadatafieldname,varargin)
      p=inputParser;
      p.addRequired( 'metadatafieldname',            @(i) ischar(i) || iscell(i));
      p.addParameter('return_empty_if_missing',false,@(i) islogical(i));
      p.addParameter('always_cell_array',      false,@(i) islogical(i));
      %if both 'default' and 'default_varargin' are given, the former takes precedence (but never over the metadata)
      p.addParameter('default',                '',   @(i) true); %anything goes
      p.addParameter('default_varargin',       {},   @(i) iscell(i) || isstruct(i));
      % parse it
      p.parse(metadatafieldname,varargin{:});
      % check if there's a default value defined
      isdefault_given=~( any(strcmp('default',p.UsingDefaults)) && any(strcmp('default_varargin',p.UsingDefaults)) );
      % vector mode
      if iscell(metadatafieldname)
        for i=1:numel(metadatafieldname)
          out.(metadatafieldname)=obj.mdget(metadatafieldname{i},varargin{:});
        end
        return
      end
      % check for existence (unless return_empty_if_missing)
      if ~p.Results.return_empty_if_missing && ~isdefault_given
        obj.mdfield_check(metadatafieldname)
      end
      % check if default is to be returned
      if ~obj.ismdfield(metadatafieldname)
        %this method can be called with varargin, useful to set defaults from the command call
        if ~isfield('default_varargin',p.UsingDefaults)
          switch class(p.Results.default_varargin)
          case 'cell'
            out=cells.varargin('get',p.Results.default_varargin,metadatafieldname);
          case 'struct'
            out=structs.get_value(p.Results.default_varargin,{metadatafieldname});
          end
        end
        if ~isfield('default',p.UsingDefaults)
          out=p.Results.default;
        end
      else
        out=obj.metadata.(metadatafieldname);
      end
      if ~iscell(out) && p.Results.always_cell_array
        out={out};
      end
    end
    function obj=mdset(obj,varargin)
      assert(mod(numel(varargin),2)==0,'Number of optional inputs must be even, debug needed!')
      for i=1:numel(varargin)/2
        obj.metadata.(varargin{2*i-1})=varargin{2*i};
      end
    end
    %% data sources (usually spits out type 'datanames')
    function out=nr_sources(obj)
      if isfield(obj.metadata,'sources')
        assert(iscell(obj.metadata.sources),['Metadata field ''sources'' must be of a cell array, not of class ',...
          class(obj.metadata.sources),'.'])
        out=numel(obj.metadata.sources);
      else
        out=0;
      end
    end
    function out=sources(obj,idx)
      assert(isnumeric(idx)&&isfinite(idx)&&isscalar(idx)&&idx>0,...
        ['input ''idx'' must be a scalar positive integer, not ''',num2str(idx),'''.'])
      assert(idx<=obj.nr_sources,...
        ['cannot get the ',str.th(idx),' data source because there are only ',num2str(obj.nr_sources),'.'])
      sources=obj.mdget('sources','always_cell_array',true);
      out=sources{idx};
      assert(isa(out,'datanames'),['source(',num2str(idx),') must be of class dataname, not ',class(out),', debug needed.'])
    end
    function out=source_list(obj)
      out=arrayfun(@(i) obj.sources(i),1:obj.nr_sources,'UniformOutput',false);
%       out=cell(1,obj.nr_sources);
%       for i=1:obj.nr_sources
%         out{i}=obj.sources(i);
%       end
    end
    function out=source_list_str(obj)
      out=cellfun(@(i) i.str,obj.source_list,'UniformOutput',false);
    end
    %% plot customization
    function out=plot_args(obj)
      %init outputs
      out={};
      %get plot parameters already defined in the metadata
      mdf=fieldnames(obj.metadata);
      %loop over all of them
      for i=1:numel(mdf)
        %skip if this is not a plot parameter
        if isempty(strfind(mdf{i},'plot_'))
          continue
        end
        %add this plot parameter to output
        out{end+1}=mdf{i}; %#ok<AGROW>
        out{end+1}=obj.metadata.(mdf{i}); %#ok<AGROW>
      end
    end
    function out=enforce_plot(obj,varargin)
      %get plot parameters already defined in the metadata and call mother routine
      out=plotting.enforce(obj.plot_args{:},varargin{:});
      % enforce fontsize and paper size
%       set(    p.Results.axis_handle,          'FontSize',p.Results.plot_fontsize_axis)
%       set(get(p.Results.axis_handle,'Title' ),'FontSize',p.Results.plot_fontsize_title);
%       set(get(p.Results.axis_handle,'XLabel'),'FontSize',p.Results.plot_fontsize_label);
%       set(get(p.Results.axis_handle,'YLabel'),'FontSize',p.Results.plot_fontsize_label);
%       set(    p.Results.fig_handle, 'Position',          p.Results.plot_size,...
%                                     'PaperUnits',        p.Results.plot_units,...
%                                     'PaperPosition',     p.Results.plot_size);
%       % enforce line properties
%       line_handles=plotting.line_handles(p.Results.axis_handle);
%       for i=1:numel(line_handles)
%         set(line_handles(i),'LineWidth',p.Results.plot_line_width)
%       end
%       % start with current axis
%       v=axis(p.Results.axis_handle);
%       % enforce x-limits
%       xl=cell(1,2);
%       for i=1:2
%         if iscell(p.Results.plot_xlimits(i))
%           xl{i}=p.Results.plot_xlimits{i};
%         else
%           xl{i}=p.Results.plot_xlimits(i);
%         end
%       end
%       %check if dates are requested
%       if p.Results.plot_xdate && ~p.Results.plot_psd
%         for i=1:2
%           if isfinite(xl{i})
%             v(i)=datenum(xl{i});
%           end
%         end
%         if ~strcmp(datestr(v(1),'yyyymmdd'),datestr(v(2),'yyyymmdd')) && ...
%             (~strcmp(datestr(v(2),'HHMMSS'),'000000') || v(2)-v(1)>1)
%           xlabel([datestr(v(1),'yyyy-mm-dd'),' to ',datestr(v(2),'yyyy-mm-dd')])
%         else
%           xlabel(datestr(v(1)))
%         end
%         datetick('x',p.Results.plot_xdateformat)
%       else
%         for i=1:2
%           if isfinite(xl{i})
%             v(i)=datenum(xl{i});
%           end
%         end
%       end
%       % enforce reasonable y-limits
%       for i=1:2
%         if isfinite(p.Results.plot_ylimits(i))
%           v(i+2)=   p.Results.plot_ylimits(i);
%         end
%       end
%       axis(p.Results.axis_handle,v);
%       %adjust the location of the legend (even if the default 'best', it may happen that it is no longer in a good position)
%       set(legend,'location',p.Results.plot_legend_location)
%       %enforce colormap
%       if ~isempty(p.Results.plot_colormap)
%         colormap(p.Results.plot_colormap)
%       end
    end
    %% level-wrapping handlers
    function out=islevel_wrapped(obj,level)
      out=obj.ismdfield(dataproduct.level_vals_str(level))&&obj.ismdfield(dataproduct.level_name_str(level));
    end
    function out=is_wrapped(obj)
      out=obj.islevel_wrapped(1);
    end
    function out=level_vals(obj,level)
      out=obj.mdget(dataproduct.level_vals_str(level));
      %need level_vals to be a cell array
      assert(iscell(out),['BUG TRAP: need metadata entry ''',dataproduct.level_vals_str(level),...
        ''' to be a cell array, not a ',class(out),'.'])
    end
    function out=level_name(obj,level)
      out=obj.mdget(dataproduct.level_name_str(level));
    end
    %% field_path expanding
    function out=field_path_expand(obj,dn_list)
      assert(iscell(dn_list),...
        ['input ''dn_list'' must be a cell array, not a ',class(dn_list),'.'])
      assert(all(cellfun(@(i) isa(i,'datanames'),dn_list)),...
        ['input ''dn_list'' must be a cell array of datanames, not ',...
        strjoin(cellfun(@class,dn_list,'UniformOutput',false),', '),'.'])
      %make room for outputs
      out=cell(size(dn_list));
      %patch expand this product
      for i=1:numel(dn_list)
        %propagate
        out{i}=obj;
        %edit dataname to point to leaf
        out{i}.dataname=out{i}.dataname.set_field_path(dn_list{i}.field_path);
        %patch sources
        for j=1:obj.nr_sources
          out{i}.metadata.sources{j}=out{i}.metadata.sources{j}.set_field_path(dn_list{i}.field_path);
        end
      end
    end
  end
end