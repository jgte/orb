classdef dataproduct
  %static
  properties(Constant)
    default_list=struct(...
      'metadata_dir',fullfile(dataproduct.scriptdir,'metadata'),...
          'plot_dir',fullfile(dataproduct.scriptdir,    'plot'),...
          'data_dir',fullfile(dataproduct.scriptdir,    'data'),...
      'metadata',struct(...
        'plot_columns',{{-1}},...
        'plot_column_names',{{}},...
        'plot_columns_together',true,...
        'plot_together','',...
        'plot_file_prefix','',...
        'plot_file_suffix','',...
        'plot_file_full_path',true,...
        'plot_legend',{{}},...
        'plot_legend_suppress',{{}},...
        'plot_ylabel','',...
        'plot_label_prefix_sat','',...
        'plot_label_prefix_field','',...
        'plot_label_prefix_level','',...
        'plot_label_prefix_type','',...
        'plot_xdate',true,...
        'plot_xdateformat','HH:MM',...
        'plot_xlimits',[-inf,inf],...
        'plot_ylimits',[-inf,inf],...
        'plot_size',200+[0,0,21,9]*50,...
        'plot_units','points',...
        'plot_visible','on',...
        'plot_fontsize_axis', 24,...
        'plot_fontsize_title',32,...
        'plot_fontsize_label',28,...
        'plot_title','',...
        'plot_title_suppress',{{}},...
        'plot_title_suffix','',...
        'plot_grid',true,...
        'plot_line_width',2,...
        'plot_autoscale',false,... %y-scale is derived from the data (in dataproduct.enforce_plot)
        'plot_automean',false,...  %middle-point of y axis is derived from the data (in dataproduct.enforce_plot)
        'plot_zeromean',false,...  %mean of data is removed before plotting (in simpledata.plot)
        'plot_colormap','',...
        'plot_outlier',0,...
        'plot_method','timeseries'...
      )...
    );
    parameter_list=struct(...
      'metadata_dir',   struct('default',dataproduct.default_list.metadata_dir,   'validation',@(i) ischar(i)),...
      'plot_dir',       struct('default',dataproduct.default_list.plot_dir,       'validation',@(i) ischar(i)),...
      'data_dir',       struct('default',dataproduct.default_list.data_dir,       'validation',@(i) ischar(i)),...
      'metadata',       struct('default',dataproduct.default_list.metadata,       'validation',@(i) isstruct(i))...
    );
    %this is the maximum depth of the structure inside on dataname, i.e. obj.data.(dataname.name)
    %this parameter is arbitrary but keep it as low as possible to prevent time-consuming searches.
    maxdepth=4;
  end
  %read-write
  properties
    dataname
    metadata
  end
  %read only
  properties(SetAccess=private)
    metadata_dir
    data_dir
    plot_dir
  end
  %internal
  properties(SetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
    name
    codename
  end
  methods(Static)
    function out=parameters
      out=fieldnames(dataproduct.parameter_list);
    end
    function out=scriptdir
      out=fileparts(which(mfilename));
      if isempty(out)
        out='.';
      end
    end
    function out=parse_commands(in)
      %sanity
      if ~ischar(in)
        error([mfilename,': can only handle strings, not ',class(in),'.'])
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
  end
  methods
    %% constructor
    function obj=dataproduct(in,varargin)
      % trivial call (all other classes are handled in datanames)
      if isa(in,'dataproduct')
        obj=in;
        return
      end
      %need to read YAML
      addpath(genpath(fullfile(dataproduct.scriptdir,'yamlmatlab')));
      %parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      %declare parameters
      for j=1:numel(dataproduct.parameters)
        %shorter names
        pn=dataproduct.parameters{j};
        %declare parameters
        p.addParameter(pn,dataproduct.parameter_list.(pn).default,dataproduct.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(varargin{:});
      % call superclass constructor
      obj.dataname=datanames(in);
      % save parameters
      for i=1:numel(dataproduct.parameters)
        %shorter names
        pn=dataproduct.parameters{i};
        if ~isscalar(p.Results.(pn))
          %vectors are always lines (easier to handle strings)
          obj.(pn)=transpose(p.Results.(pn)(:));
        else
          obj.(pn)=p.Results.(pn);
        end
      end
      % load metadata
      obj=obj.metadata_load;
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
          'use_storage_period',true,...
          'dir',obj.data_dir,...
          'ext','mat'...
        };
      case 'plot'
        out={...
          'use_storage_period',false,...
          'dir',obj.plot_dir,...
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
      p.addParameter('use_storage_period',true,    @(i) islogical(i))
      p.addParameter('discover',          false,   @(i) islogical(i))
      % parse it
      p.parse(varargin{:});
      % retrive arguments for this file type
      typeargs=obj.file_args(mode);
      % build filename
      filename=obj.dataname.file(...
        typeargs{:},varargin{:},'keeptsplaceholder',p.Results.use_storage_period...
      );
      %resolve time stamp
      if p.Results.use_storage_period
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
          if iscell(obj.metadata.plot_xlimits(1));xl=obj.metadata.plot_xlimits{1};
          else                                    xl=obj.metadata.plot_xlimits(1);
          end
          if isfinite(xl);startlist=xl;else startlist=[];end
          if iscell(obj.metadata.plot_xlimits(2));xl=obj.metadata.plot_xlimits{2};
          else                                    xl=obj.metadata.plot_xlimits(2);
          end
          if isfinite(xl);stoplist=xl;else stoplist=[];end
          return
        case {'none','direct'}
          out={};
          startlist=datetime(0,'ConvertFrom','datenum');
          stoplist =datetime(0,'ConvertFrom','datenum');
          return
        otherwise
          error([mfilename,': cannot handle metadata key ''storage_period'' with value ''',obj.mdget('storage_period'),'''.'])
        end
        %branch on requested behaviour
        if p.Results.discover
          %discover existing files
          out=file.wildcard(strrep(filename,'<TIMESTAMP>','*'),'disp',false);
          %enforce the right kind of empty
          if isempty(out);out={};end
        else
          %build list of files
          out=cell(size(startlist));
          for i=1:numel(startlist)
            out{i}=strrep(filename,'<TIMESTAMP>',datestr(startlist(i),timestamp_fmt));
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
        out=cellfun(@(i)(~isempty(dir(i))),file_list);
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
    function out=md_file(obj)
      assert(logical(exist(obj.metadata_dir,'dir')),[mfilename,': ',...
        'cannot find metadata dir ''',obj.metadata_dir,'''.'])
      % build filename and add path
      out=fullfile(obj.metadata_dir,[obj.dataname.name,'.metadata']);
    end
    function out=ismd_file(obj)
      out=~isempty(dir(obj.md_file));
    end
    function md_file_check(obj)
      assert(obj.ismd_file,[mfilename,': ',...
        'could not find the metadata for product ',obj.name,' (expecting ',obj.md_file,').'])
    end    
    function obj=metadata_load(obj)
      obj.md_file_check
      obj=obj.metadata_merge(ReadYaml(obj.md_file));
      %load sub-metadata files
      c=0;
      while true
        c=c+1;
        fieldname=['submetadata',num2str(c)];
        if obj.ismd_field(fieldname)
          submetadatafile=fullfile(obj.metadata_dir,obj.mdget(fieldname));
          %sanity
          assert(~isempty(dir(submetadatafile)),['Cannot find sub-metadatafile ',submetadatafile,'.'])
          obj=obj.metadata_merge(ReadYaml(submetadatafile));
        else
          break
        end
      end
    end
    function obj=metadata_merge(obj,new_metadata)
      % entries in new_metadata replace those already in obj.metadata
      f=fieldnames(new_metadata);
      for i=1:numel(f)
        obj.metadata.(f{i})=new_metadata.(f{i});
      end
      %translate sources as strings to datanames
      for i=1:obj.nr_sources
        obj.metadata.sources{i}=datanames(obj.metadata.sources{i});
      end
    end
    %% metadata fields
    function out=ismd_field(obj,metadatafieldname)
      out=isfield(obj.metadata,metadatafieldname);
    end
    function md_field_check(obj,metadatafieldname)
      if ~obj.ismd_field(metadatafieldname)
        error([mfilename,': cannot find field ',metadatafieldname,' in the metadata of ',obj.name,'.'])
      end
    end
    %% metadata parsing
    function out=mdget(obj,metadatafieldname,varargin)
      p=inputParser;
      p.addRequired( 'metadatafieldname',            @(i) ischar(i));
      p.addParameter('return_empty_if_missing',false,@(i) islogical(i));
      p.addParameter('always_cell_array',      false,@(i) islogical(i));
      p.addParameter('default',                '',   @true); %anything goes
      % parse it
      p.parse(metadatafieldname,varargin{:});
      % check for existence (unless return_empty_if_missing)
      if ~p.Results.return_empty_if_missing && ~isempty(p.Results.default)
        obj.md_field_check(metadatafieldname)
      end
      % check if default is to be returned
      return_default=~obj.ismd_field(metadatafieldname);
      if return_default
        out=p.Results.default;
      else
        in=obj.metadata.(metadatafieldname);
        switch class(in)
        case {'double','logical','struct'}
          out=in;
        case 'cell'
          if iscellstr(in)
            out=cellfun(@dataproduct.parse_commands,in,'UniformOutput',false);
          else
            out=in;
          end
        case 'char'
          out=dataproduct.parse_commands(in);
        otherwise
          error([mfilename,': cannot handle class ',class(in),'. Implementation needed!'])
        end
      end
      if ~iscell(out) && p.Results.always_cell_array
        out={out};
      end
    end
    %% data sources (usually spits out type 'datanames')
    function out=nr_sources(obj)
      if isfield(obj.metadata,'sources')
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
      out=cell(1,obj.nr_sources);
      for i=1:obj.nr_sources
        out{i}=obj.sources(i);
      end
    end
    %% plot customization
    function p=plot_args(obj,p,varargin)
      %get plot parameters already defined in the metadata
      mdf=fieldnames(obj.metadata);
      %loop over all of them
      for i=1:numel(mdf)
        %skip if this is not a plot parameter
        if isempty(strfind(mdf{i},'plot_'))
          continue
        end
        %declare this plot parameter as optional argument to this function
        p.addParameter(...
          mdf{i},...
          num.cell(obj.metadata.(mdf{i})),...
          @(i) iscell(i) || ischar(i) || isnumeric(i) || islogical(i)... %pretty generic validation needed to let anything in
        )
      end
      % parse arguments
      p.parse(varargin{:});
    end
    function out=plot_args_list(obj)
      %get all parameters already defined in the metadata
      mdf=fieldnames(obj.metadata);
      %filter out those that are not plot_*
      mdf=mdf(cellfun(@(i) ~isempty(strfind(i,'plot_')),mdf));
      %make room for outputs
      out=cell(size(mdf));
      %loop over all of them
      for i=1:numel(mdf)
        %append to varargout
        out{2*i-1}=mdf{i};
        out{2*i}=obj.mdget(mdf{i});
      end
    end
    function enforce_plot(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('fig_handle',  gcf,  @(i) ishandle(i));
      p.addParameter('axis_handle', gca,  @(i) ishandle(i));
      %get plot parameters already defined in the metadata and parse them
      p=obj.plot_args(p,varargin{:});
      % sanity
      if any(isfinite(p.Results.plot_ylimits)) && ( p.Results.plot_autoscale ||  p.Results.plot_automean )
        error([mfilename,': option ''ylimits'' and ''autoscale'' or ''automean'' do not work concurrently.'])
      end
      % enforce fontsize and paper size
      set(    p.Results.axis_handle,          'FontSize',p.Results.plot_fontsize_axis)
      set(get(p.Results.axis_handle,'Title' ),'FontSize',p.Results.plot_fontsize_title);
      set(get(p.Results.axis_handle,'XLabel'),'FontSize',p.Results.plot_fontsize_label);
      set(get(p.Results.axis_handle,'YLabel'),'FontSize',p.Results.plot_fontsize_label);
      set(    p.Results.fig_handle, 'Position',          p.Results.plot_size,...
                                    'PaperUnits',        p.Results.plot_units,...
                                    'PaperPosition',     p.Results.plot_size);
      % enforce line properties
      line_handles=plotting.line_handles(p.Results.axis_handle);
      for i=1:numel(line_handles)		
        set(line_handles(i),'LineWidth',p.Results.plot_line_width)		
      end		
      % enforce x-limits
      v=axis(p.Results.axis_handle);
      if p.Results.plot_xdate
        for i=1:2
          if iscell(p.Results.plot_xlimits(i))
            xl=p.Results.plot_xlimits{i};
          else
            xl=p.Results.plot_xlimits(i);
          end
          if isfinite(   xl)
            v(i)=datenum(xl);
          end
        end
        if ~strcmp(datestr(v(1),'yyyymmdd'),datestr(v(2),'yyyymmdd')) && ...
            (~strcmp(datestr(v(2),'HHMMSS'),'000000') || v(2)-v(1)>1)
          xlabel([datestr(v(1),'yyyy-mm-dd'),' to ',datestr(v(2),'yyyy-mm-dd')])
        else
          xlabel(datestr(v(1)))
        end
        datetick('x',p.Results.plot_xdateformat)
      end
      % enforce reasonable y-limits
      for i=1:2
        if isfinite(p.Results.plot_ylimits(i))
          v(i+2)=   p.Results.plot_ylimits(i);
        end
      end
      % enforce auto-scale and/or auto-mean
      if p.Results.plot_automean || p.Results.plot_autoscale || p.Results.plot_outlier>0
        %gather plotted data
        dat=cell(size(line_handles));
        for i=1:numel(line_handles)
          dat{i}=get(line_handles(i),'ydata');
        end
        dat=[dat{:}];
        dat=dat(~isnan(dat(:)));
        %remove outliers before computing axis limits (if requested)
        for c=1:p.Results.plot_outlier
          dat=simpledata.rm_outliers(dat);
        end
        dat=dat(~isnan(dat(:)));
        if ~isempty(dat) && diff(minmax(dat))~=0
          %enfore data-driven mean and/or scale
          if p.Results.plot_automean && p.Results.plot_autoscale
            v(3:4)=mean(dat)+4*std(dat)*[-1,1];
          elseif p.Results.plot_automean
            v(3:4)=mean(dat)+0.5*diff(v(3:4))*[-1,1];
          elseif p.Results.plot_autoscale
            v(3:4)=mean(v(3:4))+4*std(dat)*[-1,1];
          end
          %fix non-negative data sets
          if all(dat>=0) && v(3)<0
            v(3)=0;
          end
        end
      end
      axis(p.Results.axis_handle,v);
      %enforce colormap
      if ~isempty(p.Results.plot_colormap)
        colormap(p.Results.plot_colormap)
      end
    end
    %% level-wrapping handlers
    function out=islevel_wrapped(obj,level)
      out=obj.ismd_field(dataproduct.level_vals_str(level))&&obj.ismd_field(dataproduct.level_name_str(level));
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
    function out=field_path_expand_sources(obj)
      out=obj.field_path_expand(obj.sources);
    end
  end
end