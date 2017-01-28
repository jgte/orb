classdef dataproduct
  %static
  properties(Constant)
    default_list=struct(...
      'metadata_dir',fullfile(dataproduct.scriptdir,'metadata'),...
          'plot_dir',fullfile(dataproduct.scriptdir,    'plot'),...
          'data_dir',fullfile(dataproduct.scriptdir,    'data'),...
      'metadata',struct(...
        'plot_columns',{1},...
        'plot_prefix','',...
        'plot_suffix','',...
        'plot_xdate',true,...
        'plot_xlimits',[-inf,inf],...
        'plot_ylimits',[-inf,inf],...
        'plot_size',200+[0,0,21,9]*50,...
        'plot_units','points',...
        'plot_visible','on',...
        'plot_fontsize_axis', 24,...
        'plot_fontsize_title',32,...
        'plot_fontsize_label',28,...
        'plot_line_width',2,...
        'plot_autoscale',false,... %y-scale is derived from the data (in dataproduct.enforce_plot)
        'plot_automean',false,...  %middle-point of y axis is derived from the data (in dataproduct.enforce_plot)
        'plot_zeromean',false,...  %mean of data is removed before plotting (in simpledata.plot)
        'plot_outlier',0 ...
      )...
    );
    parameter_list=struct(...
      'metadata_dir',   struct('default',dataproduct.default_list.metadata_dir,   'validation',@(i) ischar(i)),...
      'plot_dir',       struct('default',dataproduct.default_list.plot_dir,       'validation',@(i) ischar(i)),...
      'data_dir',       struct('default',dataproduct.default_list.data_dir,       'validation',@(i) ischar(i)),...
      'metadata',       struct('default',dataproduct.default_list.metadata,       'validation',@(i) isstruct(i))...
    );
  
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
        out=in;
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
  end
  methods
    %% constructor
    function obj=dataproduct(in,varargin)
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
    %% filename handlers
    function out=file_args(obj,mode)
      switch lower(mode)
      case 'data'
        out={...
          'dir',obj.data_dir,...
          'ext','mat'...
        };
      case 'plot'
        out={...
          'dir',obj.plot_dir,...
          'ext','png'...
        };
      otherwise
        error([mfilename,': cannot handle mode ''',mode,'''.'])
      end
    end
    function out=file(obj,mode,varargin)
      %NOTICE: this procedure ALWAYS returns cell arrays!
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('start',     datetime('now'), @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',      datetime('now'), @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('use_storage_period',true,    @(i) islogical(i))
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
          timestamp_list=time.day_list(p.Results.start,p.Results.stop);
          timestamp_fmt='yyyymmdd';
        case {'month','monthly'}
          timestamp_list=time.month_list(p.Results.start,p.Results.stop);
          timestamp_fmt='yyyymm';
        case {'year','yearly'}
          timestamp_list=time.year_list(p.Results.start,p.Results.stop);
          timestamp_fmt='yyyy';
        case {'infinite','global'}
          out={strrep(filename,'.<TIMESTAMP>','')};
          return
        otherwise
          error([mfilename,': cannot handle metadata key ''storage_period'' with value ''',obj.mdget('storage_period'),'''.'])
        end
        %build list of files
        out=cell(size(timestamp_list));
        for i=1:numel(timestamp_list)
          out{i}=strrep(filename,'<TIMESTAMP>',datestr(timestamp_list(i),timestamp_fmt));
        end
      else
        out={filename};
      end
    end
    function out=isfile(obj,varargin)
      file_list=obj.file(varargin{:});
      out=cellfun(@(i)(~isempty(dir(i))),file_list);
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
    %% metadata file
    function out=md_file(obj,dir)
      if ~exist('dir','var') || isempty(dir)
        dir=obj.metadata_dir;
      end
      %get non-empty filename parts
      filename=obj.dataname.cells_clean;
      filename{end+1}='metadata';
      %assemble components
      filename=strjoin(filename,'.');
      %add path
      out=fullfile(dir,filename);
    end
    function out=ismd_file(obj)
      out=~isempty(dir(obj.md_file));
    end
    function md_file_check(obj)
      if ~obj.ismd_file
        error([mfilename,': could not find the metadata for product ',obj.dataname.name,' (expecting ',obj.md_file,').'])
      end
    end    
    function obj=metadata_load(obj)
      obj.md_file_check
      obj=obj.metadata_merge(ReadYaml(obj.md_file));
    end
    function obj=metadata_merge(obj,new_metadata)
      % entries in new_metadata replace those already in obj.metadata
      f=fieldnames(new_metadata);
      for i=1:numel(f)
        obj.metadata.(f{i})=new_metadata.(f{i});
      end
    end
    %% metadata fields
    function out=ismd_field(obj,metadatafieldname)
      out=isfield(obj.metadata,metadatafieldname);
    end
    function md_field_check(obj,metadatafieldname)
      if ~obj.ismd_field(metadatafieldname)
        error([mfilename,': cannot find field ',metadatafieldname,' in the metadata of ',obj.dataname.name,'.'])
      end
    end
    %% metadata parsing
    function out=mdget(obj,metadatafieldname)
      obj.md_field_check(metadatafieldname)
      in=obj.metadata.(metadatafieldname);
      switch class(in)
      case {'double','logical'}
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
    %% data sources
    function out=nr_sources(obj)
      if isfield(obj.metadata,'sources')
        out=numel(obj.metadata.sources);
      else
        out=0;
      end
    end
    function source_name=sources(obj,idx)
      if ~isnumeric(idx) || ~isfinite(idx) || ~isscalar(idx) 
        error([mfilename,': input ''idx'' must be a scalar integer, not ''',num2str(idx),'''.'])
      end
      if ~isfield(obj.metadata,'sources')
        error([mfilename,': product ',obj.name,' does not have any ''source'' field in it''s metadata.'])
      end
      if idx>obj.nr_sources
        error([mfilename,': cannot get the ',str.th(idx),' data source because there are only ',num2str(obj.nr_sources),'.'])
      end
      source_name=datanames(obj.metadata.sources{idx});
    end
    %% plot customization
    function enforce_plot(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('fig_handle',  gcf,  @(i) ishandle(i));
      p.addParameter('axis_handle', gca,  @(i) ishandle(i));
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
          strrep(mdf{i},'plot_',''),...
          obj.metadata.(mdf{i}),...
          @(i) iscell(i) || ischar(i) || isdouble(i) || islogical(i)... %pretty generic validation needed to let anything in
        )
      end
      % parse it
      p.parse(varargin{:});
      % sanity
      if any(isfinite(p.Results.ylimits)) && ( p.Results.autoscale ||  p.Results.automean )
        error([mfilename,': option ''ylimits'' and ''autoscale'' or ''automean'' do not work concurrently.'])
      end
      % enforce fontsize and paper size
      set(    p.Results.axis_handle,          'FontSize',p.Results.fontsize_axis)
      set(get(p.Results.axis_handle,'Title' ),'FontSize',p.Results.fontsize_title);
      set(get(p.Results.axis_handle,'XLabel'),'FontSize',p.Results.fontsize_label);
      set(get(p.Results.axis_handle,'YLabel'),'FontSize',p.Results.fontsize_label);
      set(    p.Results.fig_handle, 'Position',          p.Results.size,...
                                    'PaperUnits',        p.Results.units,...
                                    'PaperPosition',     p.Results.size);
      % enforce line properties
      line_handles=get(p.Results.axis_handle,'Children');
      for i=1:numel(line_handles)
        set(line_handles(i),'LineWidth',p.Results.line_width)
      end
      % enforce x-limits
      v=axis(p.Results.axis_handle);
      if p.Results.xdate
        for i=1:2
          if isfinite(   p.Results.xlimits(i))
            v(i)=datenum(p.Results.xlimits(i));
          end
        end
        if ~strcmp(datestr(v(1),'yyyymmdd'),datestr(v(2),'yyyymmdd')) && ...
            (~strcmp(datestr(v(2),'HHMMSS'),'000000') || v(2)-v(1)>1)
          xlabel([datestr(v(1),'yyyy-mm-dd'),' to ',datestr(v(2),'yyyy-mm-dd')])
        else
          xlabel(datestr(v(1)))
        end
      end
      % enforce reasonable y-limits
      for i=1:2
        if isfinite(p.Results.ylimits(i))
          v(i+2)=   p.Results.ylimits(i);
        end
      end
      % enforce auto-scale and/or auto-mean
      if p.Results.automean || p.Results.autoscale
        %gather plotted data
        dat=cell(size(line_handles));
        for i=1:numel(line_handles)
          dat{i}=get(line_handles(i),'ydata');
        end
        dat=[dat{:}];
        dat=dat(~isnan(dat(:)));
        %remove outliers (if requested)
        for c=1:p.Results.outlier
          dat=simpledata.rm_outliers(dat);
        end
        dat=dat(~isnan(dat(:)));
        %enfore data-driven mean and/or scale
        if p.Results.automean && p.Results.autoscale
          v(3:4)=mean(dat)+4*std(dat)*[-1,1];
        elseif p.Results.automean
          v(3:4)=mean(dat)+0.5*diff(v(3:4))*[-1,1];
        elseif p.Results.autoscale
          v(3:4)=mean(v(3:4))+4*std(dat)*[-1,1];
        end
        %fix non-negative data sets
        if all(dat>=0) && v(3)<0
          v(3)=0;
        end
      end
      axis(p.Results.axis_handle,v);
    end
  end
end