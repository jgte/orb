classdef dataproduct
  %static
  properties(Constant)
    default_metadata=struct(...
      'plot_product',false...
    );
     parameter_list={...
      'metadata_dir',   file.orbdir('metadata'), @ischar;...
      'plot_dir',       file.orbdir(    'plot'), @ischar;...
      'data_dir',       file.orbdir(    'data'), @ischar;...
      'metadata',       dataproduct.default_metadata,               @isstruct;...
      'debug_plot',     false,                                      @islogical;...
      'start',          time.zero_date,                             @isdatetime;...
      'stop',           time.inf_date,                              @isdatetime;...
      'submetadataname','submetadata',                              @isdatetime;...
    };
    %this is the maximum depth of the structure inside on dataname, i.e. obj.data.(dataname.name)
    %this parameter is arbitrary but keep it as low as possible to prevent time-consuming searches.
    maxdepth=4;
  end
  %read-write
  properties
    dataname
    metadata_dir
    plot_dir
    data_dir
    metadata
    debug_plot
    %%these are methods
    %start
    %stop
    submetadataname
  end
  %calculated only when asked for
  properties(Dependent)
    name
    codename
    mdfile
  end
  properties (GetAccess=private)
    mdfilei
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(dataproduct.parameter_list); end
      out=v.picker(varargin{:});
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
              product_patched.dataname=products_in{i}.dataname.append_field_leaf(...
                str.clean([level_name,'_',str.show(level_vals{p})],'fieldname')...
              );
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
    function out=array(in,varargin)
      assert(iscell(in),[mfilename,': cannot handle input ''in'' of class ',class(in),', expecting a cell array.'])
      if isempty(in)
        out={};
      else
        in=cells.scalar(in,'set');
        out=cellfun(@(i)dataproduct(i,varargin{:}),in,'UniformOutput',false);
      end
    end
  end
  methods
    %% constructor
    function obj=dataproduct(in,varargin)
      %sanity
      assert(~isempty(in),'cannot handle empty input ''in''.')
      assert(isscalar(in)||ischar(in),'cannot handle non-scalar inputs; consider using dataproduct.array')
      % trivial call (all other classes are handled in datanames)
      if isa(in,'dataproduct')
        obj=in;
        return
      end
      % input parsing (relevant to the parameters defined in dataproduct.parameter_list)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('field_path',    '',@(i) ischar(i) || iscell(i));
      p.addParameter('metadata_from', '',@isstruct);
      %create argument object, declare and parse parameters, save them to obj
      [~,p,obj]=varargs.wrap('sinks',{obj},'parser',p,'sources',{dataproduct.parameters('obj')},varargin{:});
      %sanity
      assert(file.exist(obj.metadata_dir),[mfilename,': ',...
        'cannot find metadata dir ''',obj.metadata_dir,'''.'])
      % call superclass constructor
      if ~isempty(p.Results.field_path)
        obj.dataname=datanames(in,p.Results.field_path);
      else
        obj.dataname=datanames(in);
      end
      %get the metadata name
      metadataname=strsplit(obj.name,file.build_element_char);
      %save mdfile
      obj.mdfile=strjoin(metadataname,file.build_element_char);
      %make sure it exists, if not start digging parent metadata names 
      %NOTICE: this facility is probably unused, not sure what it is for
      while ~obj.ismdfile && numel(metadataname)>1
        %ditch last part of metadata name
        metadataname(end)=[];
        %try this metadata name
        obj.mdfile=strjoin(metadataname,file.build_element_char);
      end
      %really need that metadata
      obj.mdfile_check;
      %inform if the metadata file had to be dug up
      if numel(strsplit(obj.name,file.build_element_char))~=numel(metadataname)
        str.say('NOTICE: could not find metadata of product ',obj.dataname.str,...
          ' but found metadata of product ',strjoin(metadataname,file.build_element_char),' (using the latter).')
      end
      % load metadata
      obj=obj.mdload(p.Results.metadata_from);
      % input parsing (relevant to the parameters defined in obj.metadata)
      [~,~,obj.metadata]=varargs.wrap('sources',{obj.default_metadata,obj.metadata},'sinks',{obj.metadata},varargin{:});
    end
    function obj=rm_data(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('rm_args',[], @ischar);
      p.addParameter('mode','data', @ischar);
      % parse it
      p.parse(varargin{:});
      %get list of files
      files=obj.file(p.Results.mode,varargin{:},'discover',true);
      %loop over all files and delete them
      for i=1:numel(files)
        if file.exist(files{i})
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
          'dir',obj.data_dir,...
          'sub_dirs','single',...
          'ext','mat'...
        };
      case 'plot'
        out={...
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
      p.addParameter('start_timestamp_only', true, @islogical);
      p.addParameter('use_storage_period',   true, @islogical);
      p.addParameter('discover',            false, @islogical);
      p.addParameter('no_extension',        false, @islogical);
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
        %first handle storage periods with infinite span
        switch lower(obj.mdget('storage_period'))
        case {'infinite','global'}
          filename=strrep(filename,'.<TIMESTAMP>','');
          startlist=time.zero_date;
          stoplist =time.inf_date;
          timestamp_fmt='none';
        case {'none'}
          filename='';
          startlist=time.zero_date;
          stoplist =time.zero_date;
          timestamp_fmt='none';
        otherwise
          %check if stop date is compatible
          assert(isfinite(p.Results.stop),'Trying to return infinite number of files, better set explicit start/stop')
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
          case {'direct'}
            startlist=p.Results.start;
             stoplist=p.Results.stop;
             timestamp_fmt='yyyymmdd';
          otherwise
            error([mfilename,': cannot handle metadata key ''storage_period'' with value ''',obj.mdget('storage_period'),'''.'])
          end
        end
        %build list of files
        out=cell(size(startlist));
        if str.none(timestamp_fmt)
          assert(numel(startlist)==1,'BUG TRAP: timestamp_fmt is ''none'' but there are multiple entries in ''startlist''.')
          out{1}=filename;
        else
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
        end
      else
        out={filename};
      end
      %enforce no extension option
      if p.Results.no_extension; out=file.remove_ext(out); end
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
          if ~file.exist(f)
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
    function out=get.mdfile(obj)
      % build filename and add path
      out=obj.mdfilei;
    end
    function obj=set.mdfile(obj,metadataname)
      %NOTICE: although metadataname is (usually) obj.dataname.name, don't
      %        set obj.dataname explicity to metadataname, e.g.:
      %        obj.dataname=datanames(metadataname);
      %        It's healthy to keep the mdfile and dataname members disconnected.
      %NOTICE: don't check for the existence of the metadata file, that is
      %        done at init
      % build filename and add path
      obj.mdfilei=fullfile(obj.metadata_dir,[metadataname,'.yaml']);
    end
    function out=ismdfile(obj)
      out=file.exist(obj.mdfile);
    end
    function mdfile_check(obj)
      assert(obj.ismdfile,[mfilename,': ',...
        'could not find the metadata for product ',obj.name,' (expecting ',obj.mdfile,').'])
    end
    function obj=mdload(obj,metadata)
      if ~exist('metadata','var') || isempty(metadata)
        %need to read YAML (load packages dir, inside which the +yaml package sits)
        addpath(fullfile(file.orbdir('packages')));
        %make sure the metadata files is there
        obj.mdfile_check
        %load metadata
        metadata=yaml.ReadYaml(obj.mdfile);
        %convert numeric cell arrays to normal arrays
        metadata=structs.fun(@cells.c2m,metadata);
      end
      %load and merge with default/existing metadata (this needs to come before handling sub-metadata files, so that these entries actually exist in the metadata)
      obj=obj.mdmerge(metadata);
      %init loop variables (for the submetadatafilelist, assume some size, doesn't really matter)
      c=0;submetadatafilelist=cell(0);
      %load sub-metadata 
      while true
        if obj.ismdfield(obj.submetadataname)
          %build the submetadata filename
          submetadatafiles=cellfun(@(i) fullfile(obj.metadata_dir,[i,'.yaml']),cells.scalar(obj.mdget(obj.submetadataname),'set'),'UniformOutput',false);
          %delete the submetadata field
          obj.metadata=rmfield(obj.metadata,obj.submetadataname);
          %loop over all 
          for i=1:numel(submetadatafiles)
            %check if this submetadatafile has not already been saved
            if ~ismember(submetadatafiles{i},submetadatafilelist)
              %sanity
              assert(file.exist(submetadatafiles{i}),['Cannot find sub-metadatafile ',submetadatafiles{i},'.'])
              %increment counter
              c=c+1;
              %add this metadata file to the list
              submetadatafilelist{c}=submetadatafiles{i};
              %load and merge this submetadatafile         
              obj=obj.mdmerge(yaml.ReadYaml(submetadatafiles{i}));
            end
          end
        else
          break
        end
      end
      %if there is metadata, do some more stuff 
      if ~isempty(submetadatafilelist)
        %re-load metadata so that any metadata in the main metadata file over-writes that in the sub-metadata files
        obj=obj.mdmerge(metadata);
        %save submetadatafilelist (deletes the original submetadata field, but updates any deep dependencies)
        obj.metadata.(obj.submetadataname)=submetadatafilelist;
      end
      %parse commands (if any)
      obj.metadata=dataproduct.parse_commands(obj.metadata);
    end
    function obj=mdmerge(obj,new_metadata)
      obj.metadata=structs.copy(new_metadata,obj.metadata);
      %need to translate sources: from char to scalar cellstr and convert them to datanames
      if isfield(obj.metadata,'sources')
        obj.metadata.sources=datanames.array(obj.metadata.sources);
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
      p.addParameter('return_empty_if_missing',false,@islogical);
      p.addParameter('always_cell_array',      false,@islogical);
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
    %% resolve some important metadata
    function out=force(obj,force_in)
      out=(force_in && ~obj.mdget('never_force','default',false)) || obj.mdget('always_force','default',false);
      if (force_in && obj.mdget('never_force','default',false))
        str.say('stack_delta',1,'WARNING: Force is true as input argument but product',obj,...
          'has the never_force property: force as input argument ignored')
      end
    end
    %% start/stop metadata wrappers
    function out=start(obj)
      if obj.ismdfield('start')
        out=obj.mdget('start');
      else
        out=dataproduct.parameters('start');
      end
    end
    function out=stop(obj)
      if obj.ismdfield('stop')
        out=obj.mdget('stop');
      else
        out=dataproduct.parameters('stop');
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
    end
    function out=source_list_str(obj)
      out=cellfun(@(i) i.str,obj.source_list,'UniformOutput',false);
    end
    %% plot customization
    function out=plot_args(obj)
      %get plot parameters already defined in the metadata
      mdf=fieldnames(obj.metadata);
      %count how many plot_* field there are
      idx=cells.strfind(mdf,'plot_');
      %init outputs
      out=cell(1,2*numel(idx));
      %loop over all of them
      for i=1:numel(idx)
        %add this plot parameter to output
        out{2*i-1}=mdf{idx(i)};
        out{2*i  }=obj.metadata.(mdf{idx(i)});
      end
    end
    function out=enforce_plot(obj,varargin)
      %get plot parameters already defined in the metadata and call mother routine
      out=plotting.enforce(obj.plot_args{:},varargin{:});
    end
    %% argument parsing
    function out=args(obj,field_list)
      %NOTICE: need to remove some fields so it doesn't conflict with varargs
      rm_field_list=intersect(fieldnames(obj.metadata),varargs.reserved_fields);
      out=structs.varargin(rmfield(obj.metadata,rm_field_list));
      %return only the requested fields, if any
      if exist('field_list','var')
        out=cells.varargin('get',out,field_list);
      end
    end
    %% level-wrapping handlers
    function out=islevel_wrapped(obj,level)
      out=obj.ismdfield(dataproduct.level_vals_str(level)) && obj.ismdfield(dataproduct.level_name_str(level));
    end
    function out=is_wrapped(obj)
      out=obj.islevel_wrapped(1);
    end
    function out=level_vals(obj,level)
      out=obj.mdget(dataproduct.level_vals_str(level));
      switch class(out)
      case 'datetime'
        out=arrayfun(@(i) {i},out);
      case 'double'
        out=cells.m2c(out);
      otherwise  
        out=cells.scalar(out,'set');
      end
%       %need level_vals to be a cell array
%       assert(iscell(out),['BUG TRAP: need metadata entry ''',dataproduct.level_vals_str(level),...
%         ''' to be a cell array, not a ',class(out),'.'])
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
