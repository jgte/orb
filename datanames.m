classdef datanames
  %static
  properties(Constant)
    separator='.';
    separator_clean='_'; %needed to construct legal field names out of obj.name
  end
  %read only
  properties(SetAccess=private)
    name
    field_path
  end
  methods(Static)
    function out=array(in,field_path)
      assert(iscell(in),[mfilename,': cannot handle input ''in'' of class ',class(in),', expecting a cell array.'])
      if isempty(in)
        out={};
      else
        if exist('field_path','var') && ~isempty(field_path)
          out=cellfun(@(i) datanames(i,field_path),in,'UniformOutput',false);
        else
          out=cellfun(@(i) datanames(i),in,'UniformOutput',false);
        end
      end
    end
    function out=common(dn_array)
      %this spits out a vertical cell array
      out=dn_array{1}.split;
      for i=2:numel(dn_array)
        out=intersect(out,dn_array{i}.split,'stable');
      end
    end
    function out=unique(dn_array)
      %this spits out a vertical cell array
      common=datanames.common(dn_array);
      out=cell(size(dn_array));
      for i=1:numel(dn_array)
        out{i}=setdiff(dn_array{i}.split,common,'stable');
      end
    end
    function out=transmute(in) 
      if iscell(in)
        out=cellfun(@datanames.transmute,in,'UniformOutput',false);
      else
        out=in.dataname;
      end
    end
  end
  methods
    %% constructor
    function obj=datanames(in,field_path)
      %sanity
      assert(~isempty(in),'cannot handle empty input ''in''.')
      if ischar(in)
        if ~isempty(strfind(in,filesep))
          in=strsplit(in,filesep);
          obj.name=in{1};
          obj=obj.set_field_path(in(2:end));
        else
          obj.name=in;
        end
      % handle non-scalar inputs
      elseif ~isscalar(in) 
        error([mfilename,': cannot handle non-scalar inputs.'])
      else
        % branch on class
        switch class(in)
        case 'datanames'
          obj=in;
        case 'cell'
          obj=datanames(in{1}); %already known to be scalar
        case 'dataproduct'
          obj=in.dataname;
        otherwise
          error([mfilename,': cannot handle input ''in'' of class ',class(in),'.'])
        end
      end
      %add field_path, if there
      if exist('field_path','var') && ~isempty(field_path)
        %propagate
        obj=obj.set_field_path(field_path);
      else
        %enforce the right kind of empty
        if isempty(obj.field_path)
          obj.field_path={};
        end
      end
    end
    %% field_path operations
    function obj=append_field_leaf(obj,field_path_leaf)
      %some trivial translation
      if isscalar(field_path_leaf) && iscell(field_path_leaf); field_path_leaf=field_path_leaf{1}; end
      %some sanity
      assert(ischar(field_path_leaf),['input ''field_path_leaf'' must be a string, not a ',class(field_path_leaf),'.'])
      obj.field_path=[obj.field_path,{field_path_leaf}];
    end
    function obj=prepend_field_root(obj,field_path_root)
      %some trivial translation
      if isscalar(field_path_root) && iscell(field_path_root); field_path_root=field_path_root{1}; end
      %some sanity
      assert(ischar(field_path_root),['input ''field_path_root'' must be a string, not a ',class(field_path_root),'.'])
      obj.field_path=[{field_path_root},obj.field_path];
    end
    function obj=set_field_path(obj,field_path)
      %handle shallow structures
      if ischar(field_path)
        field_path={field_path};
      end
      %some sanity
      assert(iscellstr(field_path),['input ''field_path'' must be a cell of strings, not a ',class(field_path),'.'])
      %assign
      obj.field_path=field_path;
    end
    %input 'field_parts' can be a cell string, with possible values for a particular part.
    function obj=edit_field_part(obj,field_parts,value)
      %translate inputs
      if ~iscell(field_parts)
        field_parts={field_parts};
      end
      %sanity on type
      assert(iscellstr(field_parts),'Can only handle cell strings.')
      %loop over all current parts
      for i=1:numel(obj.field_path)
        %check if this part shows up in any entries of 'field_parts'
        if any(strcmp(obj.field_path{i},field_parts))
          %if so, assign input 'value' to it and bail
          obj.field_path(i)={value};return
        end
      end
    end
    function out=split(obj)             %e.g. {'grace'    'calpar'    'csr'    'estim'    'AC0X'    'A'}
      out=cells.flatten({strsplit(obj.name,datanames.separator),obj.field_path});
    end
    function out=global_field_path(obj) %e.g. {'grace_calpar_csr'    'estim'    'AC0X'    'A'}
      out=[{obj.name_clean},obj.field_path];
    end
    %% names
    function out=filename(obj)       %e.g. grace.calpar.csr
%NOTICE: filename used to be like this but i don't think it makes much sense: the data should all be contained inside the
%        the file with name equalt to the product name, which is the same as the metadata name. The advantage is that a
%        product can be referred from a sub field and still the correct file is picked.
%       out=strjoin([{obj.name},obj.field_path],datanames.separator);
      out=obj.name;
    end
    function out=codename(obj)       %e.g. grace_calpar_csr.estim.AC0X.A
      out=strjoin(obj.global_field_path,'.');
    end
    function out=dotname(obj)        %e.g. grace.calpar.csr.estim.AC0X.A
      out=strjoin([{obj.name};obj.field_path(:)],'.');
    end
    function out=str(obj)            %e.g. grace.calpar.csr/estim/AC0X/A
      out=fullfile(obj.name,obj.field_path{:});
    end
    function out=field_path_str(obj) %e.g.                  estim.AC0X.A
      out=strjoin(obj.field_path,'.');
    end
    %% name operations
    function out=name_clean(obj)
      out=str.clean(obj.name,'fieldname',datanames.separator_clean);
    end
    %% filename builders
    function out=file(obj,varargin)
      %NOTICE: this procedure ALWAYS returns a string!
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('start',     datetime('now'), @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',      datetime('now'), @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('ext',       '',              @(i) ischar(i));
      p.addParameter('dir',       '',              @(i) ischar(i));
      p.addParameter('timestamp',         false,   @(i) islogical(i));
      p.addParameter('keeptsplaceholder', false,   @(i) islogical(i));
      p.addParameter('ensure_dir',        true,    @(i) islogical(i));
      p.addParameter('remove_part',       '',      @(i) ischar(i) || iscellstr(i))
      p.addParameter('prefix',            '',      @(i) ischar(i))
      p.addParameter('suffix',            '',      @(i) ischar(i))
      p.addParameter('sub_dirs',          'none',  @(i)ischar(i))
      % parse it
      p.parse(varargin{:});
      %propagate dataname to filename
      switch p.Results.sub_dirs
      case 'none'
        filename={obj.filename};
      case 'single'
        filename={fullfile(obj.filename,obj.filename)};
      case 'deep'
        filename=strsplit(obj.filename,datanames.separator);
        filename={fullfile(filename{:},obj.filename)};
      otherwise
        error(['Cannot understand input ''sub_dirs'' with value ''',p.Results.sub_dirs,'''.'])
      end
      %add prefix and suffix (if non-empty)
      if ~isempty(p.Results.prefix); filename=[{p.Results.prefix},filename]; end
      if ~isempty(p.Results.suffix); filename=[filename,{p.Results.suffix}]; end
      %add time stamp placeholder
      filename{end+1}='<TIMESTAMP>';
      %add extension
      if ~isempty(p.Results.ext)
        filename{end+1}=p.Results.ext;
      end
      %assemble components and add path
      out=fullfile(p.Results.dir,strjoin(filename,datanames.separator));
      %include detailed timestamp if requested
      if p.Results.timestamp
        out=strrep(...
          out,'<TIMESTAMP>',[...
            datestr(p.Results.start,30),'_',datestr(p.Results.stop,30)...
          ]...
        );
      elseif ~p.Results.keeptsplaceholder
        out=strrep(out,'.<TIMESTAMP>','');
      end
      %make sure dir exists, if requested
      if p.Results.ensure_dir
        dir_now=fileparts(out);
        if ~isempty(dir_now) && exist(dir_now,'dir')==0
          warning('off','MATLAB:MKDIR:DirectoryExists')
          mkdir(fileparts(out))
          warning('on','MATLAB:MKDIR:DirectoryExists')
        end
      end
    end
    function out=isfile(obj,varargin)
      out=~isempty(dir(obj.file(varargin{:})));
    end
    function filecheck(obj,varargin)
      if ~obj.isfile(varargin{:})
        error([mfilename,': for product ',obj.name,', could not find the following file: ',obj.file(varargin{:})])
      end
    end
    %% operator overloading
    function out=eq(obj1,obj2)
      out=strcmp(datanames(obj1).name,datanames(obj2).name);
    end
  end
end
