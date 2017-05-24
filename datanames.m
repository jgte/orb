classdef datanames
  %static
  properties(Constant)
    separator='.';
    separator_clean='_'; %needed to construct legal field names out of obj.name
  end
  %read-write
  properties
  end
  %read only
  properties(SetAccess=private)
    name
    field_path
  end
  %internal
  properties(SetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
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
  end
  methods
    %% constructor
    function obj=datanames(in,field_path)
      %sanity
      assert(~isempty(in),'cannot handle empty input ''in''.')
      if ischar(in)
        obj.name=in;
      % handle non-scalar inputs
      elseif ~isscalar(in) 
        error([mfilename,': cannot handle non-scalar inputs.'])
      else
        % branch on class
        switch class(in)
        case 'datanames'
          obj=in;
        case 'cell'
          obj=datanames(in{1});
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
      %some sanity
      assert(ischar(field_path_leaf),['input ''field_path_leaf'' must be a string, not a ',class(field_path_leaf),'.'])
      obj.field_path=[obj.field_path,{field_path_leaf}];
    end
    function out=global_field_path(obj)
      out=[{obj.name_clean},obj.field_path];
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
    function obj=edit_field_part(obj,field_part,value)
      idx=strcmp(obj.field_path,field_part);
      if any(idx)
        obj.field_path(idx)={value};
      end
    end
    function out=field_path_str(obj)
      out=strjoin(obj.field_path,'.');
    end
    function out=split(obj)
      out=cells.flatten({strsplit(obj.name,datanames.separator),obj.field_path});
    end
    %% names
    function out=filename(obj)
      out=fullfile(obj.global_field_path{:});
    end
    function out=codename(obj)
      out=['(''',strjoin(obj.global_field_path,''').('''),''')'];
    end
    function out=str(obj)
      out=fullfile(obj.name,obj.field_path{:});
    end
    %% name operations
    function out=name_clean(obj)
      out=strrep(obj.name,'.',datanames.separator_clean);
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
      p.addParameter('remove_part',       '',      @(i) ischar(i))
      p.addParameter('prefix',            '',      @(i) ischar(i))
      p.addParameter('suffix',            '',      @(i) ischar(i))
%       p.addParameter('full_path',         true,    @(i) islogical(i))
      % parse it
      p.parse(varargin{:});
      %propagate dataname to filename
      filename={obj.filename};
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
        if isempty(dir(fileparts(out)))
          mkdir(fileparts(out))
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
      out=strcmp(obj1.name,obj2.name);
    end
  end
end