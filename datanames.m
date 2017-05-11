classdef datanames
  %static
  properties(Constant)
    separator='.';
    separator_clean='__'; %needed to construct legal field names out of obj.name
    parts={...
      'type',...
      'level',...
      'field',...
      'sat'...
    };
  end
  %read-write
  properties
  end
  %read only
  properties(SetAccess=private)
    category
    type
    level
    field
    sat
  end
  %internal
  properties(SetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
    name
    cells
    length
    leaf_length
    trunk_length
  end
  methods(Static)
    function array(in)
      if ~iscell(in)
        error([mfilename,': cannot handle input ''in'' of class ',class(in),', expecting a cell array.'])
      end
      out=cell(size(in));
      for i=1:numel(in)
        out{i}=datanames(in{i});
      end
    end
    function out=build_name(cells)
      out=strjoin(cells,datanames.separator);
    end
    function out=build_name_clean(cells)
      out=strjoin(cells,datanames.separator_clean);
    end
    function out=common(datanames_array)
      out=datanames_array{1}.cells;
      for i=2:numel(datanames_array)
        out=intersect(out,datanames_array{i}.cells,'stable');
      end
    end
    function out=unique(datanames_array)
      com=datanames.common(datanames_array);
      out=cell(size(datanames_array));
      for i=1:numel(datanames_array)
        out{i}=setdiff(datanames_array{i}.cells,com,'stable');
      end
    end
  end
  methods
    %% constructor
    function obj=datanames(in)
      % init object
      obj=obj.clear;
      % branch on class
      switch class(in)
      case 'datanames'
        obj=datanames(in.cells);
      case 'cell'
        obj.cells=in;
      case 'char'
        obj.name=in;
      case 'dataproduct'
        obj=in.dataname;
      otherwise
        error([mfilename,': cannot handle input ''in'' of class ',class(in),'.'])
      end
    end
    function obj=clear(obj)
      obj.category='';
      obj.type='';
      obj.level='';
      obj.field='';
      obj.sat='';
    end
    %% name operations
    function out=get.name(obj)
      %get non-empty data name parts and assemble components
      out=datanames.build_name(obj.cells_clean);
    end
    function out=name_clean(obj)
      %NOTICE: the meaning of 'clean' is not the same:
      % - obj.cells_clean returns non-empty cells and
      % - obj.build_name_clean returns strings with obj.separator replaced with obj.separator_clean
      out=datanames.build_name_clean(obj.cells_clean);
    end
    function obj=set.name(obj,in)
      if isempty(strfind(in,obj.separator));
        error([mfilename,': this does not seem to be a valid data storage name: ''',in,'''.'])
      end
      obj.cells=strsplit(in,obj.separator);
    end
    %% cells operations
    function out=get.cells(obj)
      %assign outputs
      out={...
        obj.category,...
        obj.type,    ...
        obj.level,   ...
        obj.field,   ...
        obj.sat      ...
      };
    end
    function out=cells_clean(obj)
      %outputs
      out=obj.cells;
      %enforce cleaning cells
      out=out(~cellfun(@isempty,out));
    end
    function out=get.length(obj)
      out=numel(obj.cells_clean);
    end
    function obj=set.cells(obj,in)
      %sanity
      if ~iscellstr(in) || numel(in)<2
        error([mfilename,': ',...
          'input ''in'' must be a cell of strings with length at least 2, '...
          'not a ',class(in),' of ',class(in{1}),' with length ',num2str(numel(in)),'.'...
        ])
      end
      obj.category=in{1};
      obj.type=    in{2};
      if numel(in)>=3; obj.level=in{3}; end
      if numel(in)>=4; obj.field=in{4}; end
      if numel(in)>=5; obj.sat=  in{5}; end
    end
    %% leaf operations: returns the single outter most category, type, level, field or sat (down by delta)
    function out=leaf(obj,delta)
      if ~exist('delta','var') || isempty(delta)
        delta=0;
      end
      out=obj.cells;
      out=out(obj.length-delta);
    end
    function out=leaf_clean(obj)
      %NOTICE: this is similar to obj.cells_clean but the empty cells are only removed
      %        from the outtermost category, type, level, field or sat, up to the 
      %        *outtermost* non-empty cell.
      in=obj.cells;
      %outputs
      out=in;
      %trim empty cells, starting from the outtermost
      for i=numel(in):-1:1
        if isempty(in{i})
          out=in(1:i-1);
        else
          return
        end
      end
    end
    function out=get.leaf_length(obj)
      out=numel(obj.leaf_clean);
    end
    function out=leaf_type(obj,delta)
      if ~exist('delta','var') || isempty(delta)
        delta=0;
      end
      switch obj.length-delta
      case 1; out='category';
      case 2; out='datatype';
      case 3; out='level';
      case 4; out='field';
      case 5; out='sat';
      otherwise
        error([...
          mfilename,': ilegal input ''delta'' (with value ',num2str(delta),...
          '), it cannot exceed ',num2str(obj.length),'.'...
        ])
      end
    end
    %% trunk operations: returns all up to the outtermost non-empty category, type, level, field or sat (down by delta)
    function out=trunk(obj,delta)
      if ~exist('delta','var') || isempty(delta)
        delta=0;
      end
      out=obj.trunk_clean;
      out=out(1:obj.trunk_length-delta);
    end
    function out=trunk_clean(obj)
      %NOTICE: this is similar to obj.cells_clean but the empty cells are only removed
      %        from the outtermost category, type, level, field or sat, up to the 
      %        *innermost* non-empty cell.
      in=obj.cells;
      %outputs
      out=cell(0);
      %trim empty cells, starting from the outtermost
      for i=1:numel(in)
        if ~isempty(in{i})
          out=in(1:i);
        else
          return
        end
      end
    end
    function out=get.trunk_length(obj)
      out=numel(obj.trunk_clean);
    end
    function out=trunk_name(obj,delta)
      out=datanames.build_name(obj.trunk(delta));
    end
    function out=trunk_name_clean(obj,delta)
      out=datanames.build_name_clean(obj.trunk(delta));
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
      p.addParameter('full_path',         true,    @(i) islogical(i))
      % parse it
      p.parse(varargin{:});
      %get data file/path parts (and remove the specified filename part, N.B.: empty parts are handled in obj.edit)
      dataparts=obj.edit(p.Results.remove_part,'').cells_clean;
      %propagate to filename
      filename=dataparts;
      %add prefix and suffix (is non-empty)
      if ~isempty(p.Results.prefix); filename=[{p.Results.prefix},filename]; end
      if ~isempty(p.Results.suffix); filename=[filename,{p.Results.suffix}]; end
      %add time stamp placeholder
      filename{end+1}='<TIMESTAMP>';
      %add extension
      if ~isempty(p.Results.ext)
        filename{end+1}=p.Results.ext;
      end
      %assemble components and add path
      if p.Results.full_path
        out=fullfile(p.Results.dir,dataparts{:},strjoin(filename,datanames.separator));
      else
        out=fullfile(p.Results.dir,strjoin(filename,datanames.separator));
      end
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
    %% direct editing
    function obj=edit(obj,part,value)
      if ~isempty(part)
        obj.(part)=value;
      end
    end
    %% legend strings
    function out=legend(obj,prefix,suffix)
      if ~exist('prefix','var') || isempty(prefix)
        prefix='';
      end
      if ~exist('suffix','var') || isempty(suffix)
        suffix='';
      end
      switch obj.category
      case 'grace'
        switch obj.type
        case 'calpar_csr'
          out={obj.level,obj.field,['GRACE-',obj.sat]};
        otherwise
          out={obj.type,obj.level,obj.field,['GRACE-',obj.sat]};
        end
      otherwise
        out={obj.category,obj.type,obj.level,obj.field,obj.sat};
      end
      out=strjoin([{prefix},out,{suffix}],' ');
    end
    %% operator overloading
    function out=eq(obj1,obj2)
      out=all(cellfun(@(i) strcmp(obj1.(i),obj2.(i)),datanames.parts));
    end
  end
end