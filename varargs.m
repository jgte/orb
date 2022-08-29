classdef varargs < dynamicprops
  properties(Constant,GetAccess=private)
    %The name and value fields must be the first 2 defined below (they are never optional)
    template=struct(...
      'name',      [],...
      'value',     [],...
      'validation',@(i)true...
    );
    template_validation={...
      @ischar;...
      @(i) true;...
      @(i) isa(i, 'function_handle');...
    };
    nprops=numel(fieldnames(varargs.template));
  end
  properties(Constant,GetAccess=public)
    template_fields=fieldnames(varargs.template)
    reserved_fields={'parser','mandatory','sources','sinks'};
  end
  properties(GetAccess=public,SetAccess=private)
    S
    rest
  end
  methods(Static)
    %% S methods
    function out=wrongS(S)
      if ~isstruct(S)
        out=-1; return
      end
      for i=1:2
        if ~isfield(S,varargs.prop(i))
          out=i;return
        end
      end
      out=0;
    end
    function out=isS(S)
      out=varargs.wrongS(S)==0;
    end
    function assertS(S)
      switch varargs.wrongS(S)
      case 0
        %do nothing
      case -1
        error(['Need input ''S_now'' to be a structure, not a ''',class(S),'''.'])
      otherwise
        error(['Need structure ''S_now'' to have field ''', varargs.prop(wrongS(S)),'''.'])
      end
    end
    function S=initS(name,value,varargin)
      %need list of field names in the template
      template_fieldnames=fieldnames(varargs.template);
      %init parser
      p=machinery.inputParser;
      %declare (mandatory and optional) input parameters
      for i=1:numel(template_fieldnames)
        switch template_fieldnames{i}
        case {'name','value'}
          p.addRequired(template_fieldnames{i},varargs.template_validation{i});
        otherwise
          p.addParameter(...
            template_fieldnames{i},...                    %field name (e.g. 'validation')
            varargs.template.(template_fieldnames{i}),... %default value for this field
            varargs.template_validation{i});              %validation for this field
        end
      end
      %parse it!
      p.parse(name,value,varargin{:});
      %init Sout
      S=varargs.template;
      %propagate values caught by the parser
      for i=1:numel(template_fieldnames)
        S.(template_fieldnames{i})=p.Results.(template_fieldnames{i});
      end
    end
    %% is* methods
    function out=isvarargin(in)
      out=mod(numel(in),2)==0 && size(in,1)==1 && iscellstr(in(1:2:end));
    end
    function out=iscell(in)
      out=size(in,2)==varargs.nprops && iscellstr(in(:,1));
    end
    function out=isvarargs(in)
      out=...
        isempty(in)            || ...
        varargs.isS(in)        || ...
        varargs.isvarargin(in) || ...
        varargs.iscell(in)     || ...
        isstruct(in)           || ...
        isa(in,'varargs');
    end
    %% utilities
    %returns the a/all fieldname/s of varargs.template
    function out=prop(i)
      out=fieldnames(varargs.template);
      if exist('i','var') && ~isempty(i)
        out=out{i};
      end
    end
    %returns the positional index in varargin of the requested parameters (in cell array parameter_list)
    %returns 0 if not a parameter is not present
    function idx=locate(parameter_list,varargin)
      if ischar(parameter_list)
        parameter_list={parameter_list};
      end
      assert(iscellstr(parameter_list),['Input ''parameter_list'' must be a cell of strings, not a ''',...
        class(parameter_list),'''.'])
      idx=zeros(size(parameter_list));
      for i=1:numel(parameter_list)
        m=strcmp(varargin,parameter_list{i});
        if any(m)
          idx(i)=(find(m)+1)/2;
        else
          idx(i)=0;
        end
      end
    end
    %checks if a field name is reserved
    function out=isreserved(name)
      out=any(ismember(varargs.reserved_fields,name));
    end
    %% general wrapper for functions with only optional arguments
    function [v,p,sinks]=wrap(varargin)
      %The fields of v, p.Results and sinks{:} are defined from, in this order:
      % - 'sources': using (a cell array of) structs/cell matrices compatible with varargs' constructor
      %              (agregated and also returned in variable v);
      % - 'parser' : using matlab's parser object (also returned in variable p);
      % - varargin : cleaned of the parameters relevant to this method (see the code).
      %
      %NOTICE: only the parameters defined in 'sources' or 'parser' make their way into v; all other parameters in varargin are
      %ignored.
      %
      %NOTICE: If the same parameter is in 'sources' and varargin, the value defined in the latter is the one kept.
      %
      %NOTICE: the default parameters in 'parser' (specifically those that are not in varargin) are now passed to v (those
		  % that are defined in varargin are ignored).
      %
      %Both p and v outputs have the same information (except what concerns the 3rd notice above), which can be retrieved as:
      % - p.Results.(parameter_name) (as usual)
      % - v.Results.(parameter_name), v.value(parameter_name) or simply v.(parameter_name)
      %
      %The variables 'sinks' is a cell array with the objects passed in 'sinks' (in the same order), with their fields/methods
      %set with the (possible) values of any parameter passed in 'parser', 'sources' or varargin, as implemented in the
      %varargs.save methods.
      %this is the parser for this method (there's an additional parser going in and out if this method)
      pn=machinery.inputParser;
      pn.addParameter('parser',    inputParser, @(i) isa(i,'inputParser'));
      pn.addParameter('mandatory', {}, @iscell);
      pn.addParameter('sources',   {}, @(i) all(cellfun(@(j) varargs.isvarargs(j),i)));
      pn.addParameter('sinks',     {}, @iscell);
      pn.parse(varargin{:});
      %clean varargin
      varargin=cells.vararginclean(varargin,pn.Parameters);
      %retrieve (possible) sinks
      sinks=pn.Results.sinks;
      %retrieve sources
      sources=pn.Results.sources;
      %patch empty sources with whatever is in varargin
      if isempty(sources);sources={varargin};end
      %create first varargs obj
      v=varargs({});
      %loop over all sources and join to v (no entries are ignored, sources over-writes common entries in v)
      for i=1:numel(sources);v.join(sources{i});end
      %retrieve (possible) external parser
      p=pn.Results.parser; p.KeepUnmatched=true; p.PartialMatching=false; p.CaseSensitive=true;
      %the parameters already defined in the external parser are not to be saved to the 'sinks' (this only concerns the 'sinks')
      external_parameters=p.Parameters;
      %declare parameters in parser
      [p,v]=v.declare(p);
      %parse it
      p.parse(pn.Results.mandatory{:},varargin{:});
      %find p.UsingDefaults that are already part of v
      todelete=setdiff(p.UsingDefaults,setdiff(p.UsingDefaults,v.Parameters));
      %join parsed parameters (except those using default values already in v) into argument object v
      v=v.join(rmfield(p.Results,todelete));
      %go over all sinks (if there)
      for i=1:numel(sinks)
        %save all non-external parse parameters to this sink
        sinks{i}=v.save(sinks{i},external_parameters);
      end
      %convert to scalar if only one sink (isscalar({})=false)
      if isscalar(sinks)
        sinks=sinks{1};
      end
      %save rest: parameter pairs that were not save to object
      v.rest=cells.vararginclean(varargin,v.Parameters);
    end
  end
  methods
    %% constructor
    % Argument 'in' can be:
    % - structure scalar:
    %   - struct('name1','value1','name2','value2',...), as produced by obj.Results (no validation possible!);
    % - structure array:
    %   - S(:)=struct('name',name,'value',value,'validation',validation), as produced by obj.S;
    % - cell matrix:
    %   - {name1,value1,validation1;name2,value2,validation1;...}, as produced by obj.cell;
    % - cell array:
    %   - {name1,value1,name2,value2,...}, as produced by obj.varargin;
    function obj=varargs(in)
      assert(varargs.isvarargs(in),['Cannot handle input of class ''',class(in),'''.'])
      %handle empty inputs (defer to later)
      if isempty(in)
        return
      end
      %branch on class
      switch class(in)
      case 'varargs'
        %just propagate
        obj=in;
      case 'struct'
        if varargs.isS(in)
          obj.S=in;
          obj.update_dynamic_fields;
        else
          %parameter names are given in the field names of this structure
          fn=fieldnames(in);
          n=numel(fn);
          for i=1:n
            %propagate the name and value of this field
            obj.set(varargs.initS(fn{i},in.(fn{i})));
          end
        end
      case 'cell'
        if varargs.isvarargin(in)
          %need some reshaping (ony parameter name and value are retained)
          in=transpose(reshape(in,[2,numel(in)/2]));
        end
        n=size(in,1);m=size(in,2);
        %need at least two columns
        assert(m>=2,['If a cell array, input ''in'' must have at least two columns, not ',num2str(m),'.'])
        %need list of field names in the template
        template_fieldnames=fieldnames(varargs.template);
        %the paramenter names are given in the first column, the value in the second, etc... (as defined in varargs.template)
        for i=1:n
          %create cell array with optional fields of the template
          optargs=cell(1,m-2);
          for j=3:min([m,varargs.nprops])
            j1=(j-2)*2-1;
            j2=j1+1;
            optargs(j1)=template_fieldnames(j);
            optargs(j2)=in(i,j);
          end
          %propagate the name and value of this field, along with additional optional arguments (if not empty)
          obj.set(varargs.initS(in{i,1},in{i,2},optargs{:}));
        end
      end
    end
    %% representations
    function out=cell(obj)
      out=cell(obj.length,varargs.nprops);
      for i=1:obj.length
        for j=1:varargs.nprops
          out{i,j}=obj.get(i).(varargs.prop(j));
        end
      end
    end
    function out=Results(obj)
      for i=1:obj.length
        out.(obj.S(i).name)=obj.S(i).value;
      end
    end
    %function out=S(obj) (already a field)
    function out=varargin(obj,hide_rest_flag)
      if ~exist('hide_rest_flag','var') || isempty(hide_rest_flag)
        hide_rest_flag=false;
      end
      out=obj.cell;
      out=reshape(transpose(out(:,[1 2])),[1 obj.length*2]);
      if ~hide_rest_flag
        out=[out,obj.rest];
      end
    end
    %% utilities
    function out=size(obj,varargin)
      out=size(obj.S,varargin{:});
    end
    function out=length(obj)
      out=length(obj.S);
    end
    function out=isempty(obj)
      out=isempty(obj.S);
    end
    function l=tab(obj,fieldname)
      if ~exist('fieldname','var') || isempty(fieldname)
        fieldname='name';
      end
      switch fieldname
      case 'max'
        l=0;
        for i=1:varargs.nprops
          l=max([l,obj.tab(varargs.template_fields{i})]);
        end
      case varargs.template_fields
        %get length of longest parameter name/value/validation
        l=0;
        for i=1:obj.length
          l=max([length(obj.get(i).(fieldname)),l]);
        end
      otherwise
        error(['Cannot handle input ''fieldname'' with value ''',fieldname,'''.'])
      end
    end
    function [idx,c]=sort_idx(obj)
      c=obj.cell;
      [~,idx]=sort(c(:,1));
    end
    function out=str(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=obj.tab('max');
      end
      [sort_idx,c]=obj.sort_idx;
      out=cell(obj.length+1,1);
      out{1}=str.tablify(tab,varargs.prop);
      for i=1:obj.length
        out{i+1}=str.tablify(tab,c(sort_idx(i),:));
      end
      out=strjoin(out,newline);
    end
    function out=show(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=obj.tab('name');
      end
      %outputs
      out=cell(obj.length,1);
      sort_idx=obj.sort_idx;
      for i=1:obj.length
        p=obj.get(sort_idx(i));
        out{i}=str.show({str.tabbed(p.name,tab,true),p.value},'join_char',' : ');
      end
      out=strjoin(out,newline);
    end
    function out=Parameters(obj)
      %keep this a column vector, makes is easy to append other vectors
      out=cell(numel(obj.S),1);
      for i=1:numel(out)
        out{i}=obj.S(i).name;
      end
    end
    function out=picker(obj,mode,name)
      if ~exist('mode','var') || isempty(mode)
        mode='obj';
      end
      if iscellstr(mode)
        out=cell(size(mode));
        if exist('name','var')
          for i=1:numel(mode)
            out{i}=obj.picker(mode{i},name{i});
          end
        else
          out=cellfun(@obj.picker,mode,'UniformOutput',false);
        end
      else
        switch mode
          case 'length';    out=obj.length;
          case 'list';      out=obj.Parameters;
          case 'obj';       out=obj;
          case 'get';       out=obj.get(name);
          %these also make sense when name (the input var) is an idx (integer)
          case 'value';     out=obj.get(name).value;
          case 'name';      out=obj.get(name).name;
          case 'validation';out=obj.get(name).validation;
          otherwise;        out=obj.(mode);
        end
      end
    end
    function out=name(obj,parameters_to_skip)
      if ~exist('parameters_to_skip','var') || isempty(parameters_to_skip)
        parameters_to_skip={};
      end
      out=str.say('say_skip_stack',true,'say_join_char','.',dup(obj).delete(parameters_to_skip).varargin);
    end
    %% get methods
    %abstracts string or index
    function out=idx(obj,name)
      switch class(name)
      case 'char'
        out=find(strcmp(obj.Parameters,name));
      case 'double'
        assert(name>0,'If double, input ''name'' must be positive.')
        if name>obj.length
          out=[];
        else
          out=name;
        end
      end
    end
    function out=isparameter(obj,name)
      out=~isempty(obj.idx(name));
    end
    function out=get(obj,name)
      i=obj.idx(name);
      if isempty(i)
        out=varargs.template;
      else
        out=obj.S(i);
      end
    end
    function out=subset(obj,name_part)
      %dup the object
      out=dup(obj);
      %get the parameters
      p=out.Parameters;
      %get the indexes of the parameters to keep
      idx_to_keep=cellfun(@(i) str.contains(i,name_part),p);
      %check if there's anything to delete
      if sum(~idx_to_keep)>0
        %trim S
        out.S=out.S(idx_to_keep);
        %get idx of properties to delete
        idx_to_delete=find(~idx_to_keep);
        %delete the properties
        for i=1:numel(idx_to_delete)
          delete(findprop(obj,p{idx_to_delete(i)}));
        end
      end
    end
    %% set methods
    function obj=set(obj,Snew)
      %input checking
      varargs.assertS(Snew);
      %get the index of this parameters
      idx=obj.idx(Snew.name);
      %check if this is a new parameter
      new_parameter=isempty(idx);
      %if so append to the end
      if new_parameter
        idx=obj.length+1;
      end
      %need to init the structure in case this is the first field
      if obj.isempty
        obj.S=Snew;
      else
        obj.S(idx)=Snew;
      end
      %update dynamic fields
      obj.update_dynamic_fields(idx,new_parameter);
    end
    function out=isdynamic_defined(obj,name)
      try
        %if this works, then this name is already defined
        obj.(name);
        %set output
        out=true;
      catch
        %set output
        out=false;
      end
    end
    function obj=dynamic_set(obj,value,name)
      obj.S(obj.idx(name)).value=value;
    end
    function value=dynamic_get(obj,name)
      value=obj.S(obj.idx(name)).value;
    end
    function obj=update_dynamic_fields(obj,idx,new_parameter)
      if ~exist('idx','var') || isempty(idx)
        idx=1:obj.length;
      end
      if ~exist('new_parameter','var') || isempty(new_parameter)
        new_parameter=false;
      end
      %add dynamic names
      for i=1:numel(idx)
        name=obj.S(idx(i)).name;
        if new_parameter
          assert(~obj.isdynamic_defined(name),['Cannot handle parameter with name ''',name,...
            ''' because it conflcits with an already-defined (possibly inherited) method or property of this class.'])
        end
        if ~isprop(obj,name)
          h=obj.addprop(name);
          h.SetMethod=@(obj,i) obj.dynamic_set(i,name);
          h.GetMethod=@(obj  ) obj.dynamic_get(  name);
        end
      end
    end
    function obj=rm_empty(obj)
      to_delete=cell(1,obj.length);
      for i=1:obj.length
        if isempty(obj.S(i).value)
          to_delete{i}=obj.S(i).name;
        end
      end
      obj=obj.delete(cells.rm_empty(to_delete));
    end
    %% edit methods
    function obj=delete(obj,varargin) %deletes parameters given in varargin
      if numel(varargin)==1
        parameters=cells.scalar(cells.flatten(varargin{1}),'set');
      else
        parameters=varargin;
      end
      %sanity
      assert(iscellstr(parameters),['Can only handle cell array of strings, not ',class(varargin{1}),'.'])
      %check if there's nothing to delete ( intersect({},in.Parameters) is always empty )
      if isempty(intersect(parameters,obj.Parameters))
        return
      end
      %get the parameters that are not to be deleted
      parameters_to_keep=setdiff(obj.Parameters,parameters);
      %get the corresponding indexes
      idx_to_keep=cellfun(@(i) obj.idx(i),parameters_to_keep);
      %trim S
      obj.S=obj.S(idx_to_keep);
      %get the properties that are in obj
      idx=find(cellfun(@(i) isprop(obj,i),parameters));
      %remove the dynamic properties associated with the deleted parameters (if any)
      for i=1:numel(idx)
        delete(findprop(obj,parameters{idx(i)}));
      end
    end
    function obj=isolate(obj,varargin) %deletes all parameters except those given in varargin
      if numel(varargin)==1
        parameters=cells.scalar(cells.flatten(varargin{1}),'set');
      else
        parameters=varargin;
      end
      %sanity
      assert(iscellstr(parameters),['Can only handle cell array of strings, not ',class(varargin{1}),'.'])
      %get the parameters that are to be plucked
      parameters_to_delete=setdiff(obj.Parameters,parameters);
      %get the indexes of the parameters to keep
      idx_to_keep=cell2mat(cellfun(@(i) obj.idx(i),parameters,'UniformOutput',false));
      %trim S
      obj.S=obj.S(idx_to_keep);
      %remove the dynamic properties associated with the deleted parameters (if any)
      for i=1:numel(parameters_to_delete)
        delete(findprop(obj,parameters_to_delete{i}));
      end
    end
    function obj=rename(obj,old_name,new_name)
      assert(obj.isparameter(old_name),...
        ['Cannot renamed parameter ''',old_name,''' to ''',new_name,''' because the former does not exist.'])
      assert(~obj.isparameter(new_name),...
        ['Cannot renamed parameter ''',old_name,''' to ''',new_name,''' because the latter already exists.'])
      S_now=obj.get(old_name);
      obj=obj.delete({old_name});
      S_now.name=new_name;
      obj=obj.set(S_now);
    end
    function obj=rename_silent(obj,old_name,new_name)
      if obj.isparameter(old_name) && ~obj.isparameter(new_name)
        obj=obj.rename(old_name,new_name);
      end
    end
    %% multi-object methods
    function out=dup(in)
      out=varargs({});
      out.S=in.S;
      out.update_dynamic_fields;
    end
    %obj is updated with the common entries in obj_new. New entries in obj_new are ignored
    function obj=merge(obj,obj_new)
      if ~isa(obj_new,'varargs');obj_new=varargs(obj_new).delete(varargs.reserved_fields{:});end
      for i=1:obj_new.length
        name=obj_new.S(i).name;
        if obj.isparameter(name)
          obj.set(obj_new.get(i));
        end
      end
    end
    %obj receives the new entries in obj_new. Common entries are ignored.
    function obj=append(obj,obj_new)
      if ~isa(obj_new,'varargs');obj_new=varargs(obj_new).delete(varargs.reserved_fields{:});end
      for i=1:obj_new.length
        name=obj_new.S(i).name;
        if ~obj.isparameter(name)
          obj.set(obj_new.get(i));
        end
      end
    end
    %obj receives all entries in obj_new. Nothing is ignored.
    function obj=join(obj,obj_new)
      if ~isa(obj_new,'varargs');obj_new=varargs(obj_new).delete(varargs.reserved_fields{:});end
      for i=1:obj_new.length
        obj.set(obj_new.get(i));
      end
    end
    %checks if two objs are equal, overloads with ==
    function [out,msg]=eq(obj1,obj2)
      %defaults
      out=false;
      %check equality
      if ~cells.isequal(obj1.template_fields,obj2.template_fields)
        msg='template_fields differ';
        return
      end
      if ~cells.isequal(obj1.reserved_fields,obj2.reserved_fields)
        msg='reserved_fields differ';
        return
      end
      if ~numel(obj1.S)==numel(obj2.S)
        msg='obj length differs';
        return
      end
      if any(arrayfun(@(i,j) ~strcmp(i.name,j.name),obj1.S,obj2.S))
        msg='obj field names differs';
        return
      end
      if any(arrayfun(@(i,j) ~cells.isequal(i.value,j.value),obj1.S,obj2.S))
        msg='obj field values differs';
        return
      end
      if any(arrayfun(@(i,j) ~strcmp(func2str(i.validation),func2str(j.validation)),obj1.S,obj2.S))
        msg='obj field validation differs';
        return
      end
      %equality verified
      out=true;
      msg='';
    end
    %% parser
    function [p,obj]=declare(obj,p)
      p.PartialMatching=false;
      for i=1:obj.length
        p.addParameter(obj.S(i).name,obj.S(i).value,obj.S(i).validation)
      end
    end
    %% propagate to other object
    function o=save(obj,o,skip_list)
      if ~exist('skip_list','var')
        skip_list={};
      end
      for i=1:obj.length
        name=obj.S(i).name;
        if cells.isincluded(skip_list,name); continue; end
        if any(strcmp(properties(o),name)) || isfield(o,name) % isprop(o,n) doesn't work
          v=obj.S(i).value;
          try
            if isscalar(v)
              o.(name)=v;
            else
              %This is useful to handle certain class of values in a robust way
              switch class(v)
              case 'char'
                o.(name)=transpose(v(:));
              otherwise
                o.(name)=v;
              end
            end
          catch e
            warning(['Could not save parameter ''',name,''' because: ',e.message])
          end
        end
      end
    end
  end
end