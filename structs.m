classdef structs
  methods(Static)
    function test
      a=struct( ...
        'A', 1, ...
        'B','2',...
        'C',struct( ...
          'AA',{{10,11}},...
          'BB', '20 21', ...
          'CC',struct(...
            'AAA',[100,101],...
            'BBB',{{200,201}},...
            'CCC','300 301'...
          )...
        )...
      );
      field_path={...
        {'A'},...
        {'C','BB'},...
        {'C','CC'}...
      };
      disp('---- structs.get_value ----')
      disp(structs.str(a,field_path))
      disp('---- structs.set_value ----')
      a=structs.set_value(a,field_path{2},-1);
      disp(structs.str(a,field_path))
      disp('---- structs.field_list ----')
      disp(structs.str(a))
    end
    %% utils
    %creates a string with the contents of the deep structure (if 'varname' is '_', the output is appropriate for a filename)
    function out=str(S,field_list,varname,show_class)
      if ~exist('field_list','var') || isempty(field_list)
        field_list=structs.field_list(S);
      end
      if ~exist('varname','var') || isempty(varname)
        varname='struct';
      end
      if ~exist('show_class','var') || isempty(show_class)
        show_class=true;
      end
      out=cell(size(field_list));
      for i=1:numel(field_list)
        v=structs.get_value(S,field_list{i});
        if strcmp(varname,'_')
          out{i}=[strjoin(field_list{i},','),str.show(v,'',varname)];
        else
          out{i}=[varname,'.',strjoin(field_list{i},'.'),'=',str.show(v)];
          if show_class
            out{i}=[out{i},', class=',class(v)];
          end
        end
      end
      if strcmp(varname,'_')
        out=strjoin(out,'_');
      else
        out=strjoin(out,char(10));
      end
    end
    %% get/set values
    %'field_path' is cell array with the sub-field path to the value to be retrieved from structure in 'in' (i.e. a cell of cells)
    function out=get_value(S,field_path,search_flag)
      if ~exist('search_flag','var')||isempty(search_flag)
        %do not search by default, otherwise isleaf fails (possibly others)
        search_flag=false;
      end
      %vector mode
      %NOTICE: this illustrates that (the scalar mode of) this routine is formally inconsistent with the idea that
      %        field_path is a cells of cells (the vector mode was added after defining this method)
      %TODO: needs revision
      if cells.iscellofcells(field_path,1)
        out=cellfun(@(i) structs.get_value(S,i,search_flag),field_path,'UniformOutput',false);
        return
      end
      if ~exist('field_path','var') || isempty(field_path)
        out=S;
        return
      end
      %check if this field exists
      if ~isfield(S,field_path{1}) && ~isprop(S,field_path{1})
        if search_flag
          %climb down the field_path list, maybe a valid field appears later
          out=structs.get_value(S,field_path(2:end),true);
        else
           out=[];
        end
      else
        out=structs.get_value(S.(field_path{1}),field_path(2:end),search_flag);
      end
    end
    %'field_path' is cell array with the sub-field path to the value to be set in structure in 'in' (i.e. a cell of cells)
    function S=set_value(S,field_path,value)
      if ~exist('field_path','var')
        field_path={};
      end
      %trivial call
      if isempty(field_path)
        S=value;
        return
      end
      %initiate this field, if needed
      if ~isfield(S,field_path)
        S(1).(field_path{1})={};
      end
      if numel(field_path)==1
        S.(field_path{1})=value;
      else
        S.(field_path{1})=structs.set_value(S.(field_path{1}),field_path(2:end),value);
      end
    end
    %like get_value but for all structure entries
    function out=get_value_all(S)
      out=structs.get_value(S,structs.field_list(S));
    end
    %% field list
    %creates a cell array of "field_path"s, each being a cell array with field names, to be used with structs.get/set_value
    function out=field_list(S,parents)
      if ~exist('parents','var')
        parents={};
      end
      %check if this is a structure
      if ~isstruct(S)
        %if not, then this is a leaf field (no structures inside this field)
        leaf=true;
      else
        %try to get the field names of the structure in this field
        try
          fn=fieldnames(S);
          %this is not a leaf name
          leaf=false;
        catch
          %failed to get names of fields, so this is a leaf name
          leaf=true;
        end
      end
      %check if this is a leaf field
      if leaf
        %if so, return the parents (if any)
        if isempty(parents)
          out={};
        else
          out={parents};
        end
      else
        %if not, then loop over all field names and make a recursive call
        out=cell(0);c=0;
        for i=1:numel(fn)
          sfn=structs.field_list(S.(fn{i}),[parents,fn(i)]);
          for j=1:numel(sfn)
            c=c+1;out{c}=sfn{j};
          end
        end
      end
    end
    %returns true if field lists are the same in both structure
    function [out,fl1,fl2]=iseq_field_list(S1,S2,parents1,parents2)
      if ~exist('parents1','var');parents1={};end
      if ~exist('parents2','var');parents2={};end
      fl1=structs.field_list(S1,parents1);
      fl2=structs.field_list(S2,parents2);
      fl1f=cells.flatten(fl1);
      fl2f=cells.flatten(fl2);
      out=numel(fl1f) == numel(fl2f) && all(strcmp(fl1f,fl2f));
    end
    %like field_list but accepts '*' entries in field_path
    function out=field_list_glob(S,field_path)
      %sanity on type
      assert(iscellstr(field_path),['input ''field_path'' must be a cell of strings, not a ',class(field_path),'.'])
      %get the top-most list of datanames
      top_list=structs.field_list(S);
      %init outputs
      out=cell(0);c=0;
      %loop over all list entries
      for i=1:numel(top_list)
        %assume this field path is to keep
        keep=true;
        %loop over the depth of this field path entry
        for j=1:numel(top_list{i})
          %check if this is a globbed field
          if strcmp(field_path{j},'*')
            %continue searching
            continue
          elseif ~strcmp(top_list{i}{j},field_path{j})
            %discard it
            keep=false;
            break
          end
        end
        %check if this field path is to be saved
        if keep
          c=c+1;out{c}=top_list{i};
        end
      end
    end
    %% member class interfaces
    %checks if a method/field/property exists for a general object
    function out=respondto(S,method)
      out=...
                    isfield(S, method) || ...
      any(strcmp(   methods(S),method))|| ...
      any(strcmp(properties(S),method));
    end
    %applies ''method'' to the object at the leafs of structure ''S''
    function S=objmethod(method,S,varargin)
      %get field list
      fl=structs.field_list(S);
      %maybe this is not a structure
      if isempty(fl)
        S=S.(method);
      else
        %loop over all fields
        for i=1:numel(fl)
          %get the object
          o=structs.get_value(S,fl{i});
          if ~isempty(o)
            %apply the method
            o=o.(method)(varargin{:});
          end
          %save the operated object back in the structure
          S=structs.set_value(S,fl{i},o);
        end
      end
    end
    %applies 'method' to the object at the leafs of structure 'Sout', considering the corresponding leafs of 'Sin'
    %NOTICE: empty entries in Sout are assigned with the corresponding entry in Sin, i.e. the 'method' implicitly becomes '='
    function Sout=objmethod2(method,Sout,Sin,varargin)
      %get field list
      fl=structs.field_list(Sout);
      %loop over all fields
      for i=1:numel(fl)
        %get the input object
        Oin=structs.get_value(Sin,fl{i});
        %get the output object
        Oout=structs.get_value(Sout,fl{i});
        %apply the method, if Oout is not empty
        if ~isempty(Oout)
          Oout=Oout.(method)(Oin,varargin{:});
        else
          Oout=structs.get_value(Sin,fl{i});
        end
        %save the operated object back in the structure
        Sout=structs.set_value(Sout,fl{i},Oout);
      end
    end
    %% utils
    function out=isleaf(S,field_path,non_empty)
      if ~exist('field_path','var')
        field_path={};
      end
      if ~exist('non_empty','var')||isempty(non_empty)
        non_empty=false;
      end
      val=structs.get_value(S,field_path,false);
      out=~isstruct(val);
      if non_empty
        out=out && ~isempty(val);
      end
    end
    %returns a new structure with all fields (and corresponding values) that either:
    % - contain the substring 'fieldname_part' or
    % - are present in the cellstr 'fieldname_part'
    function out=filter(S,fieldname_part)
      fn=fieldnames(S);
      if numel(fn)==0
        out=struct([]);
      else
        if iscellstr(fieldname_part)
          out=structs.copy(S,struct([]),fieldname_part,true);
        elseif ischar(fieldname_part)
          for i=1:numel(fn)
            if strfind(fn{i},fieldname_part)
              out.(fn{i})=S.(fn{i});
            end
          end
        else
          error(['Cannot handle input ''fieldname_part'' of class ''',class(fieldname_part),'''.'])
        end
      end      
    end
    %renames field 'field_old' to 'field_new' of structure 'S'
    function S=rename(S,field_old,field_new)
      if iscell(field_old) && iscell(field_new) && numel(field_old)==numel(field_new)
        for i=1:numel(field_old)
          S=structs.rename(S,field_old{i},field_new{i});
        end
        return
      end
      if strcmp(field_old,field_new) || ~isfield(S,field_old)
        return
      end
      S.(field_new)=S.(field_old);
      S=rmfield(S,field_old);
    end
    %returns a new structure with all fields (and corresponding values) without 'fieldname_part'
    function out=fieldname_strip(S,fieldname_part)
      %get all fieldnames
      fn_old=fieldnames(S);
      %remove field_part from field names
      fn_new=strrep(fn_old,fieldname_part,'');
      %rename structure
      out=structs.rename(S,fn_old,fn_new);
    end
    %applies function 'f' to all fields (top level) of structure 'S' (this is NOT the same as matlab's structfun)
    function S=fun(f,S)
      fn=fieldnames(S);
      for i=1:numel(fn)
        S.(fn{i})=f(S.(fn{i}));
      end
    end
    %if S does not have 'fieldname', adds 'field' to it; if S has 'fieldname', use 'method' on existing 'field'
    function S=build(S,fieldname,field,method)
      %check for vector mode
      if iscellstr(fieldname) && iscell(field)
        str.sizetrap(fieldname,field)
        for i=1:numel(fieldname)
          S=structs.build(S,fieldname{i},field{i},method);
        end
      else
        if ~isempty(field)
          if isempty(S) || ~isfield(S,fieldname)
            S.(fieldname)=field;
          else
            S.(fieldname)=structs.objmethod2(method,S.(fieldname),field);
          end
        end
      end
    end
    %% multi-object manipulation
    function Sout=copy(Sin,Sout,field_path,copy_empty_values)
      if ~exist('field_path','var') || isempty(field_path)
        field_path=fieldnames(Sin);
      end
      if ~exist('copy_empty_values','var') || isempty(copy_empty_values)
        copy_empty_values=false;
      end
      for i=1:numel(field_path)
        if ischar(field_path{i})
          %translate simple vectors of fieldnames to field paths
          fp=field_path(i);
        else
          %you're on your own, better be a cell of cells
          fp=field_path{i};
        end
        value=structs.get_value(Sin,fp,false);
        if copy_empty_values || ~isempty(value)
          Sout=structs.set_value(Sout,fp,value);
        end
      end
    end
  end
end