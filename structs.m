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
    %creates a string with the contents of the deep structure
    function out=str(S,field_list,varname)
      if ~exist('field_list','var')
        field_list=structs.field_list(S);
      end
      if ~exist('varname','var')
        varname='struct';
      end
      out=cell(size(field_list));
      for i=1:numel(field_list)
        v=structs.get_value(S,field_list{i});
        out{i}=[varname,'.(',strjoin(field_list{i},').('),')=',str.show(v),', class=',class(v)];
      end
      out=strjoin(out,char(10));
    end
    %'field_path' is cell array with the sub-field path to the value to be retrieved from structure in 'in'
    function out=get_value(S,field_path)
      if ~exist('field_path','var') || isempty(field_path)
        out=S;
        return
      end
      %check if this field exists
      if ~isfield(S,field_path{1})
        out=[];
      else
        out=structs.get_value(S.(field_path{1}),field_path(2:end));
      end
    end
    %'field_path' is cell array with the sub-field path to the value to be set in structure in 'in'
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
      %loop over all fields
      for i=1:numel(fl)
        %get the object
        o=structs.get_value(S,fl{i});
        %apply the method
        o=o.(method)(varargin{:});
        %save the operated object back in the structure
        S=structs.set_value(S,fl{i},o);
      end
    end
    %applies 'method' to the object at the leafs of structure 'Sout', considering the corresponding leafs of 'Sin'
    function Sout=objmethod2(method,Sout,Sin,varargin)
      %get field list
      fl=structs.field_list(Sout);
      %loop over all fields
      for i=1:numel(fl)
        %get the input object
        Oin=structs.get_value(Sin,fl{i});
        %get the output object
        Oout=structs.get_value(Sout,fl{i});
        %apply the method
        Oout=Oout.(method)(Oin,varargin{:});
        %save the operated object back in the structure
        Sout=structs.set_value(Sout,fl{i},Oout);
      end
    end
    function out=isleaf(S,field_path,non_empty)
      if ~exist('field_path','var')
        field_path={};
      end
      if ~exist('non_empty','var')||isempty(non_empty)
        non_empty=false;
      end
      val=structs.get_value(S,field_path);
      out=~isstruct(val);
      if non_empty
        out=out && ~isempty(val);
      end
    end
  end
end