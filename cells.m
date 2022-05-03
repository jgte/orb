classdef cells
    methods(Static)
    function y = flatten(x)
      % https://github.com/ronw/ronw-matlab-tools/blob/master/celltools/flatten.m
      if ~iscell(x)
        y = {x};
      else
        y = {};
        for n = 1:length(x)
          tmp = cells.flatten(x{n});
          y = [y(:);tmp(:)];
        end
      end
    end
    function out=isequal(c1,c2)
      %ensures all correponding arbitrary class entries are the same
      out=false;
      if ~iscell(c1)
        out=cells.isequal(c2,{c1});
        return
      end
      if numel(c1)~=numel(c2); return; end
      for i=1:numel(c1)
        if ~isa(c1{i},class(c2{i}))
          return
        end
        if isnumeric(c1{i}) || islogical(c1{i})
          if ~isequaln(c1{i},c2{i})
            return
          else
            continue
          end
        end
        if ischar(c1{i})
          if ~strcmp(c1{i},c2{i})
            return
          else
            continue
          end
        end
        if iscell(c1{i})
          if ~cells.isequal(c1{i},c2{i})
            return
          else
            continue
          end
        end
        if isa(c1{i},'function_handle')
          if ~strcmp(func2str(c1{i}),func2str(c2{i}))
            return
          else
            continue
          end
        end
        if isdatetime(c1{i})
          if c1{i} ~= c2{i}
            return
          else
            continue
          end
        end
        try
          %TODO: replace obj.isequal with obj.eq and use the overloaded == operator here
          if ~c1{i}.isequal(c2{i})
            return
          else
            continue
          end
        catch ME
          error(['Cannot handle data of class ',class(c1{i}),'.'])
        end
      end
      out=true;
    end
    function out=isequalstr(c1,c2)
      %ensures all correponding string entries are the same
      %http://stackoverflow.com/questions/3231580/matlab-comparison-of-cell-arrays-of-string
      out=isempty(setxor(c1,c2));
    end
    function out=allequal(in)
      %ensures all entries are exactly the same
      in=cells.flatten(in);
      for i=2:numel(in)
        if any(in{1} ~= in{i})
          out=false;
          return
        end
      end
      out=true;
    end
    function out=iscellofcells(in,depth,depth_now)
      if ~exist('depth_now','var')
        depth_now=0;
      end
      if ~iscell(in)
        out=false;
        return
      end
      if depth_now<depth %&& ~isempty(in) %an empty cell returns true with depth 1
        for i=1:numel(in)
          if ~cells.iscellofcells(in{i},depth,depth_now+1)
            out=false;
            return
          end
        end
        out=~isempty(in);
      else
        out=iscell(in);
      end
    end
    function out=isempty(in)
      if iscell(in)
        out=cellfun(@isempty,in);
        if isempty(out); out=true;end
      else
        out=isempty(in);
      end
    end
    function out=rm_empty(in)
      out=in(~cells.isempty(in));
    end
    function out=rm_duplicates(in)
      %ignore empty entries
      out=cells.rm_empty(in);
      %repeat the search often enough
      for c=1:numel(out)-1
        %loop over all elements of the cell array starting form the second
        for i=2:numel(out)
          %adapt to smaller out
          if numel(out)<i;break;end
          %if this is the same as the first, delete it
          if out{1}==out{i}; out(i)=[]; end
        end
        %we're done if there's only one element
        if numel(out)==1;break;end
      end
    end
    function out=ismethod(cellstrin,method)
     out=cellfun(@(i) any(strcmp(   methods(i),method)),cellstrin);
    end
    function out=isprop(cellstrin,method)
      out=cellfun(@(i) isfield(i, method) || ...
         any(strcmp(properties(i),method)),cellstrin);
    end
    %checks if a method/field/property exists for a general object
    function out=respondto(cellstrin,method)
      out=cells.ismethod(cellstrin,method) | cells.isprop(cellstrin,method);
    end
    function out=allrespondto(cellstrin,method)
      out=cellfun(@(i) all(cells.respondto(cellstrin,i)),cells.scalar(method,'set'));
    end
    %% overloading strfind
    function io=cellstr(io)
      io=io(cellfun(@ischar,io));
    end
    function out=iscellstr(in)
      out=iscell(in) && any(cellfun(@ischar,in));
    end
    function out=isstrfind(cellstrin,strin) %this used to be called iscellstrempty
      if cells.iscellstr(strin) && ischar(cellstrin)
        %switch it around
        tmp=strin;
        strin=cellstrin;
        cellstrin=tmp;
      elseif cells.iscellstr(strin) && cells.iscellstr(cellstrin)
        %expand
        out=cellfun(@(i) cells.isstrfind(cellstrin,i),strin,'UniformOutput',false);
        return
      end
      if isempty(strin)
        out=cellfun(@isempty,cellstrin);
      elseif isempty(cellstrin)
        out={};
      else
        %pick the elements of cellstrin that are chars
        cellstrwork=cellstrin(cellfun(@ischar,cellstrin));
        %check if strin is in cellstrin
        out=contains(cellstrwork,strin);
        %old implementation:
        %out=~cellfun(@isempty,strfind(cellstrin(cellfun(@ischar,cellstrin)),strin));
      end
    end
    function out=strfind(cellstrin,strin) %this used to be called cellstrfind
      if isempty(cellstrin)
        out={};
      else
        out=find(cells.isstrfind(cellstrin,strin));
      end
    end
    function out=isstrequal(cellstrin,strin) %this is similar to isstrfind but the strings have to be an exact match
      if cells.iscellstr(strin) && ischar(cellstrin)
        %switch it around
        tmp=strin;
        strin=cellstrin;
        cellstrin=tmp;
      end
      if isempty(strin)
        out=cellfun(@isempty,cellstrin);
      elseif isempty(cellstrin)
        out={};
      else
        out=cellfun(@(i) all(size(i)==size(strin)) && all(i==strin),cellstrin);
      end
    end
    function out=strequal(cellstrin,strin)
      if isempty(cellstrin)
        out={};
      else
        out=find(cells.isstrequal(cellstrin,strin));
      end
    end
    function out=strcount(cellstrin,strin) 
      if isempty(cellstrin)
        out=0;
      else
        out=sum(cells.isstrfind(cellstrin,strin));
      end
    end
    function io=rm_strfind(io,strin)
      idx=cells.strfind(io,strin);
      if isempty(idx); return; end
      for i=1:numel(idx)
        io{idx(i)}=[];
      end
      io=cells.rm_empty(io);
    end    
    function out=isincluded(cellstrin,strin)
      %N.B. this is not the same as ismember (isincluded returns true of strin is a sub-string
      %of any entries in cellstrin, unlike ismember where the whole entry must be the same)
      out=cells.c2m(cells.isstrfind(cellstrin,strin));
      out=~isempty(out) && any(out(:));
    end
    function out=cellstrget(cellstrin,strin)
      out=cellstrin(cells.isstrfind(cellstrin,strin));
    end
    %% overloading strcmp
    function out=isstrcmp(cellstrin,strin) %this is mainly to remind me the operational difference between strcmp and strfind
      if iscellstr(strin) && ischar(cellstrin)
        %switch it around
        tmp=strin;
        strin=cellstrin;
        cellstrin=tmp;
      end
      out=strcmp(cellstrin,strin);
    end
    function out=strcmp(cellstrin,strin)
      out=find(cells.isstrcmp(cellstrin,strin));
    end
    %% overloading cell2mat
    %this is similar to cell2mat but any type of objects is supported
    %NOTICE: cell strings are passed through unchanged, use cells.num to convert char to num
    function out=c2m(in)
      %trivial call
      if ~iscell(in) || iscellstr(in)
        out=in;
        return
      end
      assert(ndims(in)<=2,'Can only handle matrices or vectors')
      if cells.allequal(cells.size(in))
        out=[];
        for i=1:size(in,1)
          out_row=[];
          for j=1:size(in,2)
            if ischar(in{i,j})
              try %#ok<TRYNC>
                in{i,j}=str.num(in{i,j});
              end
            end
            out_row=[out_row,in{i,j}]; %#ok<AGROW>
          end
          out=[out;out_row]; %#ok<AGROW>
        end
      else
        %can't do much here, because the cells do not have the same size
        out=in;
      end
    end
    %this is similar to num2cell but any type of objects is supported (and there's special handling for strings)
    function out=m2c(in)
      if iscell(in)
        out=in;
      else
        if ischar(in)
          if size(in,1)==1
            out=transpose(strsplit(in,newline));
          else
            out=cell([size(in,1),1]);
            for i=1:size(in,1)
              out{i}=strtrim(in(i,:));
            end
          end
        else
          out=cell(size(in));
          for i=1:numel(out)
            out{i}=in(i);
          end
        end
      end
    end
    %% overloading misc stuff
    %returns a cell array with the sizes of the contents
    function out=size(in)
      out=cellfun(@(i) size(i),in,'UniformOutput',false);
    end
    %% handling numeric cells
    %returns true if all entries in 'in' are numeric of strings that can be converted to numeric
    %NOTICE: always returns a cell array
    function [out,in]=isnum(in)
      if isnumeric(in)
        out=true;
        in=num2cell(in);
      elseif ischar(in)
        [out,in]=cells.isnum({in});
      elseif iscell(in)
        for i=1:numel(in)
          try 
            in{i}=str.num(in{i});
          catch
            out=false;in=[];return
          end
        end
        out=true;
      else
        out=false;
        in=[];          
      end
    end
    function in=num(in)
      [~,in]=cells.isnum(in);
    end
    %% misc stuff
    %checks if cell array 'in' has size equal to 's':
    % - if so, do nothing;
    % - if not and scalar, create cell array of size 's' the value content of 'in';
    % - if not, issue error;
    function out=deal(in,s)
      %sanity
      assert(isnumeric(s),'Expecting input ''s'' to be numeric, representing the size of ''out''.')
      %enforce cell class
      out=cells.m2c(in);
      if any(size(out)~=s)
        tmp=cell(s);
        if isscalar(out)
          tmp(:)=out;
        elseif isempty(out)
          %do nothing, already have empty cells with the right size
        else
          assert(numel(tmp)==numel(out),...
            ['Cannot deal array ''in'' with size [',num2str(size(in)),...
            '] into an array of size [',num2str(s),'].'])
          tmp(:)=out(:);
        end
        out=tmp;
      end
    end
    function [io,changed]=scalar(io,direction)
      %[io,changed]=scalar(io,direction)
      %depending on the value of 'direction' (defaults to 'get'):
      % - if get: checks if there is only one cell entry, if so return it (changed is true);
      % otherwise nothing changes (changed is false)
      % - if set: checks if it's a cell array, if so return it (changed is
      % false); otherwise make it a cell and return io (changed is true)
      if ~exist('direction','var') || isempty(direction)
        direction='get';
      end
      changed=false;
      switch lower(direction)
      case 'set'; if ~iscell(io);                 io={io};  changed=true; end
      case 'get'; if  iscell(io) && isscalar(io); io=io{1}; changed=true; end
      otherwise; error(['Cannot handle input ''direction'' with value ''',direction,'''.'])
      end
    end
    %returns a cell array with the contents of 'fieldname' of all entries of the cell array of structures S
    function out=deal_struct(S,fieldname)
      out=cell(size(S));
      for i=1:numel(S)
        out{i}=S{i}.(fieldname);
      end
    end
    function io=patch_empty(io,patch)
      if any(isempty(io));io(isempty(io))=patch(isempty(io));end
    end
    %% stuff to handle varargin
    %removes from the varargin-like cell array 'in' the parameters with names defined in the cell array 'parameters'.
    %e.g.: varargin{{'a',1,'b','2','c',{3}},{'a','c'}) => {'b','2'}
    function out=vararginclean(in,parameters)
      %trivial call
      if ~exist('parameters','var')
        out=in;
        return
      end
      if ischar(parameters)
        parameters={parameters};
      end
      out=in;
      for i=1:numel(parameters)
        idx=0;
        for j=1:2:numel(out)
          if strcmp(out{j},parameters{i})
            idx=j;
            break
          end
        end
        if idx>0
          out=out([1:j-1,j+2:end]);
        end
      end
    end
    %retrieves from the varargin-like cell array 'in' the parameters with names defined in the cell array 'parameters'.
    %e.g.: varargin{{'a',1,'b','2','c',{3}},{'a','c'}) => {'a',1,'c',{3}}
    function out=vararginget(in,parameters)
      %trivial call
      if ~exist('parameters','var')
        out=in;
        return
      end
      if ischar(parameters)
        parameters={parameters};
      end
      out=cell(1,numel(parameters)*2);
      for i=1:numel(parameters)
        idx=0;
        for j=1:2:numel(in)
          if strcmp(in{j},parameters{i})
            idx=j;
            break
          end
        end
        if idx>0
          out(i*2-1:i*2)=in(idx:idx+1);
        end
      end
      out=cells.rm_empty(out);
    end
    %wrapper for vararginclean and vararginget (if parameters are empty, get returns nothing and clean returns the same)
    function out=varargin(mode,in,parameters)
      %trivival call
      if ~exist('parameters','var')
        out=in;
        return
      end
      switch mode
      case 'get';   out=cells.vararginget(  in,parameters);
      case 'clean'; out=cells.vararginclean(in,parameters);
      otherwise;    error(['Cannot handle mode ''',mode,'''.'])
      end  
    end
    %% first/last/nth wrapper
    function out=ith(in,i)
      out=in{i};
    end
    function out=first(in)
      out=cells.ith(in,1);
    end
    function out=last(in)
      out=cells.ith(in,numel(in));
    end
  end
end