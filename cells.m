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
      %http://stackoverflow.com/questions/3231580/matlab-comparison-of-cell-arrays-of-string
      out=isempty(setxor(c1,c2));
    end
    function out=allequal(in)
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
          if ~cells.iscellofcells(in{i},depth,depth_now+1);
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
        out=isempty(i);
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
    %% overloading strfind
    function out=isstrfind(cellstrin,strin) %this used to be called iscellstrempty
      if iscellstr(strin) && ischar(cellstrin)
        %switch it around
        tmp=strin;
        strin=cellstrin;
        cellstrin=tmp;
      end
      if isempty(strin)
        out=cellfun(@isempty,cellstrin);
      else
        out=~cellfun(@isempty,strfind(cellstrin,strin));
      end
    end
    function out=strfind(cellstrin,strin) %this used to be called cellstrfind
      out=find(cells.isstrfind(cellstrin,strin));
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
      out=any(cells.isstrfind(cellstrin,strin));
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
            out=transpose(strsplit(in,char(10)));
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
    function [out,in]=isnum(in)
      if isnumeric(in);                   out=true; in=out; return; end
      if   ~iscell(in) && ~iscellstr(in); out=false; in=[]; return; end
      for i=1:numel(in)
        try 
          in{i}=str.num(in{i});
        catch
          out=false;in=[];return
        end
      end
      out=true;
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
    %depending on the value of 'direction' (defaults to 'set'):
    % - if get: checks if there is only one cell entry, if so return it; otherwise nothing changes
    % - if set: checks if it's a cell array, if so return it; otherwise make it a cell and return in
    function io=scalar(io,direction)
      if ~exist('direction','var') || isempty(direction)
        direction='get';
      end
      switch lower(direction)
      case 'set'; if ~iscell(io);                 io={io};  end
      case 'get'; if  iscell(io) && isscalar(io); io=io{1}; end
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
    %wrapper for vararginclean and vararginget
    function out=varargin(mode,in,parameters)
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