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
      out=cellfun(@isempty,in);
    end
    function out=rm_empty(in)
      out=in(~cells.isempty(in));
    end
    function out=iscellstrempty(cellstrin,strin)
      if iscellstr(strin) && ischar(cellstrin)
        %switch it around
        tmp=strin;
        strin=cellstrin;
        cellstrin=tmp;
      end
      out=~cellfun(@isempty,strfind(cellstrin,strin));
    end
    function out=cellstrfind(cellstrin,strin)
      out=find(cells.iscellstrempty(cellstrin,strin));
    end
    %returns a cell array with the sizes of the contents
    function out=size(in)
      out=cellfun(@(i) size(i),in,'UniformOutput',false);
    end
    %this is similar to cell2mat but any type of objects is supported
    function out=c2m(in)
      %trivial call
      if ~iscell(in) || iscellstr(in)
        out=in;
        return
      end
      assert(ndims(in)<=2,'Can only handle matrices or vectors')
      assert(cells.allequal(cells.size(in)),'Need all entries of ''in'' to have the same size')
      out=[];
      for i=1:size(in,1)
        out_row=[];
        for j=1:size(in,2)
          out_row=[out_row,in{i,j}]; %#ok<AGROW>
        end
        out=[out;out_row]; %#ok<AGROW>
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
    %checks if cell array 'in' has size equal to 's':
    % - if not and scalar, create cell array of size 's' the value content of 'in';
    % - if not, issue error;
    % - if so, do nothing;
    function out=deal(in,s)
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
    %returns a cell array with the contents of 'fieldname' of all entries of the cell array of structures S
    function out=deal_struct(S,fieldname)
      out=cell(size(S));
      for i=1:numel(S)
        out{i}=S{i}.(fieldname);
      end
    end
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
  end
end