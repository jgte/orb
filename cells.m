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
    function out=iscellofcells(in,depth,depth_now)
      if ~exist('depth_now','var')
        depth_now=0;
      end
      out=true;
      if ~iscell(in)
        out=false;
        return
      end
      if depth_now<depth
        for i=1:numel(in)
          if ~cells.iscellofcells(in{i},depth,depth_now+1);
            out=false;
            return
          end
        end
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
    function out=iscellstrfind(cellstrin,strin)
      if iscellstr(strin) && ischar(cellstrin)
        tmp=strin;
        strin=cellstrin;
        cellstrin=tmp;
      end
      out=~cellfun(@isempty,strfind(cellstrin,strin));
    end
  end
end