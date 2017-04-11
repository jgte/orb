classdef plotting
  methods(Static)
    function out=line_handles(axis_handle)
      tmp=get(axis_handle,'Children');
      c=0;out=[];
      for i=1:numel(tmp)
        if isprop(tmp(i),'LineWidth')
          c=c+1;
          out(c)=tmp(i); %#ok<AGROW>
        end
      end
    end
    function v=common_axis_limits(axis_handle)
      %set output
      v=[inf -inf inf -inf];
      %get line handles
      lh=plotting.line_handles(axis_handle);
      %loop over all lines
      for i=1:numel(lh)
        %get x-data extremeties
        mm=minmax(get(lh(i),'XData'));
        %accumulate inclusive x-axis limits
        v(1)=min([v(1),mm(1)]);
        v(2)=max([v(2),mm(2)]);
        %get y-data extremeties
        mm=minmax(get(lh(i),'YData'));
        %accumulate inclusive x-axis limits
        v(3)=min([v(3),mm(1)]);
        v(4)=max([v(4),mm(2)]);
      end
    end
  end
end