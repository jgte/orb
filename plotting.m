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
    function line_color(mode,axis_handle)
      % PLOT_LINE_COLOR assigns different colors to all the plotted lines.
      %
      %   By default PLOT_LINE_COLOR works on the current graphics and uses a
      %   colormap that shows colors that are also different inr black and white
      %   visualizations.
      %
      %   PLOT_LINE_COLOR(CLR_STR) sets the color of each line to the colors
      %   specified in the cell array containing single characters CLR_STR.
      %   Recognized characters are the ones also recognized by the plot
      %   function. See 'help plot' for more information.
      %
      %   PLOT_LINE_COLOR('reversed') uses the sames colors as default but in
      %   reverse order.
      %
      %   PLOT_LINE_COLOR('multi') uses the colors specified at:
      %   http://blogs.mathworks.com/pick/2008/08/15/colors-for-your-multi-line-plots/#comment-12842
      %
      %   PLOT_LINE_COLOR(COLORMAP_NAME) uses the a COLORMAP_NAME to color the
      %   lines in the plot, see the documentation of the colormap function.
      %
      %   PLOT_LINE_COLOR(CLR_STR,H) as before, but now will work on the figures
      %   with handles H.

      % TODO: add copyright
    
      % handle inputs
      if ~exist('mode','var') || isempty(mode)
          mode='spiral';
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      % translate mode
      switch strrep(lower(mode),'-reversed','')
      case {'flat','boring'}
        branch=0;
      case{'multi','many'}
        branch=1;
      case{'reversed','spiral'}
        branch=2;
      case{'jet','hsv','hot','cool','spring','summer','autumn','winter','gray','bone','copper','pink','lines'}
        branch=3;
      otherwise
        error(['Unknown mode ''',mode,'''.'])
      end
      %loop over all axis
      for h = axis_handle(:)'
        %getting lines in this axis
        lines=findobj(h,'Type','line');
        %trivial call
        if isempty(lines)
          continue
        end
        %get colororder according to branch
        switch (branch)
        case{0}
          colororder=[
            0     0     1 %blue
            1     0     0 %red
            0     1     0 %green
            0     1     1 %cyan
            1     0     1 %magenta
            1     1     0 %yellow
            0     0     0 %black
          ];
        case{1}
          % http://blogs.mathworks.com/pick/2008/08/15/colors-for-your-multi-line-plots/#comment-12842
          colororder = [
              0.00  0.00  1.00
              0.00  0.50  0.00
              1.00  0.00  0.00
              0.00  0.75  0.75
              0.75  0.00  0.75
              0.75  0.75  0.00
              0.25  0.25  0.25
              0.75  0.25  0.25
              0.95  0.95  0.00
              0.25  0.25  0.75
              0.75  0.75  0.75
              0.00  1.00  0.00
              0.76  0.57  0.17
              0.54  0.63  0.22
              0.34  0.57  0.92
              1.00  0.10  0.60
              0.88  0.75  0.73
              0.10  0.49  0.47
              0.66  0.34  0.65
              0.99  0.41  0.23
          ];
        case{2}
          % color spiral map size
          n=length(lines);
          %number of points in the extremeties to remove
          n_extrem=ceil(n*0.3);
          % http://bsp.pdx.edu/Publications/2006/SPM_McNames.pdf
          nc=n+2*n_extrem;
          np=2.5;
          % algorithm
          wn = sqrt(3/8)*[0;triang(nc-2);0];  % Triangularwindow function
          a12 = asin(1/sqrt(3));              % First rotation angle
          a23 = pi/4;                         % Second rotation angle
          t = linspace(sqrt(3),0,nc).';       % Independent variable
          r0 = t;                             % Initial red values
          g0 = wn.*cos(((t-sqrt(3)/2)*np*2*pi/sqrt(3))); % Initial green values
          b0 = wn.*sin(((t-sqrt(3)/2)*np*2*pi/sqrt(3))); % Initial blue values
          [ag,rd] = cart2pol(r0,g0);          % Convert RG to polar
          [r1,g1] = pol2cart(ag+a12,rd);      % Rotate& convert back
          b1 = b0;
          [ag,rd] = cart2pol(r1,b1);          % Convert RB to polar
          [r2,b2] = pol2cart(ag+a23,rd);      % Rotate & convert back
          g2 = g1;
          r = max(min(r2,1),0);               % Ensure finite precision
          g = max(min(g2,1),0);               % effects don't exceed
          b = max(min(b2,1),0);               % unit cube boundaries
          colororder = [r g b];               % Final colormap matrix
          %cropping extremeties (avoids too dark or too light colors)
          colororder([1:n_extrem,end-n_extrem+1:end],:)=[];
        case{3}
          % colormaps name
          switch mode
          %some colormaps include the white color, need to remove it
          case{'gray','bone','pink','hot'}
            colororder=eval([colormap_name,'(',num2str(length(lines)+1),')']);
            colororder=colororder(1:end-1,:);
          otherwise
            colororder=eval([colormap_name,'(',num2str(length(lines)),')']);
          end
        end
        %transform lines of colororder to cell array
        clr=cell(size(colororder,1),1);
        for i=1:size(colororder,1)
            clr{i}=colororder(i,:);
        end
      end
      %enforce reverse order of colours
      if ~isempty(strfind(mode,'-reversed'))
        clr=flipud(clr);
      end
      %matlab orders lines in the reverse way
      fix_idx=numel(lines):-1:1;
      %loop over all lines or clr entries
      for i=1:min([numel(lines),numel(clr)])
        set(lines(fix_idx(i)),'Color',clr{i})
      end
    end
    function line_width(w,axis_handle)
      % handle inputs
      if ~exist('w','var') || isempty(w)
          w=2;
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      %loop over all axis
      for ah = axis_handle(:)'
        lh=plotting.line_handles(ah);
        for i=1:numel(lh)		
          set(lh(i),'LineWidth',w)
        end
      end
    end
    function font_size(fs,axis_handle)
      % handle inputs
      if ~exist('fs','var') || isempty(fs)
          fs=18;
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      %loop over all axis
      for h = axis_handle(:)'
        set(    h,          'FontSize',fs);
        set(get(h,'Title' ),'FontSize',round(fs*1.3));
        set(get(h,'XLabel'),'FontSize',round(fs*1.1));
        set(get(h,'YLabel'),'FontSize',round(fs*1.2));
      end
    end
    function size(s,fig_handle)
      % handle inputs
      if ~exist('s','var') || isempty(s)
        s=200+[0,0,21,9]*50;
      end
      if ~exist('fig_handle','var') || isempty(fig_handle)
          fig_handle = gcf;
      end
      set(fig_handle, 'Position',      s,...
                      'PaperUnits',    'points',...
                      'PaperPosition', s);
    end
  end
end