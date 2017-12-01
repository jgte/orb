classdef plotting
  properties(Constant,GetAccess=private)
    default={...
        'plot_file_prefix',    '',        @(i) ischar(i);...
        'plot_file_suffix',    '',        @(i) ischar(i);...
        'plot_legend',         {},        @(i) iscellstr(i);...
        'plot_legend_suppress',{},        @(i) iscellstr(i);...
        'plot_legend_location','best',    @(i) ischar(i);...
     'plot_legend_fontname','FixedWidth', @(i) ischar(i);...
        'plot_ylabel',         '',        @(i) ischar(i);...
        'plot_xlabel',         '',        @(i) ischar(i);...
        'plot_xdate',          false,     @(i) islogical(i);...
        'plot_xdateformat',    '',        @(i) ischar(i);...
        'plot_xlimits',        [-inf,inf],@(i) isnumeric(i) && numel(i)==2;...
        'plot_ylimits',        [-inf,inf],@(i) isnumeric(i) && numel(i)==2;...
        'plot_size',    200+[0,0,21,9]*50,@(i) isnumeric(i) && numel(i)==4;...
        'plot_units',          'points',  @(i) ischar(i);...
        'plot_visible',        'on',      @(i) islogical(i);...
        'plot_fontsize_axis',   24,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_fontsize_title',  32,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_fontsize_label',  28,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_fontsize_legend', 28,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_title',           '',       @(i) ischar(i);...
        'plot_title_suppress',  {},       @(i) iscellstr(i);...
        'plot_title_suffix',    '',       @(i) ischar(i);...
        'plot_title_prefix',    '',       @(i) ischar(i);...
        'plot_grid',            true,     @(i) islogical(i);......
        'plot_line_width',      2,        @(i) isnumeric(i);...
        'plot_colormap',        '',       @(i) ischar(i) || ishandle(i);...
    };
  end
  methods(Static)
    %% primitives
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
    function out=markersize(n)
      out=round(50./log10(n.^2));
    end
    %% utility
    function fig_handle=figure(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default},varargin{:});
      %outputs
      fig_handle=figure('visible',str.logical(v.plot_visible,'onoff'));
      %plot size
      set(fig_handle, 'Position',          v.plot_size,...
                      'PaperUnits',        v.plot_units,...
                      'PaperPosition',     v.plot_size);
    end
    function out=enforce(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle',  gca,  @(i) ~isempty(i) && ishandle(i);...
      }},varargin{:});
    
      %outputs
      out.axis_handle=v.axis_handle;
            
      % enforce line properties
      line_handles=plotting.line_handles(out.axis_handle);
      for i=1:numel(line_handles)
        set(line_handles(i),'LineWidth',v.plot_line_width)
      end
      
      % start with current axis
      a=axis(out.axis_handle);
      %check if dates are requested
      if v.plot_xdate
        % enforce (possible) requested x-limits
        for i=1:2
          if isfinite(v.plot_xlimits(i))
            a(i)=datenum(v.plot_xlimits(i));
          end
        end
        % set auto x-label, unless one is explicity given
        if isempty(v.plot_xlabel)
          if ~strcmp(datestr(a(1),'yyyymmdd'),datestr(a(2),'yyyymmdd')) && ...
            (~strcmp(datestr(a(2),'HHMMSS'),'000000') || a(2)-a(1)>1)
            out.xlabel_handle=xlabel(out.axis_handle,...
              out.axis_handle,[datestr(a(1),'yyyy-mm-dd'),' to ',datestr(a(2),'yyyy-mm-dd')]...
            );
          else
            out.xlabel_handle=xlabel(out.axis_handle,out.axis_handle,datestr(a(1)));
          end
        end
        %enforce requested x date tick format
        if ~isempty(v.plot_xdateformat)
          datetick(out.axis_handle,'x',v.plot_xdateformat)
        end
      else
        % enforce (possible) requested x-limits
        for i=1:2
          if isfinite(v.plot_xlimits(i))
            a(i)=v.plot_xlimits(i);
          end
        end
      end
      % enforce requested y-limits
      for i=1:2
        if isfinite(v.plot_ylimits(i))
          a(i+2)=   v.plot_ylimits(i);
        end
      end
      % set axis limits (can be the ones matlab so wisely set)
      axis(out.axis_handle,a);
      
      %enforce labels
      if ~isempty(v.plot_xlabel)
        out.xlabel_handle=xlabel(out.axis_handle,v.plot_xlabel);
      elseif ~isfield(out,'xlabel_handle')
        out.xlabel_handle=[];
      end
      if ~isempty(v.plot_ylabel)
        out.ylabel_handle=ylabel(out.axis_handle,v.plot_ylabel);
      elseif ~isfield(out,'ylabel_handle')
        out.ylabel_handle=[];
      end
      
      %enforce title
      if ~isempty(v.plot_title)
        %suppress some parts, if requested
        title_str=setdiff(strsplit(v.plot_title,{' ','.'}),v.plot_title_suppress,'stable');
        %add prefix and suffix
        title_str=strjoin([{v.plot_title_prefix};title_str(:);{v.plot_title_suffix}],' ');
        %clean up the tile put it there
        out.title_handle=title(out.axis_handle,str.clean(title_str,'title'));
      else
        out.title_handle=[];
      end
      
      % enforce fontsize and paper size
      set(    out.axis_handle,          'FontSize',v.plot_fontsize_axis)
      set(get(out.axis_handle,'Title' ),'FontSize',v.plot_fontsize_title);
      set(get(out.axis_handle,'XLabel'),'FontSize',v.plot_fontsize_label);
      set(get(out.axis_handle,'YLabel'),'FontSize',v.plot_fontsize_label);

      %enforce grid
      if v.plot_grid
        grid(out.axis_handle,'on')
      end
      
      %enforce colormap
      if ~isempty(v.plot_colormap)
        out.colormap_handle=colormap(out.axis_handle,v.plot_colormap);
      else
        out.colormap_handle=[];
      end
      
      %enforce legend
      out.legend_handle=plotting.legend(varargin{:});
    end
    function legend_handle=legend(varargin)
      scale_legend_str=' x ';
      scale_legend_len=length(scale_legend_str);
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle',  gca,  @(i) ~isempty(i) && ishandle(i);...
      }},varargin{:});
      % check if any legend is given
      if ~isempty(v.plot_legend)
        %clean legend of _ and remove whatever is given in plot_legend_suppress
        v.plot_legend=cellfun(@(i) str.clean(i,[v.plot_legend_suppress(:);{'title'}]),v.plot_legend,'UniformOutput',false);
        %align scale and mean, if there
        legend_out=cell(numel(v.plot_legend),3);
        %determine locations of the scales
        scale_idx=strfind(v.plot_legend,scale_legend_str);
        %patch empty entries
        for i=1:numel(scale_idx)
          if isempty(scale_idx{i})
            scale_idx{i}=length(v.plot_legend{i});
          else
            scale_idx{i}=scale_idx{i}+scale_legend_len;
          end
        end
        %determine the locatios of biases
        bias_idx=cell(size(scale_idx));
        %patch empty entries
        for i=1:numel(bias_idx)
          %determine locations of the biases
          bias_idx{i}=max(strfind(strtrim(v.plot_legend{i}(1:(scale_idx{i}-scale_legend_len))),' '));
          if isempty(bias_idx{i})
            bias_idx{i}=scale_idx{i};
          else
            bias_idx{i}=bias_idx{i}+1;
          end
        end
        for i=1:numel(v.plot_legend)
          legend_out{i,3}=strtrim(v.plot_legend{i}(scale_idx{i}:end));
          legend_out{i,2}=strtrim(v.plot_legend{i}( bias_idx{i}:(scale_idx{i}-scale_legend_len)));
          legend_out{i,1}=strtrim(v.plot_legend{i}(        1:(bias_idx{i}-1)));
        end
        %get maximum length of each legend section and pad with blanks (right-align)
        for i=1:size(legend_out,2)
          section_max_len=max(cellfun(@(x) length(x),legend_out(:,i)));
          for j=1:size(legend_out,1)
            legend_out{j,i}=str.tabbed(legend_out{j,i},section_max_len,true);
          end
        end
        %assign column-aligned legend
        for i=1:size(legend_out,1)
          v.plot_legend{i}=strjoin(legend_out(i,:),' ');
        end
        %set the legend text
        legend_handle=legend(v.axis_handle,v.plot_legend);
        %adjust the location of the legend (even if the default 'best', it may happen that it is no longer in a good position)
        set(legend_handle,'location',v.plot_legend_location,'FontSize',v.plot_fontsize_legend,'FontName',v.plot_legend_fontname)
      else
        legend(v.axis_handle,'off')
        legend_handle=[];
      end
    end
    %% nice plotting stuff
    function out=dhist(x,y,z,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'x', @(i) isnumeric(i));
      p.addRequired( 'y', @(i) isnumeric(i) && numel(i)==numel(x));
      p.addRequired( 'z', @(i) isnumeric(i) && numel(i)==numel(x));
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap(...
        'parser',p,...
        'mandatory',{x,y,z},...
        'sources',{plotting.default,{...
          'xlabel','x',@(i) ischar(i);...
          'ylabel','y',@(i) ischar(i);...
          'zlabel','z',@(i) ischar(i);...
          'normalization','probability',@(i) ischar(i);...
          'xhistogramargs',{},@(i) iscell(i);...
          'yhistogramargs',{},@(i) iscell(i);...
          'zhistogramargs',{},@(i) iscell(i);...
          'markersize',    plotting.markersize(numel(x)),@(i) isnumeric(i) && isscalar(i);...
        }},varargin{:}...
      );
      %create figure
      out.fig_handle=plotting.figure(varargin{:});
      %plot histogram of x
      axis_handle=subplot(2,2,1);
      plot_handle=histogram(x(:),v.xhistogramargs{:},'Normalization',v.normalization);
      plot_handle.FaceColor=plot_handle.EdgeColor;
      out.xhist=plotting.enforce(v.varargin{:},...
        'axis_handle',axis_handle,...
        'plot_xdate',false,...
        'plot_xlabel',v.xlabel,...
        'plot_ylabel',get(plot_handle,'Normalization')...
      );
      out.xhist.plot_handle=plot_handle;
      %plot histogram of y
      axis_handle=subplot(2,2,2);
      plot_handle=histogram(y(:),v.yhistogramargs{:},'Normalization',v.normalization);
      plot_handle.FaceColor=plot_handle.EdgeColor;
      out.yhist=plotting.enforce(v.varargin{:},...
        'axis_handle',axis_handle,...
        'plot_xdate',false,...
        'plot_xlabel',v.ylabel,...
        'plot_ylabel',get(plot_handle,'Normalization')...
      );
      out.yhist.plot_handle=plot_handle;
      %plot histogram of z
      axis_handle=subplot(2,2,3);
      plot_handle=histogram(z(:),v.zhistogramargs{:},'Normalization',v.normalization);
      plot_handle.FaceColor=plot_handle.EdgeColor;
      out.zhist=plotting.enforce(v.varargin{:},...
        'axis_handle',axis_handle,...
        'plot_xdate',false,...
        'plot_xlabel',v.zlabel,...
        'plot_ylabel',get(plot_handle,'Normalization')...
      );
      out.zhist.plot_handle=plot_handle;
      %plot scatter
      axis_handle=subplot(2,2,4);
      plot_handle=scatter(x(:),y(:),v.markersize,z(:),'filled');
      out.scatter=plotting.enforce(v.varargin{:},...
        'axis_handle',axis_handle,...
        'plot_xdate',false,...
        'plot_xlabel',v.xlabel,...
        'plot_ylabel',v.ylabel...
      );
      out.scatter.plot_handle=plot_handle;
      out.scatter.cb_handle=colorbar;
      ylabel(out.scatter.cb_handle,[v.zlabel,' ',get(out.zhist.plot_handle,'Normalization')])
    end
  end
end
