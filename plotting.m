classdef plotting
  properties(Constant)
    %NOTE: access particular fields of the default parameters with:
    %out=varargs(plotting.default).pluck({'parameter1','parameter2'}).cell;
    default={...
'plot_legend',                 {},        @iscellstr;...
'plot_legend_suppress',        {},        @iscellstr;...
'plot_legend_replace',         {},        @iscellstr;...
'plot_legend_location',    'best',        @ischar;...    %'none' turns off the legend (can't be practically done with plot_legend)
'plot_legend_fontname','FixedWidth',      @ischar;...
'plot_legend_align',           '',        @ischar;...    %empty means no alignment
'plot_legend_align_scale_bias',false,     @str.islogical;... %this is done after plot_legend_align, destroying it
'plot_scale_legend_str',    ' x ',        @ischar;...
'plot_legend_align_str',       {},        @iscellstr;...
'plot_legend_align_right_just',true,      @str.islogical;...
'plot_legend_box',           true,        @str.islogical;...
        'plot_ylabel',         '',        @ischar;...
        'plot_xlabel',         '',        @ischar;...
        'plot_xdate',          false,     @str.islogical;...
        'plot_xdateformat',    '',        @ischar;...
        'plot_xlimits',        [-inf,inf],@(i) ( isnumeric(i) || isdatetime(i) || iscell(i) ) && numel(i)==2;...
        'plot_ylimits',        [-inf,inf],@(i) ( isnumeric(i) || iscell(i) ) && numel(i)==2;...
        'plot_set_axis_limits'  true,     @str.islogical;...
        'plot_size',    200+[0,0,21,9]*50,@(i) isnumeric(i) && numel(i)==4;...
        'plot_units',          'points',  @ischar;...
        'plot_visible',         true,     @str.islogical;...
        'plot_fontsize_axis',   18,       @num.isscalar;...
        'plot_fontsize_title',  24,       @num.isscalar;...
        'plot_fontsize_label',  20,       @num.isscalar;...
        'plot_fontsize_legend', 20,       @num.isscalar;...
        'plot_title',           '',       @ischar;...
        'plot_title_default',   '',       @ischar;...
        'plot_title_suppress',  {},       @iscellstr;...
        'plot_title_suffix',    '',       @ischar;...
        'plot_title_prefix',    '',       @ischar;...
        'plot_title_replace',   {},       @iscellstr;...
        'plot_grid',          true,       @str.islogical;...
        'plot_line_width',       2,       @isnumeric;...
        'plot_line_style',  'none',       @(i) ischar(i) || iscellstr(i);...
        'plot_line_color',  'none',       @ischar;...
        'plot_line_color_order', 0,       @isnumeric;...
        'plot_colormap',        '',       @(i) (isnumeric(i) && size(i,2)==3) || ischar(i) || ishandle(i);...
        'plot_psd',          false,       @str.islogical;...
        'plot_autoscale',    false,       @str.islogical;... %y-scale is derived from the data (in plotting.enforce)
        'plot_autoscale_factor', 4,       @num.isscalar;...
        'plot_automean',     false,       @str.islogical;... %middle-point of y axis is derived from the data (in plotting.enforce)
        'plot_zeromean',     false,       @str.islogical;... %mean of data is removed before plotting (in simpledata.plot)
        'plot_outlier_iter',     0,       @num.isscalar;...
        'plot_caxis',   [-inf,inf],       @(i) isnumeric(i) && numel(i)==2;...
        'plot_logy',         false,       @str.islogical;... 
        'plot_logx',         false,       @str.islogical;... 
        'plot_pause_on_save',false,       @str.islogical;... 
        'plot_save_fig',     false,       @str.islogical;... 
        'plot_save_mkdir',    true,       @str.islogical;... 
        'plot_save_ext',     'png',       @ischar;... 
        'plot_dir',file.orbdir('plot'),   @ischar;...
        'plot_save_data',    'yes',       @str.islogical;... %NOTICE: can be 'force' to recompute plot data (discards existing plot data)
        'plot_force',        false,       @str.islogical;... %NOTICE: this turns plot_save_data to 'force' if true (and force re-plotting, but needs to be implemented where needed)
    };
  end
  methods(Static)
    %% primitives
    function out=line_handles(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      tmp=get(axis_handle,'Children');
      c=0;out=[];
      for i=1:numel(tmp)
        if isprop(tmp(i),'LineWidth')
          c=c+1;
          out(c)=tmp(i); %#ok<AGROW>
        end
      end
    end
    function v=common_lim(axis_handle,mode)
      %convert to easier names
      switch lower(mode)
        case 'x'; mode='XData';
        case 'y'; mode='YData';
        case {'xdata','ydata'} %do nothing
        otherwise; error(['Cannot handle mode ''',mode,'''.'])
      end
      %set output
      v=[inf -inf];
      %get line handles
      lh=plotting.line_handles(axis_handle);
      %assume no date x data
      dates=false;
      %loop over all lines
      for i=1:numel(lh)
        %get x-data
        d=get(lh(i),mode);
        %clean up infs
        d(:)=d(isfinite(d(:)));
        if isdatetime(d)
          d=datenum(d);
          dates=true;
        end
        %clean up NaNs (datenum above converts NaTs to NaNs)
        d(:)=d(~isnan(d(:)));
        %get x-data extremeties
        mm=minmax(d);
        %accumulate inclusive x-axis limits
        v(1)=min([v(1),mm(1)]);
        v(2)=max([v(2),mm(2)]);
      end
      %return datetime type
      if dates
        v=datetime(v,'ConvertFrom','datenum');
      end
    end
    function out=line_color_map(mode,n)
      % handle inputs
      if ~exist('mode','var') || isempty(mode) || str.none(mode)
          mode='spiral';
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
        colormap_name=strrep(lower(mode),'-reversed','');
        branch=3;
      otherwise
        error(['Unknown mode ''',mode,'''.'])
      end
      %get colororder according to branch
      switch (branch)
      case 0
        colororder=[
          0     0     1 %blue
          1     0     0 %red
          0     1     0 %green
          0     1     1 %cyan
          1     0     1 %magenta
          1     1     0 %yellow
          0     0     0 %black
        ];
      case 1
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
      case 2
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
      case 3
        % colormaps name
        switch strrep(mode,'-reversed','')
        %some colormaps include the white color, need to remove it
        case{'gray','bone','pink','hot'}
          colororder=eval([colormap_name,'(',num2str(n+1),')']);
          colororder=colororder(1:end-1,:);
        otherwise
          colororder=eval([colormap_name,'(',num2str(n),')']);
        end
      end
      %make sure we got the right number of colors
      while size(colororder,1) < n
        colororder=[colororder;colororder]; %#ok<AGROW>
      end
      if size(colororder,1) > n
        colororder=colororder(1:n,:);
      end
      %transform colororder to cell array
      out=cell(size(colororder,1),1);
      for i=1:size(colororder,1)
          out{i}=colororder(i,:);
      end
      %enforce reverse order of colours
      if contains(mode,'-reversed')
        out=flipud(out);
      end
      %matlab orders lines in the reverse way
      out=flipud(out);
    end
    function line_color(mode,order,axis_handle)
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
      if ~exist('order','var') || isempty(order)
          order = -1;
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      %getting lines in this axis
      lines=findobj(axis_handle,'Type','line');
      %enforcing
      if isscalar(order)
        order=1:numel(lines);
      end
      %flip order upside down because matlab orders lines in the reverse way
      order=flipud(order(:));
      %get requested color scheme
      cm=plotting.line_color_map(mode,numel(lines));
      %loop over all lines or clr entries
      for i=1:min([numel(lines),numel(cm)])
        set(lines(i),'Color',cm{order(i)})
      end
    end
    function out=line_color_set_last(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
        axis_handle = gca;
      end
      lines=findobj(axis_handle,'Type','line');
      set(lines(2),'Color',get(lines(1),'Color'))
      if nargout>0; out=lines(1); end
    end
%     function line_color_reorder(order,axis_handle)
%       if ~exist('axis_handle','var') || isempty(axis_handle)
%         axis_handle = gca;
%       end
%       lines=findobj(axis_handle,'Type','line');
%       n=min([numel(order),numel(lines)]);
%       colors=flipud(arrayfun(@(i) get(i,'Color'),lines,'UniformOutput',false));
%       for i=1:n
%         set(lines(i),'Color',colors{order(i)})
%       end
%     end
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
    function line_markersize(w,axis_handle)
      % handle inputs
      if ~exist('w','var') || isempty(w)
          w=12;
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      %loop over all axis
      for ah = axis_handle(:)'
        lh=plotting.line_handles(ah);
        for i=1:numel(lh)		
          set(lh(i),'MarkerSize',w)
        end
      end
    end
    function line_style(line_str,axis_handle)
      % PLOT_LINE_LINESTYLE(LINE_STR,AXIS_HANDLE) sets the linestyle of the lines in
      % the specified plot.

      % Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

      %NOTICE: marker sizes must come first (since they are always only one char).
      %If there's no marker, use blank space.

      if ~exist('line_str','var') || isempty(line_str)
          line_styles={'-','-.','--'};
          line_ticks={'o','+','d','p','s','x'};
          line_str=cell(length(line_styles)*(length(line_ticks)+1),1);
          line_str(1:3)={' -',' -.',' --'};
          c=3;
          for i=1:length(line_styles)
              for j=1:length(line_ticks)
                  c=c+1;
                  line_str{c}=[line_ticks{j},line_styles{i}];
              end
          end
      end
      if ~iscellstr(line_str)
          error([mfilename,': expecting <line_str> to be a cell string, not a ',class(line_str),'.'])
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end

      lines=get(axis_handle,'Children');

      %matlab orders lines in the reverse way
      fix_idx=numel(lines):-1:1;

      for i=1:min([numel(lines),numel(line_str)])
          if ~isempty(line_str{i}) && ~isempty(line_str{i}(2:end))
              set(lines(fix_idx(i)),'LineStyle',line_str{i}(2:end))
              if line_str{i}(1) ~= ' '
                  try
                      set(lines(fix_idx(i)),'marker',line_str{i}(1))
                  catch
                      error([mfilename,': setting the marker failed. Be sure to input blank first character to lines without markers.'])
                  end
              end
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
    function font_size_xticks(fs,axis_handle)
      % handle inputs
      if ~exist('fs','var') || isempty(fs)
          fs=18;
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      %loop over all axis
      for h = axis_handle(:)'
        % https://www.mathworks.com/matlabcentral/answers/306225-how-can-i-change-the-font-size-of-the-tick-labels-without-changing-the-font-size-of-the-axis-labels
        hl = get(h,'XLabel');
        hlfs = get(hl,'FontSize');
        set(get(h,'XAxis'),'FontSize', fs)
        set(hl, 'FontSize', hlfs);
      end
    end
    function font_size_yticks(fs,axis_handle)
      % handle inputs
      if ~exist('fs','var') || isempty(fs)
          fs=18;
      end
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      %loop over all axis
      for h = axis_handle(:)'
        % https://www.mathworks.com/matlabcentral/answers/306225-how-can-i-change-the-font-size-of-the-tick-labels-without-changing-the-font-size-of-the-axis-labels
        hl = get(h,'YLabel');
        hlfs = get(hl,'FontSize');
        set(get(h,'YAxis'),'FontSize', fs)
        set(hl, 'FontSize', hlfs);
      end
    end
    function no_ticks(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
        axis_handle = gca;
      end
      for i={'x','y'}
        for j={'tick','ticklabel'}
          set(axis_handle,[i{1},j{1}],[])
        end
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
    function out=xaxis_tight(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      out=plotting.common_lim(axis_handle,'x');
      xlim(axis_handle,out);
    end
    function out=yaxis_tight(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      out=plotting.common_lim(axis_handle,'y');
      ylim(axis_handle,out);
    end
    function out=xlim(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      out=xlim(axis_handle);
      if any(~isfinite(out))
        out=plotting.common_lim(axis_handle,'x');
        delta=(out(2)-out(1))*0.05;
        out(1)=out(1)-delta;
        out(2)=out(2)+delta;
      end
    end
    function out=ylim(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      out=ylim(axis_handle);
      if any(~isfinite(out))
        out=plotting.common_lim(axis_handle,'y');
        delta=(out(2)-out(1))*0.05;
        out(1)=out(1)-delta;
        out(2)=out(2)+delta;
      end
    end
    function out=text(x,y,msg,varargin)
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle', gca, @(i) ~isempty(i) && ishandle(i);...
        'lines_y',      15, @(i) @num.isscalar;...
        'lines_x',      15, @(i) @num.isscalar;...
        'fontsize',     12, @(i) @num.isscalar;...
      }},varargin{:});
      %get axis
      a=axis;
      %need this to show text
      line_x   = (a(2)-a(1))/v.lines_x;
      corner_x = a(1)+line_x;
      line_y   = (a(4)-a(3))/v.lines_y;
      corner_y = a(4)-line_y;
      out=text(corner_x+(x-1)*line_x,corner_y-(y-1)*line_y,msg,'FontSize',v.fontsize,'fontname','FixedWidth');
    end
    function out=stats(dat,varargin)
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'format',  '%.3g', @(i) ~isempty(i) && ischar(i);...
        'x_pos',        1, @(i) @num.isscalar;...
        'y_pos',        0, @(i) @num.isscalar;...
        'stat_modes',{'std','mean','rms','min','max','numel'}, @iscellstr;...
      }},varargin{:});
      out=cell(size(v.stat_modes));
      good_idx=~isnan(dat);
      for i=1:numel(v.stat_modes)
        fh=str2func(v.stat_modes{i});
        out{i}=plotting.text(v.x_pos,v.y_pos+i,...
          [str.tabbed(v.stat_modes{i},5,true),' : ',num2str(fh(dat(good_idx)),v.format)]...
        );
      end
    end
    %% utility
    function fig_handle=figure(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default},varargin{:});
      %outputs
      fig_handle=figure('visible',str.logical(v.plot_visible,'onoff'));
      %plot size
%       set(fig_handle, 'Position',          v.plot_size,...
%                       'PaperUnits',        v.plot_units,...
%                       'PaperPosition',     v.plot_size);
      set(fig_handle, 'Position',          v.plot_size,...
                      'PaperUnits',        v.plot_units,...
                      'PaperSize',         [v.plot_size(3)-v.plot_size(1),v.plot_size(4)-v.plot_size(2)],...
                      'PaperPositionMode','auto',...
                      'color','w');
    end
    function out=enforce(varargin)
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap('sources',{
        plotting.default,...
        {...
          'axis_handle',  gca,  @(i) ~isempty(i) && ishandle(i);...
        }...
      },varargin{:});
      %resolve logicals
      v.plot_autoscale=str.logical(v.plot_autoscale);
      v.plot_automean =str.logical(v.plot_automean);
      % sanity
      if any(isfinite(v.plot_ylimits)) && ( v.plot_autoscale ||  v.plot_automean )
        error([mfilename,': option ''ylimits'' and ''autoscale'' or ''automean'' do not work concurrently.'])
      end
    
      %outputs
      out.axis_handle=v.axis_handle;

      % enforce fontsize and paper size
      set(    out.axis_handle,          'FontSize',v.plot_fontsize_axis)
      set(get(out.axis_handle,'Title' ),'FontSize',v.plot_fontsize_title);
      set(get(out.axis_handle,'XLabel'),'FontSize',v.plot_fontsize_label);
      set(get(out.axis_handle,'YLabel'),'FontSize',v.plot_fontsize_label);
      
      %enforce legend
      out.legend_handle=plotting.legend(varargin{:});

      
      % enforce line properties
      line_handles=plotting.line_handles(out.axis_handle);
      for i=1:numel(line_handles)
        set(line_handles(i),'LineWidth',v.plot_line_width)
      end
      if ~str.none(v.plot_line_style)
        plotting.line_style(v.plot_line_style,out.axis_handle)
      end
      if ~str.none(v.plot_line_color)
        plotting.line_color(v.plot_line_color,v.plot_line_color_order,out.axis_handle)
      end
      
      % start with x axis
      ax=plotting.xlim(out.axis_handle);
      % ensure numeric values
      v.plot_xlimits=cells.c2m(v.plot_xlimits);
      %check if dates are requested
      if (str.logical(v.plot_xdate) && ~str.logical(v.plot_psd)) || ...
          strcmp(v.plot_xlabel,'time') ||  ...
          strcmp(get(get(gca,'Xlabel'),'String'),'time')
        % enforce (possible) requested x-limits
        for i=1:2
          if isfinite(v.plot_xlimits(i))
            %NOTICE: this is needed because earlier versions of matlab
            %could not handle datetime in the x labels
            switch class(ax)
            case 'datetime'
              ax(i)=v.plot_xlimits(i);
            case 'double'
              ax(i)=datenum(v.plot_xlimits(i));
            otherwise
              error(['Cannot handle variable "ax" with class ',class(ax),...
                '; implementation needed'])
            end
          end
        end
        % set auto x-label, unless one is explicity given
        if strcmp(v.plot_xlabel,'time')
          if ~strcmp(datestr(ax(1),'yyyymmdd'),datestr(ax(2),'yyyymmdd')) && ...
            (~strcmp(datestr(ax(2),'HHMMSS'),'000000') || ax(2)-ax(1)>1)
            out.xlabel_handle=xlabel(out.axis_handle,...
              [datestr(ax(1),'yyyy-mm-dd'),' to ',datestr(ax(2),'yyyy-mm-dd')]...
            );
          else
            out.xlabel_handle=xlabel(out.axis_handle,datestr(ax(1)));
          end
        end
        %enforce requested x date tick format
        if ~str.none(v.plot_xdateformat)
          if ~isempty(v.plot_xdateformat)
            datetick(out.axis_handle,'x',v.plot_xdateformat)
          else
            datetick(out.axis_handle,'x')
          end
        end
      else
        % enforce (possible) requested x-limits
        for i=1:2
          if isfinite(v.plot_xlimits(i)) && isnumeric(v.plot_xlimits(i))
            ax(i)=v.plot_xlimits(i);
          end
        end
      end
      if ~str.none(v.plot_set_axis_limits)     
        % set axis limits (can be the ones matlab so wisely set)
        xlim(out.axis_handle,ax);
      end
      
      % enforce requested y-limits
      ay=plotting.ylim(out.axis_handle);
      % ensure numeric values
      v.plot_xlimits=cells.c2m(v.plot_xlimits);
      for i=1:2
        if isfinite(v.plot_ylimits(i))
          ay(i)=   v.plot_ylimits(i);
        end
      end
      % enforce auto-scale and/or auto-mean
      % TODO: add detrending here (to go along plot_outlier_iter, as is done in simpledata.outlier)
      if v.plot_automean || v.plot_autoscale || v.plot_outlier_iter>0
        %gather plotted data
        dat=cell(size(line_handles));
        for i=1:numel(line_handles)
          tmp=get(line_handles(i),'ydata');
          dat{i}=transpose(tmp(:));
        end
        dat=[dat{:}];
        dat=dat(~isnan(dat(:)));
        %remove outliers before computing axis limits (if requested)
        for c=1:v.plot_outlier_iter
          dat=simpledata.rm_outliers(dat);
        end
        dat=dat(~isnan(dat(:)));
        if ~isempty(dat) && diff(minmax(dat))~=0
          %enforce data-driven mean and/or scale
          if v.plot_automean && v.plot_autoscale
            ay=mean(dat)+v.plot_autoscale_factor*std(dat)*[-1,1];
          elseif v.plot_automean
            ay=minmax(dat);
          elseif v.plot_autoscale
            ay=mean(ay)+v.plot_autoscale_factor*std(dat)*[-1,1];
          end
          %fix non-negative data sets
          if all(dat>=0) && ay(2)<0
            ay(2)=0;
          end
        end
      end
      % set axis limits (can be the ones matlab so wisely set)
      ylim(out.axis_handle,ay);
      
      %enforce labels
      if isempty(v.plot_xlabel)
        %do nothing
      elseif ~str.none(v.plot_xlabel)
        out.xlabel_handle=xlabel(out.axis_handle,v.plot_xlabel);
      else
        % if ~isfield(out,'xlabel_handle')
          out.xlabel_handle=[];
        % end
        xlabel('')
      end
      if isempty(v.plot_ylabel)
        %do nothing
      elseif ~str.none(v.plot_ylabel)
        out.ylabel_handle=ylabel(out.axis_handle,v.plot_ylabel);
      else
        % if ~isfield(out,'ylabel_handle')
          out.ylabel_handle=[];
        % end
        ylabel('')
      end
      
      %enforce title
      if str.none(v.plot_title)
        out.title_handle=[];
        title('')
      else
        if str.default(v.plot_title)
          title_str=v.plot_title_default;
        else
          title_str=v.plot_title;
        end
        %operate and put it there
        out.plot_title=plotting.title_replace_clean(varargin{:},'plot_title',title_str);
        out.title_handle=title(out.plot_title);
      end

      %enforce grid
      if str.logical(v.plot_grid)
        grid(out.axis_handle,'on')
      end
      
      %enforce colormap
      if str.none(v.plot_colormap)
        out.colormap=[];
      else
        %enforce caxis limits
        if any(isfinite(v.plot_caxis))
          ca=caxis;
          for i=1:numel(ca)
            if isfinite(v.plot_caxis(i)); ca(i)=v.plot_caxis(i);end
          end
          caxis(ca);
        end
        if ischar(v.plot_colormap) && ~isempty(v.plot_colormap)
          %get list of options for the colormap: these can be added to the names of generic colormaps, e.g.:
          % 'jetzero', 'optbone'
          clist={'opt','zero'};
          %determine if any option is being used
          cuse=cellfun(@(i) str.contains(v.plot_colormap,i),clist);
          %clean the colormap of options
          if any(cuse); v.plot_colormap=str.clean(v.plot_colormap,clist(cuse)); end
          %get the requested colormap
          out.colormap=eval([v.plot_colormap,'(1000)']);
        else
          %any fine-tuned colormap will be distorted by the caxis command above
          if any(isfinite(v.plot_caxis))
            disp('WARNING: numeric plot_colormap doesn''t play well with plot_caxis for non-generic colormaps; for generic ones, you might as well pass the colormap names');
          end
          clist={'unknown'};
          cuse=false;
          %fallback to a defaule
          if isempty(v.plot_colormap)
            out.colormap='jet';
          else
            out.colormap=v.plot_colormap;
          end
        end
        %apply options
        for c=1:numel(clist)
          if ~cuse(c); continue; end
          switch clist{c}
          case 'opt'
            out.colormap=cb.opt(out.colormap,out.axis_handle);
          case 'zero' 
            clim=caxis;
            if clim(1)<0 && clim(2)>0; czero=0;
            elseif all(clim<=0);        czero=max(clim);
            elseif all(clim>=0);        czero=min(clim);
            else
              error('BUG TRAP: debug needed, this should not happen')
            end
            out.colormap=cb.zero_white(czero,out.colormap);
          end
        end
        %apply the modified colormap
        out.colormap_handle=colormap(out.axis_handle,out.colormap);
      end

      %enforce axis type
      if str.logical(v.plot_logy)
        set(gca, 'YScale', 'log')
      end
      if str.logical(v.plot_logx)
        set(gca, 'XScale', 'log')
      end
      
    end
    function out=title_replace_clean(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
      }},varargin{:});
      %suppress some parts, if requested
      out=strsplit(v.plot_title,{' '});
      %add prefix and suffix
      out=strjoin([{v.plot_title_prefix};out(:);{v.plot_title_suffix}],' ');
      %replace explicit strings
      if ~isempty(v.plot_title_replace)
        out=str.rep(out,v.plot_title_replace{:});
      end
      %clean up the tile 
      out=str.clean(out,[v.plot_title_suppress(:);{'title'}]);
    end
    function out=legend_replace_clean(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
      }},varargin{:});
      %replace explicit strings (keep this before the legend-cleaning bit so the strings make sense outside)
      if ~isempty(v.plot_legend_replace)
        out=cellfun(@(i) str.rep(i,v.plot_legend_replace{:}),v.plot_legend,'UniformOutput',false);
      else
        out=v.plot_legend;
      end
      %clean legend of _ and remove whatever is given in plot_legend_suppress
      out=cellfun(@(i) str.clean(i,[v.plot_legend_suppress(:);{'title'}]),out,'UniformOutput',false);
    end
    function legend_handle=legend(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle',  gca,  @(i) ~isempty(i) && ishandle(i);...
      }},varargin{:});
      % check if any legend is given
      if str.none(v.plot_legend_location)
        legend(v.axis_handle,'off')
        legend_handle=[];
        return
      else
        v.plot_legend=plotting.legend_replace_clean(varargin{:});
        %maybe there's the need to align words
        if ~isempty(v.plot_legend_align)
          %keywords
          switch lower(v.plot_legend_align)
          case {'space','spaces','blank','blanks'}
            v.plot_legend_align=' ';
          end
          v.plot_legend=str.columns(v.plot_legend,'right',v.plot_legend_align);
        end
        %maybe there's the need to make the legend nicely aligned
        if numel(v.plot_legend_align_str)>0
          I_align=numel(v.plot_legend);
          J_align=numel(v.plot_legend_align_str);
          %handle default right-justify
          switch numel(v.plot_legend_align_right_just)
          case 1
            v.plot_legend_align_right_just=true(1,J_align) && str.logical(v.plot_legend_align_right_just);
          case J_align+1
            v.plot_legend_align_right_just=transpose(str.logical(v.plot_legend_align_right_just(:)));
          otherwise
            error(['plot_legend_align_right_just must have either 1 or ',num2str(numel(J_align)+1),...
              ' entries, not ',num2str(numel(v.plot_legend_align_right_just)),'.'])
          end
          %align the legend entries
          legend_out=cell(I_align,J_align+1);
           align_idx=cell(I_align,J_align);
          %get the align indexes
          for j=1:J_align
            %determine locations of this align string
            align_idx(:,j)=strfind(v.plot_legend,v.plot_legend_align_str{j});
            %patch empty entries and correct non-empty
            for i=1:I_align
              if isempty(align_idx{i,j})
                align_idx{i,j}=length(v.plot_legend{i});
              else
                align_idx{i,j}=align_idx{i,j}-1;
              end
            end
          end
          %slice the legend
          for i=1:I_align
            for j=1:J_align+1
              switch j
              case 1;         i_start=1;                 i_stop=align_idx{i,j};
              case J_align+1; i_start=align_idx{i,j-1}+1;i_stop=length(v.plot_legend{i});
              otherwise;      i_start=align_idx{i,j-1}+1;i_stop=align_idx{i,j};
              end
              legend_out{i,j}=strtrim(v.plot_legend{i}(i_start:i_stop));
            end
          end
          %get maximum length of each legend section and pad with blanks (right-align)
          for j=1:J_align+1
            section_max_len=max(cellfun(@(x) length(x),legend_out(:,j)));
            for i=1:I_align
              legend_out{i,j}=str.tabbed(legend_out{i,j},section_max_len,v.plot_legend_align_right_just(j));
            end
          end
          %assign column-aligned legend
          for i=1:I_align
            v.plot_legend{i}=strjoin(legend_out(i,:),' ');
          end
        end
        %TODO: see if the code above can be used for the purpose below
        %maybe there's the need to make the scale and bias nicely aligned
        if str.logical(v.plot_legend_align_scale_bias)
          %align scale and mean, if there
          legend_out=cell(numel(v.plot_legend),3);
          %determine locations of the scales
          scale_idx=strfind(v.plot_legend,v.plot_scale_legend_str);
          %patch empty entries
          for i=1:numel(scale_idx)
            if isempty(scale_idx{i})
              scale_idx{i}=length(v.plot_legend{i});
            else
              scale_idx{i}=scale_idx{i}+length(v.plot_scale_legend_str);
            end
          end
          %determine the locatios of biases
          bias_idx=cell(size(scale_idx));
          %patch empty entries
          for i=1:numel(bias_idx)
            %determine locations of the biases
            bias_idx{i}=max(strfind(strtrim(v.plot_legend{i}(1:(scale_idx{i}-length(v.plot_scale_legend_str)))),' '));
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
        end
      end
      %set the legend text (some plots don't have lines, so legends don't make sense)
      warning('off','MATLAB:legend:PlotEmpty')
      legend_handle=legend(v.axis_handle,v.plot_legend);
      warning('on','MATLAB:legend:PlotEmpty')
      %adjust the location of the legend (even if the default 'best', it may happen that it is no longer in a good position)
      set(legend_handle,...
        'location',v.plot_legend_location,...
        'box',str.logical(v.plot_legend_box,'onoff'),...
        'FontSize',v.plot_fontsize_legend,...
        'FontName',v.plot_legend_fontname...
        )
    end
    function save(filename,varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'fig_handle',  gcf,  @(i) ~isempty(i) && ishandle(i);...
      }},varargin{:});
      if str.logical(v.plot_pause_on_save); keyboard; end
      [p,n,e]=fileparts(filename);
      if strcmp(e,'.') || isempty(e)
        e=['.',strrep(v.plot_save_ext,'.','')];
      end
      if str.logical(v.plot_save_mkdir)
        file.mkdir(p);
      end
      if str.logical(v.plot_save_fig)
        saveas(v.fig_handle,fullfile(p,n),'fig')
      end
      saveas(v.fig_handle,fullfile(p,[n,e]))
      str.say('Plotted',fullfile(p,[n,e]))
    end
    function [loaddata,savedata]=forcing(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
      }},varargin{:});
      %don't save data by default
      savedata=false;
      %override plot_save_data when force is true
      if v.plot_force; v.plot_save_data='force'; end
      %check if loading the data is possible
      try
        loaddata=str.logical(v.plot_save_data);
      catch
        switch lower(v.plot_save_data) 
        case 'force'
          loaddata=false;
          savedata=true;
        otherwise
          error(['Cannot handle parameter ''plot_save_data'' with value ''',v.plot_save_data,'''.'])
        end
      end
    end
    function [axis_handle,l,w]=subplot(i,n,portrait_flag)
      if ~exist('portrait_flag','var') || isempty(portrait_flag)
        portrait_flag=true;
      end
      l=floor(sqrt(n)); w=ceil(n/l);
      if portrait_flag
        if l<w
          tmp=w;w=l;l=tmp;
        end
      else %landscape
        if l>w
          tmp=w;w=l;l=tmp;
        end
      end
      if 1<=i && i<=n
        axis_handle=subplot(l,w,i);
      else
        axis_handle=[];
      end
    end
    %% nice plotting stuff
    function out=hist(x,varargin)
      p=machinery.inputParser;
      p.addRequired( 'x', @isnumeric);
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap(...
        'parser',p,...
        'mandatory',{x},...
        'sources',{plotting.default,{...
          'normalization','probability',@ischar;...
          'edgecolor','black',@ischar;...
          'facecolor','black',@ischar;...
          'histogramargs',{},@iscell;...
          'markersize',    plotting.markersize(numel(x)),@num.isscalar;...
        }},varargin{:}...
      );
      %create figure
      fig_handle=plotting.figure(varargin{:});
      %plot histogram of x
      plot_handle=histogram(x(:),...
        v.histogramargs{:},...
        'Normalization',v.normalization,...
        'edgecolor',v.edgecolor,...
        'facecolor',v.facecolor);
      %enforece plot and save returned struct to outputs
      out=plotting.enforce(v.varargin{:},...
        'axis_handle',gca,...
        'plot_xdate',false,...
        'plot_ylabel',get(plot_handle,'Normalization')...
      );
      %add stats (if requested, as handled in plotting.stats)
      out.stats=plotting.stats(x,varargin{:});
      %append handles
      out.plot_handle=plot_handle;
      out.axis_handle=gca;
      out.fig_handle=fig_handle;
    end
    function out=dhist(x,y,z,varargin)
      p=machinery.inputParser;
      p.addRequired( 'x', @isnumeric);
      p.addRequired( 'y', @(i) isnumeric(i) && numel(i)==numel(x));
      p.addRequired( 'z', @(i) isnumeric(i) && numel(i)==numel(x));
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap(...
        'parser',p,...
        'mandatory',{x,y,z},...
        'sources',{plotting.default,{...
          'xlabel','x',@ischar;...
          'ylabel','y',@ischar;...
          'zlabel','z',@ischar;...
          'normalization','probability',@ischar;...
          'xhistogramargs',{},@iscell;...
          'yhistogramargs',{},@iscell;...
          'zhistogramargs',{},@iscell;...
          'markersize',    plotting.markersize(numel(x)),@num.isscalar;...
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
