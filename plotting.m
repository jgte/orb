classdef plotting
  properties(Constant)
    %NOTE: access particular fields of the default parameters with:
    %out=varargs(plotting.default).pluck({'parameter1','parameter2'}).cell;
    default={...
        'plot_file_prefix',    '',        @(i) ischar(i);...
        'plot_file_suffix',    '',        @(i) ischar(i);...
        'plot_legend',         {},        @(i) iscellstr(i);...
        'plot_legend_suppress',{},        @(i) iscellstr(i);...
        'plot_legend_replace', {},        @(i) iscellstr(i);...
        'plot_legend_location','best',    @(i) ischar(i);...
     'plot_legend_fontname','FixedWidth', @(i) ischar(i);...
        'plot_legend_align',   '',        @(i) ischar(i);...    %empty means no alignment
     'plot_legend_align_scale_bias',false,@(i) islogical(i);... %this is done after plot_legend_align, destroying it
     'plot_scale_legend_str',' x ',       @(i) ischar(i);...
        'plot_ylabel',         '',        @(i) ischar(i);...
        'plot_xlabel',         '',        @(i) ischar(i);...
        'plot_xdate',          false,     @(i) islogical(i);...
        'plot_xdateformat',    '',        @(i) ischar(i);...
        'plot_xlimits',        [-inf,inf],@(i) isnumeric(i) && numel(i)==2;...
        'plot_ylimits',        [-inf,inf],@(i) isnumeric(i) && numel(i)==2;...
        'plot_size',    200+[0,0,21,9]*50,@(i) isnumeric(i) && numel(i)==4;...
        'plot_units',          'points',  @(i) ischar(i);...
        'plot_visible',         true,     @(i) islogical(i);...
        'plot_fontsize_axis',   18,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_fontsize_title',  24,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_fontsize_label',  20,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_fontsize_legend', 20,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_title',           '',       @(i) ischar(i);...
        'plot_title_suppress',  {},       @(i) iscellstr(i);...
        'plot_title_suffix',    '',       @(i) ischar(i);...
        'plot_title_prefix',    '',       @(i) ischar(i);...
        'plot_title_replace',   {},       @(i) iscellstr(i);...
        'plot_grid',            true,     @(i) islogical(i);...
        'plot_line_width',      2,        @(i) isnumeric(i);...
        'plot_line_style',  'none',       @(i) ischar(i) || iscellstr(i);...
        'plot_colormap',        '',       @(i) ischar(i) || ishandle(i);...
        'plot_psd',          false,       @(i) islogical(i);...
        'plot_autoscale',    false,       @(i) islogical(i);... %y-scale is derived from the data (in plotting.enforce)
        'plot_automean',     false,       @(i) islogical(i);... %middle-point of y axis is derived from the data (in plotting.enforce)
        'plot_zeromean',     false,       @(i) islogical(i);... %mean of data is removed before plotting (in simpledata.plot)
        'plot_outlier',          0,       @(i) isnumeric(i) && isscalar(i); ...
        'plot_caxis',   [-inf,inf],       @(i) isnumeric(i) && numel(i)==2;...
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
      v=[-inf inf -inf inf];
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
        colormap_name=lower(mode);
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
      v1=axis(axis_handle);
      axis(axis_handle,'tight');
      v2=axis;
      axis(axis_handle,[v2(1:2),v1(3:4)]);
      if nargout>0
        out=axis;
      end
    end
    function out=yaxis_tight(axis_handle)
      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gca;
      end
      v1=axis(axis_handle);
      axis(axis_handle,'tight');
      v2=axis;
      axis(axis_handle,[v1(1:2),v2(3:4)]);
      if nargout>0
        out=axis;
      end
    end
    function out=text(x,y,msg,varargin)
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle', gca, @(i) ~isempty(i) && ishandle(i);...
        'lines_y',      15, @(i) isnumeric() && isscalar(i);...
        'lines_x',      15, @(i) isnumeric() && isscalar(i);...
        'fontsize',     12, @(i) isnumeric() && isscalar(i);...
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
        'x_pos',        1, @(i) isnumeric() && isscalar(i);...
        'y_pos',        0, @(i) isnumeric() && isscalar(i);...
        'stat_modes',{'std','mean','rms','min','max','numel'}, @(i) iscellstr(i);...
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
      set(fig_handle, 'Position',          v.plot_size,...
                      'PaperUnits',        v.plot_units,...
                      'PaperPosition',     v.plot_size);
    end
    function out=enforce(varargin)
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle',  gca,  @(i) ~isempty(i) && ishandle(i);...
      }},varargin{:});
      % sanity
      if any(isfinite(v.plot_ylimits)) && ( v.plot_autoscale ||  v.plot_automean )
        error([mfilename,': option ''ylimits'' and ''autoscale'' or ''automean'' do not work concurrently.'])
      end
    
      %outputs
      out.axis_handle=v.axis_handle;
            
      % enforce line properties
      line_handles=plotting.line_handles(out.axis_handle);
      for i=1:numel(line_handles)
        set(line_handles(i),'LineWidth',v.plot_line_width)
      end
      if ~strcmpi(v.plot_line_style,'none')
        plotting.line_style(v.plot_line_style,out.axis_handle)
      end
      
      % start with current axis
      a=axis(out.axis_handle);
      %check if dates are requested
      if (v.plot_xdate && ~v.plot_psd) || ...
          strcmp(v.plot_xlabel,'time') ||  ...
          strcmp(get(get(gca,'Xlabel'),'String'),'time')
        % enforce (possible) requested x-limits
        for i=1:2
          if isfinite(v.plot_xlimits(i))
            a(i)=datenum(v.plot_xlimits(i));
          end
        end
        % set auto x-label, unless one is explicity given
        if isempty(v.plot_xlabel) || strcmp(v.plot_xlabel,'time')
          if ~strcmp(datestr(a(1),'yyyymmdd'),datestr(a(2),'yyyymmdd')) && ...
            (~strcmp(datestr(a(2),'HHMMSS'),'000000') || a(2)-a(1)>1)
            out.xlabel_handle=xlabel(out.axis_handle,...
              [datestr(a(1),'yyyy-mm-dd'),' to ',datestr(a(2),'yyyy-mm-dd')]...
            );
          else
            out.xlabel_handle=xlabel(out.axis_handle,datestr(a(1)));
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
      % enforce auto-scale and/or auto-mean
      if v.plot_automean || v.plot_autoscale || v.plot_outlier>0
        %gather plotted data
        dat=cell(size(line_handles));
        for i=1:numel(line_handles)
          tmp=get(line_handles(i),'ydata');
          dat{i}=transpose(tmp(:));
        end
        dat=[dat{:}];
        dat=dat(~isnan(dat(:)));
        %remove outliers before computing axis limits (if requested)
        for c=1:v.plot_outlier
          dat=simpledata.rm_outliers(dat);
        end
        dat=dat(~isnan(dat(:)));
        if ~isempty(dat) && diff(minmax(dat))~=0
          %enforce data-driven mean and/or scale
          if v.plot_automean && v.plot_autoscale
            a(3:4)=mean(dat)+4*std(dat)*[-1,1];
          elseif v.plot_automean
            a(3:4)=minmax(dat);
          elseif v.plot_autoscale
            a(3:4)=mean(a(3:4))+4*std(dat)*[-1,1];
          end
          %fix non-negative data sets
          if all(dat>=0) && a(3)<0
            a(3)=0;
          end
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
      switch lower(v.plot_title)
      case {'','clear'}
        out.title_handle=[];
      otherwise
        %suppress some parts, if requested
        title_str=setdiff(strsplit(v.plot_title,{' ','.'}),v.plot_title_suppress,'stable');
        %add prefix and suffix
        title_str=strjoin([{v.plot_title_prefix};title_str(:);{v.plot_title_suffix}],' ');
        %propagate
        v.plot_title=title_str;
        %replace explicit strings
        if ~isempty(v.plot_title_replace)
          v.plot_title=str.rep(v.plot_title,v.plot_title_replace{:});
        end
        %clean up the tile put it there
        out.title_handle=title(out.axis_handle,str.clean(v.plot_title,'title'));
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
      
      %enforce caxis limits
      if any(isfinite(v.plot_caxis))
        ca=caxis;
        for i=1:numel(ca)
          if isfinite(v.plot_caxis(i)); ca(i)=v.plot_caxis(i);end
        end
        caxis(ca);
      end
      
      %enforce legend
      out.legend_handle=plotting.legend(varargin{:});
    end
    function legend_handle=legend(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{plotting.default,{...
        'axis_handle',  gca,  @(i) ~isempty(i) && ishandle(i);...
      }},varargin{:});
      % check if any legend is given
      if isempty(v.plot_legend) || (...
          ~isempty(v.plot_legend{1}) && (...
            strcmpi(v.plot_legend{1},'clean') || ...
            strcmpi(v.plot_legend{1},'clear') || ...
            strcmpi(v.plot_legend{1},'none')...
        ))
        legend(v.axis_handle,'off')
        legend_handle=[];
        return
      else
        %clean legend of _ and remove whatever is given in plot_legend_suppress
        v.plot_legend=cellfun(@(i) str.clean(i,[v.plot_legend_suppress(:);{'title'}]),v.plot_legend,'UniformOutput',false);
        %replace explicit strings
        if ~isempty(v.plot_legend_replace)
          v.plot_legend=cellfun(@(i) str.rep(i,v.plot_legend_replace{:}),v.plot_legend,'UniformOutput',false);
        end
        %maybe there's the need to align words
        if ~isempty(v.plot_legend_align)
          %keywords
          switch lower(v.plot_legend_align)
          case {'space','spaces','blank','blanks'}
            v.plot_legend_align=' ';
          end
          v.plot_legend=str.columns(v.plot_legend,'right',v.plot_legend_align);
        end
        %maybe there's the need to make the scale and bias nicely aligned
        if v.plot_legend_align_scale_bias
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
      %set the legend text
      legend_handle=legend(v.axis_handle,v.plot_legend);
      %adjust the location of the legend (even if the default 'best', it may happen that it is no longer in a good position)
      set(legend_handle,'location',v.plot_legend_location,'FontSize',v.plot_fontsize_legend,'FontName',v.plot_legend_fontname)
    end
    %% nice plotting stuff
    function out=hist(x,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'x', @(i) isnumeric(i));
      % add input arguments to collection of parameters 'v'
      v=varargs.wrap(...
        'parser',p,...
        'mandatory',{x},...
        'sources',{plotting.default,{...
          'normalization','probability',@(i) ischar(i);...
          'edgecolor','black',@(i) ischar(i);...
          'facecolor','black',@(i) ischar(i);...
          'histogramargs',{},@(i) iscell(i);...
          'markersize',    plotting.markersize(numel(x)),@(i) isnumeric(i) && isscalar(i);...
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
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'x', @(i) isnumeric(i));
      p.addRequired( 'y', @(i) isnumeric(i) && numel(i)==numel(x));
      p.addRequired( 'z', @(i) isnumeric(i) && numel(i)==numel(x));
      % add input arguments to collection of parameters 'v'
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
