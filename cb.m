classdef cb
  methods(Static)
    function label(str,axis_handle)
      % PLOT_COLORBAR_LABEL(STR) sets the colorbar label to the specified string.
      %
      %   PLOT_COLORBAR_LABEL(STR,H) applies the label to axess H.

      % Created by P.Inacio <p.m.g.inacio@tudelft.nl>

      if ~exist('axis_handle','var') || isempty(axis_handle)
          axis_handle = gcf;
      elseif ~isfigure(axis_handle)
          error('plot_colorbar_label:inputs','''axis_handle'' must be a valid figure handle')
      end

      % Get colorbar handle
      cb_handle=findobj(axis_handle,'tag','Colorbar');
      if ~isempty(cb_handle)
          % find colorbar position
          if any(strcmpi(get(cb_handle,'Location'),{'north','south','northoutside','southoutside'}))
              xlabel(cb_handle,str)
          else
              ylabel(cb_handle,str)
          end
      end
    end
    function new=zero_white(WHITEVAL,old,gap)
      % PLOT_COLORBAR_ZERO_WHITE(WHITEVAL) sets the WHITEVAL in a colormar to
      % white.
      %
      %   PLOT_COLORBAR_ZERO_WHITE(WHITEVAL,OLD) allows you to specify which
      %   colormap to work on.
      %
      %   PLOT_COLORBAR_ZERO_WHITE(WHITEVAL,OLD,GAP) additionally allows the size
      %   of the gap to be is controled. GAP*100 is the number of color shades to
      %   be whitened, in each direction. Default value is 0.1, meaning that, if
      %   colorbar has 100 tiles, then there will be 20 white-ish color shades.
      %   Values > 0.5 do not make sense.

      % Created by J.Encarnacao, <J.G.deTeixeiradaEncarnacao@tudelft.nl>

      if ~exist('WHITEVAL','var') || isempty(WHITEVAL)
          WHITEVAL=0;
      end
      if ~exist('old','var') || isempty(old)
          old=colormap;
      end

      %getting min/max values of current colorscale
      tmp=caxis;
      MINVAL=tmp(1);
      MAXVAL=tmp(2);
      NELEM=size(old,1);

      if (WHITEVAL < MINVAL) || (WHITEVAL > MAXVAL)
          disp([mfilename,':WARNING: requested white zero is out range.'])
          if nargout == 0
              %setting colormap
              colormap(old)
          else
              new=old;
          end
          return
      end

      %% parameters

      %this is the length of the gap around white. This gap is produced by
      %whitening the colors around <WHITEVAL>. <gap*100> is the number of color
      %shades to be whitened, in each direction. If <NELEM> is 100, then there
      %will be 20 white-ish color shades. Values > 0.5 do not make sense.
      if ~exist('gap','var') || isempty(gap)
          gap=0.1;
      end

      %% calcs

      %propagating
      new=old;

      %old colorbar domain
      x_old=linspace(MINVAL,MAXVAL,NELEM);

      %searching for colormap index that is closest to zero
      zero_idx=find(abs(x_old-WHITEVAL)==min(abs(x_old-WHITEVAL)),1,'first'); %This puts white on zero
                                                           % without stretching
                                                           % to colormap
      %making it white
      new(zero_idx,:)=1;

      %% whitening the gap

      if gap > 0
          %number of gap indexes
          gap_idx_nr=ceil(gap*NELEM);

          %building index vectors
          gap_idx_lower=zero_idx-gap_idx_nr;
          gap_idx_upper=zero_idx+gap_idx_nr;

          for i=1:3
              if gap_idx_lower > 0
                 new(gap_idx_lower:zero_idx,i)=interp1([gap_idx_lower,zero_idx],[new(gap_idx_lower,i), 1],...
                     gap_idx_lower:zero_idx,'linear');
              end
              if gap_idx_upper > 0
                 new(zero_idx:gap_idx_upper,i)=interp1([zero_idx,gap_idx_upper],[1,new(gap_idx_upper,i)],...
                     zero_idx:gap_idx_upper,'linear');
              end
          end
      end

      %% outputs

      if nargout == 0
          %setting colormap
          colormap(new)
      end
    end
    function nan(color)
      % PLOT_COLORBAR_NAN(COLOR) sets the NaN's in the current figure to be of the
      % specified color. COLOR can be either a color character or an RGB triplet.
      %
      %   NaNs im plots take the color of the minimum value in the color scale.
      %   When NaNs are the same color as the minimum values of the plotted data
      %   visualization becomes confusing.
      %   To address this problem, another color a new color level is added to
      %   the colormap with a different color specifically for NaNs.

      % Created by P.Inacio <p.m.g.inacio@tudelft.nl>

      % Check inputs
      if ~exist('color','var') || isempty(color)
          % default color is white
          color = [ 1 1 1 ];
      end

      % Translate to color if char
      if ischar(color)
          color = plot_color2RGB(color);
      end

      % get the colormap and current axes
      cm=colormap;
      ca=caxis;

      % get the step size of the colorbar
      dc=diff(ca)/length(cm);

      % add one step to the minimim color level
      ca(1)=ca(1)-dc;

      % add color to colormap
      cm = [ color ; cm ];

      % apply new colormap and color axis
      colormap(cm);
      caxis(ca);
    end
    function new=manual(domain,old)
      %sets the <old> colormap (if ignored, the current one is considered) to be
      %rearranged to fit the requested <domain>. If <domain> is -10:1:10, then
      %the data is shown in [-10,10] in a linear fashion, i.e. there are as many
      %colors in [-10,-9] as in [3,4]. If <domain> is [-10:-1,-0.9:0.1:0.9,1:10],
      %then there is an (almost) equal number of colors in [-1,1] as in the rest
      %of the domain.
      %
      %This is useful to highlight a particular region of the data.
      %
      %J.G.deTeixeiradaEncarnacao@tudelft.nl

      if ~exist('old','var') || isempty(old)
          old=colormap;
      end

      if min(size(domain)) > 1
          error([mfilename,': input <domain> must be a vector.'])
      end
      if ~issorted(domain)
          error([mfilename,': input <domain> must be sorted.'])
      end
      if size(domain,1) > 1
          domain=domain';
      end

      %getting min/max values of current colorscale
      MINVAL=min(domain);
      MAXVAL=max(domain);

      %dealing with homogeneous domains
      if size(domain,2) == 2
          domain=[MINVAL,mean([MINVAL,MAXVAL]),MAXVAL];
      end

      %setting colormap scale
      caxis([MINVAL,MAXVAL])

      NELEM=size(old,1);

      %% building colormap to fit requested domain

      %determining old domain
      tmp=caxis;
      x_old=linspace(tmp(1),tmp(2),NELEM);

      x_h=mean([domain(1:end-1);domain(2:end)]);
      h=1./diff(domain);

      %resampling
      x_h_new=linspace(MINVAL,MAXVAL,NELEM);
      h=interp1(x_h,h,x_h_new,'linear','extrap');

      %weights are non-dimensional histogram
      w=h/sum(h);
      %removing zeros
      w(w==0)=[];
      %calculating interval size
      a=w*(MAXVAL-MINVAL);
      %calculating new domain
      x_new=[0,cumsum(a)]+MINVAL;

      %resampling
      new=zeros(length(x_new),3);
      for i=1:3
          new(:,i)=interp1(x_old,old(:,i),x_new,'linear','extrap');
      end

      new(new>1)=1;
      new(new<0)=0;

      if nargout == 0
          %setting colormap
          colormap(new)
      end
    end
    function new=center(WHITEVAL,old,gap)
      % NEW = PLOT_COLORBAR_CENTER(WHITEVAL), returns a new colormap with a gap
      % in the colormap around WHITEVAL value in the colormap currently in use.
      % This is useful to set the predominant colorlevel of a plot to white.
      %
      %   PLOT_COLORBAR_CENTER(WHITEVAL,OLD) uses OLD colormap instead of the
      %   default one.
      %
      %   PLOT_COLORBAR_CENTER(WHITEVAL,OLD,GAP) additionally allows the size of
      %   the gap to be is controled. GAP*100 is the number of color shades to be
      %   whitened, in each direction. If <NELEM> is 100, then there will be 20
      %   white-ish color shades. Values > 0.5 do not make sense.

      % Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>

      if ~exist('WHITEVAL','var') || isempty(WHITEVAL)
          WHITEVAL=0;
      end
      if ~exist('old','var') || isempty(old)
          old=colormap;
      end

      %getting min/max values of current colorscale
      tmp=caxis;
      MINVAL=tmp(1);
      MAXVAL=tmp(2);
      NELEM=size(old,1);

      %old colorbar domain
      x_old=linspace(MINVAL,MAXVAL,NELEM);

      %% parameters

      %this is the length of the gap around white. This gap is produced by
      %whitening the colors around <WHITEVAL>. <gap*100> is the number of color
      %shades to be whitened, in each direction. If <NELEM> is 100, then there
      %will be 20 white-ish color shades. Values > 0.5 do not make sense.
      if ~exist('gap','var') || isempty(gap)
          gap=0.1;
      end

      %% adding white zero

      %number of gap indexes
      gap_idx_nr=ceil(gap*NELEM);

      %searching for colormap index that is closes to zero
      %N1_old=find(abs(x_old)==min(abs(x_old)),1,'first'); %This puts white on zero
                                                           % without stretching to colormap
      N1_old=floor(length(x_old)/2); %this puts white on zero while stretching the
                                     %colormap, so that it's center is shifted to
                                     %zero
      N2_old=N1_old+1;

      %finding new center for white
      N1=floor((WHITEVAL-MINVAL)/(MAXVAL-MINVAL)*NELEM)-1;
      N2=N1+3;

      gap_idx_nr=min([2*gap_idx_nr,N1_old-1,NELEM-N2_old+1,N1-1,NELEM-N2+1])-gap_idx_nr;

      %  if (1 >= N1_old-gap_idx_nr) || (N2_old+gap_idx_nr-1 >= NELEM) || ...
      %     (1 >= N1-gap_idx_nr)     || (N2+gap_idx_nr-1 >= NELEM)
      if gap_idx_nr < 0
          disp([mfilename,':WARNING: requested white zero is out range. Not changing colormap.'])
          if nargout == 0
              %setting colormap
              colormap(old)
          else
              new=old;
          end
          return
      end

      % if (1 >= N1_old-gap_idx_nr)
      %     disp([mfilename,':WARNING: requested zero is out range.'])
      %     N1_old=2+gap_idx_nr;
      % end
      % if (N2_old+gap_idx_nr-1 >= NELEM)
      %     disp([mfilename,':WARNING: requested zero is out range.'])
      %     N2_old=NELEM-gap_idx_nr;
      % end
      % if (1 >= N1-gap_idx_nr)
      %     disp([mfilename,':WARNING: requested zero is out range.'])
      %     N1=2+gap_idx_nr;
      % end
      % if (N2+gap_idx_nr-1 >= NELEM)
      %     disp([mfilename,':WARNING: requested zero is out range.'])
      %     N2=NELEM-gap_idx_nr;
      % end

      if N2 < N1
          N2=N1+2;
      end

      %building index vectors
      idx_neg_old=1:(N1_old-gap_idx_nr+1);
      idx_neg_rescaled=linspace(1,(N1-gap_idx_nr+1),length(idx_neg_old));
      idx_neg=1:(N1-gap_idx_nr+1);

      idx_pos_old=(N2_old+gap_idx_nr):NELEM;
      idx_pos_rescaled=linspace((N2+gap_idx_nr),NELEM,length(idx_pos_old));
      idx_pos=(N2+gap_idx_nr):NELEM;

      %performing re-centering of colormap into requested zero
      new=ones(NELEM,3);
      for i=1:3
          if ~isempty(idx_neg_old) && ~isempty(idx_neg)
              %negative range
              new(idx_neg,i)=interp1(idx_neg_rescaled,old(idx_neg_old,i),...
                  idx_neg,'linear','extrap');
          end
          if ~isempty(idx_pos_old) && ~isempty(idx_pos)
              %positive range
              new(idx_pos,i)=interp1(idx_pos_rescaled,old(idx_pos_old,i),...
                  idx_pos,'linear','extrap');
          end
      end

      %% whitening the gap
      for i=1:3
          if max(idx_neg)+1<=N1
             new(max(idx_neg)+1:N1,i)=interp1([max(idx_neg),N1+1],[new(max(idx_neg),i), 1],...
                 max(idx_neg)+1:N1,'linear');
          end
          if N2<=min(idx_pos)-1
             new(N2:min(idx_pos)-1,i)=interp1([N2-1,min(idx_pos)],[1, new(min(idx_pos),i)],...
                 N2:min(idx_pos)-1,'linear');
          end
      end

      %% filtering singularities

      new(new>1)=1;
      new(new<0)=0;

      %% outputs

      if nargout == 0
          %setting colormap
          colormap(new)
          %cleaning up
          clear new
      end
    end
    function resampled=opt(old,axis_h)
      % PLOT_COLORBAR_OPT(OLD,AXIS_H) compresses the color range of the <old>
      % colormap around the data regions where there is highest data density. If
      % ignored, the input <old> defaults to the current colormap.
      %
      %   If the standard deviation of the data is too small compared to the data
      %   range, then a less concentrated colormap is generated. The trigger to
      %   this alternative method is controled by the internal variable
      %   <min_std_f>.

      % Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>
      % List of changes:
      %
      %   P.I <P.M.G.Inacio@tudelft.nl>, when no colormap is given, then always
      %   use the current one.

      if ~exist('old','var') || isempty(old)
          % P.Inacio - if no colormap is given, use the current one.
      %     if nargout == 0
      %         old=jet;
      %     else
              old=colormap;
      %     end
      else
          if size(old,1) < 3
              error([mfilename,': input colormap is too small.'])
          end
      end

      if ~exist('axis_h','var') || isempty(axis_h)
          axis_h=gca;
      end

      %getting min/max values of current colorscale
      tmp=caxis;
      MINVAL=tmp(1);
      MAXVAL=tmp(2);
      NELEM=size(old,1);

      %% parameters

      %minimum std(data)/(MINVAL-MAXVAL) factor
      min_std_f=0.05;

      %% building colormap that shows variability best

      %getting data
      data=[];
      h=get(gca,'Children');
      for i=1:length(h)
         switch get(h(i),'Type')
             case {'image','surface'}
                 disp(['Found data of ',get(h(i),'Type'),'.'])
                 data=get(h(i),'CData');
                 break
         end
      end
      %bug trap
      if isempty(data)
          error([mfilename,': could not find useful data in the current plot to make colormap.'])
      end
      %filtering out NaNs
      data=data(isfinite(data(:)));

      %determining old domain
      x_old=linspace(MINVAL,MAXVAL,NELEM);

      %avoiding automatic procedure for data with extremely low std relative to
      %the plot domain
      if std(data) < min_std_f*abs(mean(data))
          %extremly low data std detected, artificially creating histogram
          h=pdf('Normal',x_old,mean(data),min_std_f*(MAXVAL-MINVAL));
          disp([mfilename,':warning: data std is very low (',num2str(std(data)),...
              ') so generating histogram considering std of ',num2str(min_std_f*(MAXVAL-MINVAL))])
      else
          %getting histogram of the data
          [h,x_h]=hist(data,NELEM);
          %need the histogram in the old x domain
          h=interp1(x_h,h,x_old,'linear','extrap');
      end

      %weights are non-dimensional histogram
      w=h/sum(h);
      %removing zeros
      w(w==0)=[];
      %calculating interval size
      a=w*(MAXVAL-MINVAL);
      %calculating new domain
      x_new=cumsum(a)+MINVAL;

      %resampling
      resampled=zeros(length(a),3);
      for i=1:3
          %negative range
          resampled(:,i)=interp1(x_old,old(:,i),x_new,'linear','extrap');
      end
      resampled(resampled<0)=0;
      resampled(resampled>1)=1;

      if nargout == 0
          %setting colormap
          colormap(axis_h,resampled)
          %cleaning up
          clear resampled
      end
    end
  end
end