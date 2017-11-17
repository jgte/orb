classdef simplefreqseries < simpletimeseries
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list={...
      'psd_method',     'periodogram',@(i) ischar(i);...
      'bandpass_method','fft',        @(i) ischar(i);...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={};
  end
  %NOTE: edit this if you add a new parameter (if read only)
  properties(SetAccess=private)
    nyquist
    bandpass_method
    psd_method
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    psdi
  end
  %calculated only when asked for
  properties(Dependent)
    psd
    f
  end
  methods(Static)
    function out=parameters(i,method)
      persistent v parameter_names
      if isempty(v)
        v=varargs(simplefreqseries.parameter_list);
        parameter_names=v.Parameters;
      end
      if ~exist('i','var') || isempty(i)
        if ~exist('method','var') || isempty(method)
          out=parameter_names(:);
        else
          switch method
          case 'obj'
            out=v;
          otherwise
            out=v.(method);
          end
        end
      else
        if ~exist('method','var') || isempty(method)
          method='name';
        end
        if strcmp(method,'name') && isnumeric(i)
          out=parameter_names{i};
        else
          switch method
          case varargs.template_fields
            out=v.get(i).(method);
          otherwise
            out=v.(method);
          end
        end
      end
    end
    function out=transmute(in)
      if isa(in,'simplefreqseries')
        %trivial call
        out=in;
      else
        %transmute into this object
        if isprop(in,'t')
          out=simplefreqseries(in.t,in.y,in.metadata{:});
        elseif isprop(in,'x')
          out=simplefreqseries(in.x,in.y,in.metadata{:});
        else
          error('Cannot find ''t'' or ''x''. Cannot continue.')
        end
      end
    end
    function idx=get_freq_idx(ff,Wn,mode)
      switch mode
      case 'lower'
        if Wn <= 0 || max(ff) < Wn
          idx=0;
        else
          idx=find(ff<Wn,1,'last');
        end
      case 'upper'
        if inf<=Wn || min(ff(ff>0))<Wn
          idx=0;
        else
          idx=find(ff>Wn,1,'first');
        end
      otherwise
        error([mfilename,': unknown mode ''',mode,'''.'])
      end
      if isempty(idx)
        error([mfilename,': empty output ''idx'', debug needed!'])
      end
    end
    %general test for the current object
    function out=test_parameters(field,l,w)
      switch field
      case 'y-sin'
        %here, <l> is assumed to refer to the time domain (in integers from 1 to l, or thereabout)
        t=l';
        %reassign traditional meaning of <l>
        l=numel(t);
        % number of frequencies in the signal
        nf=5;
        %build frequency content of signal
        f=0.1*l./logspace(0,1,nf);
        %make room for output
        out=zeros(l,w);
        %build signal
        for i=1:nf
          out=out+ones(l,1)*randn(1,w)*10; %random bias
          out=out+sin(2*pi*t*f(i))*(1+10*rand(1,w)); %sinusoidal signal (random amplitudes)
          out=out+randn(l,w)*0.1; %noise
        end
      case 'Wn'
        out=[0,4e-7];
      otherwise
        out=simpledata.test_parameters(field,l,w);
      end
    end
    function test(l,w)
      
      if ~exist('l','var') || isempty(l)
        l=1000;
      end
      if ~exist('w','var') || isempty(w)
        w=3;
      end
      
      %test current object
      args=simplefreqseries.test_parameters('args',l,w);
      now=juliandate(datetime('now'),'modifiedjuliandate');
      
      t=[1:l/2,round(l/2+l/10):l];
      a=simplefreqseries(...
        now+t,...  % t
        simplefreqseries.test_parameters('y-sin',t,w),...           % y
        'no-mask',simplefreqseries.test_parameters('mask',l,w),...
        args{:},...
        'format','modifiedjuliandate'...
      );
    
      %TODO: fft_bandpass doesn't work if there are explicit gaps in the time domain
      %(which is something uncommon anyway)
      a=a.fill;

      %despiking
      
      [d,s]=a.despike(5e-7,'nSigma',1);
      figure
      subplot(2,1,1)
      a.plot('columns',1)
      d.plot('columns',1)
      s.plot('columns',1)
      legend('original','despiked','spikes')
      
      subplot(2,1,2)
      a.plot_psd('columns',1)
      d.plot_psd('columns',1)
      s.plot_psd('columns',1)
      legend('original','despiked','spikes')
      title('despiking')
      return
      %PSD-smoothing
      
      figure
      a.psd_refresh('smooth',false).plot_psd('columns',1)
      a.psd_refresh('smooth',true ).plot_psd('columns',1)
      legend('original PSD','smoothed PSD')
      title('PSD smoothing')
     
      %different PSD computation methods

      methods={'periodogram','fft'};
      figure
      for i=1:numel(methods)
        a=a.psd_refresh('method',methods{i}); hold on
        a.plot_psd('columns',1)
      end
      legend(methods)
      title('different PSD computation algorithms')
      
      %band-pass filtering
      
      b=a.fft_bandpass(simplefreqseries.test_parameters('Wn',t,w));
      figure
      subplot(2,1,1)
      a.plot('columns',1);
      b.plot('columns',1);
      legend('original','filtered')
      
      subplot(2,1,2)
      a.plot_psd('columns',1)
      b.plot_psd('columns',1)
      legend('original','filtered')
      title('band-pass filtering')

    end
  end
  methods
    %% constructor
    function obj=simplefreqseries(t,y,varargin)
      % input parsing
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y',     @(i) simpledata.valid_y(i));
      %create argument object, declare and parse parameters, save them to obj
      v=varargs.wrap('parser',p,'sources',{simplefreqseries.parameters([],'obj')},'mandatory',{t,y},varargin{:});
      %call superclass
      obj=obj@simpletimeseries(t,y,varargin{:});
      % save the arguments v into this object
      obj=v.save(obj);
      %initialize internal records
      obj.psdi=[];
      %save delta frequency
      obj.nyquist=2/obj.step_num;
      %no need to sanitize, since psdi is still empty
    end
    function obj=assign(obj,y,varargin)
      %pass it upstream
      obj=assign@simpletimeseries(obj,y,varargin{:});
      %update local records
      obj=obj.psd_init;
      obj.nyquist=2/obj.step_num;
    end
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in,[simplefreqseries.parameters;more_parameters(:)]);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      out=metadata@simpletimeseries(obj,[simplefreqseries.parameters;more_parameters(:)]);
    end
    %% psd methods
    function out=get.f(obj)
      out=1/obj.step_num/2*linspace(0,1,obj.psd_nfft/2);
    end
    function out=get.psd(obj)
      obj=psd_refresh_if_empty(obj);
      out=obj.psdi;
    end
    function obj=set.psd(obj,in)
      % NOTICE: no sanity is performed on this method: if used carelessly,
      % the PSD might not be in agreement with the timeseries, so only used
      % it when you are sure the agreement is kept, e.g. when operating the
      % PSD coefficients.
      obj.psdi=in;
      %sanitize (in very very general terms)
      obj.check_sf
    end
    function out=plot_psd(obj,varargin)
      obj=psd_refresh_if_empty(obj);
      out=obj.psd.plot(varargin{:});
      set(gca,'xscale','log','yscale','log')
      %outputs
      if nargout == 0
        clear out
      end
    end
    function obj=psd_init(obj)
      obj.psdi=[];
    end
    %% delayed constructor method
    function obj=psd_refresh_if_empty(obj)
      if isempty(obj.psdi)
        obj=psd_refresh(obj);
      end
    end
    function n=psd_nfft(obj)
      n = 2^nextpow2(obj.length);
    end
    function obj=psd_refresh(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addParameter('method',  obj.psd_method,@(i) ischar(i));
      p.addParameter('resample',false,@(i)islogical(i) && isscalar(i));
      p.addParameter('detrend', true, @(i)islogical(i) && isscalar(i));
      p.addParameter('despike', false,@(i)islogical(i) && isscalar(i));
      p.addParameter('onesided',true, @(i)islogical(i) && isscalar(i));
      p.addParameter('smooth',  true, @(i)islogical(i) && isscalar(i));
      % parse it
      p.parse(varargin{:});
      %duplicate
      working=obj;
      %interpolate gaps for computing the PSD if requested
      if p.Results.detrend
        working=working.resample;
      end
      %remove trend if requested
      if p.Results.detrend
        working=working.detrend;
      end
      %remove spikes if requested
      if p.Results.despike
        %not implemented yet!
        working=working.despike;
      end
      %handle empty data
      if isempty(working.y_masked)
        f0=obj.f;
        psd0=zeros(numel(f0),obj.width);
      else
        % compute PSD, according to requested method (y_masked is needed in
        % order to avoid NaNs going into PSD computation methods)
        switch p.Results.method
        case 'periodogram'
          if p.Results.onesided
            freqrange='onesided';
          else
            freqrange='twosided';
          end
          [psd0,f0]=periodogram(working.y_masked,[],2^nextpow2(obj.length),1/obj.step_num,freqrange);
        case 'fft'
          %build long filter domain
          f0=obj.f;
          m=numel(f0);
          %compute Fourier coefficients
          X=fft(working.y_masked,obj.psd_nfft);
          %compute power spectra
          psd0=X(1:m,:).*conj(X(1:m,:))/obj.psd_nfft*obj.step_num;
          %one-sidedeness
          if p.Results.onesided
            psd0=2*psd0;
          end
        otherwise
          error([mfilename,': unknown method ''',method,'''.'])
        end
      end
      if any(isnan(psd0(:)))
        error([mfilename,': detected NaNs in psd computation algorithm.'])
      end
      y_units=cell(size(obj.y_units));
      for i=1:numel(obj.y_units)
        y_units{i}=[obj.y_units{i},'/(Hz)^{1/2}'];
      end
      obj.psdi=simpledata(f0,psd0,...
        'y_units',y_units,...
        'x_units','Hz',...
        'labels' ,obj.labels,...
        'descriptor',obj.descriptor...
      );
      if p.Results.smooth
        obj=obj.smooth;
      end
      %sanitize
      obj.check_sf
    end
    %% management methods
    function check_sf(obj)
      %check for negative PSD entries
      if ~isempty(obj.psdi)
        if any(obj.psdi.x<0)
          error([mfilename,': the frequency domain must be positive.'])
        end
        if any(obj.psdi.y<0)
          error([mfilename,': the power must be positive.'])
        end
      end
    end
    %% PSD operations
    function obj=smooth(obj,varargin)
      obj=psd_refresh_if_empty(obj);
      dws=max([round(obj.length*1e-3),2]);
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addParameter('window_size',    dws,@(i)isnumeric(i) && isscalar(i) && round(i)==i);
      p.addParameter('iter',           2,  @(i)isnumeric(i) && isscalar(i) && round(i)==i);
      p.addParameter('window_stretch', 20, @(i)isnumeric(i) && isscalar(i) && round(i)==i);
      % parse it
      p.parse(varargin{:});
      % trivial call
      if p.Results.window_size < 2 || p.Results.iter < 1
          %nothing to do
          return
      end
      % add additional dependent parameters
      p.KeepUnmatched=false;
      p.addParameter('min_window_size',max([1,round(p.Results.window_size/100)]),@(i)isnumeric(i) && isscalar(i) && round(i)==i);
      % parse it again
      p.parse(varargin{:});
      %easier names
      n=obj.psd.length;
      %sanity
      if p.Results.window_size > n
          disp([mfilename,':WARNING: window size too big for input data (',...
              num2str(p.Results.window_size),'), setting at maximum value (',num2str(n),').'])
          p.Results.window_size=n;
      end
      %building window size domain
      if (p.Results.window_stretch==0)
          wz=ones(1,n)*p.Results.window_size;
      else
          min_win_size=max([1,round(p.Results.window_size/100)]);
          wz=round(min_win_size+(p.Results.window_size-min_win_size)*((1:n)/n).^p.Results.window_stretch);
      end
      %building lower/upper index
      lower_idx=max([  ones(1,n);(1:n)-wz]);
      upper_idx=min([n*ones(1,n);(1:n)+wz]);
      %propagating
      out=obj.psd.y;
      %smoothing
      for j=1:p.Results.iter
        for i = 1:n
          out(i,:)=sum(out( lower_idx(i):upper_idx(i),: ),1)/(upper_idx(i)-lower_idx(i)+1);
        end
      end
      %recovering
      obj.psd=obj.psd.assign(out);
      %sanitize
      obj.check_sf
    end
    %% band-pass filtering
    function [obj,filter_response]=fft_bandpass(obj,Wn,varargin)
      obj=psd_refresh_if_empty(obj);
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addRequired('Wn',                  @(i) isnumeric(i) && numel(i)==2);
      p.addParameter('gaps',      'zeroed',@(i) ischar(i));
      p.addParameter('detrend',   true,    @(i) islogical(i) && isscalar(i));
      p.addParameter('debug_plot',false,   @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(Wn,varargin{:});
      %handle gaps
      switch p.Results.gaps
      case 'trunc'
        %truncating bad data (not a good idea)
        data_in=obj.y_masked;
      case 'zeroed'
        %zeroing bad data
        data_in=obj.y;
        data_in(~obj.mask,:)=0;
      otherwise
        error([mfilename,': unknown gap handling mode ''',p.Results.gaps,'''.'])
      end
      %sanity
      if any(isnan(data_in(:)))
          error([mfilename,': found NaNs in the input data.'])
      end
      disp(['FFT filter: [',num2str(Wn),']'])
      %computational length
      n = 2^nextpow2(size(data_in,1));
      %build long filter domain
      ff=1/obj.step_num/2*linspace(0,1,n/2);
      fP=zeros(size(ff));
      %assign filter factors
      fP(Wn(1)<=ff&ff<=Wn(2))=1;
      fP(Wn(1)<=ff&ff<=Wn(2))=1;
      %parameters (min is needed in case the pass-band is very wide)
      smooth_radius=min([ceil(sum(~fP)*0.1),ceil(sum(fP)*0.1)]); %data points
      %smooth transitions
      idx=[...
        simplefreqseries.get_freq_idx(ff,Wn(1),'lower'),...
        simplefreqseries.get_freq_idx(ff,Wn(2),'upper')...
      ];
      % figure
      % semilogx(ff,fP), hold on
      for i=1:2
        if idx(i)>0
          idx_out=(idx(i)-smooth_radius):(idx(i)+smooth_radius+1);
          idx_in =[idx_out(1),idx_out(end)];
          fP(idx_out)=spline(ff(idx_in),[0 fP(idx_in) 0],ff(idx_out));
        end
      end
      % semilogx(ff,fP), hold on
      % keyboard
      %mirror the filter
      fP=[fP,fliplr(fP)];
      %apply the filter
      fX=fft(data_in,n).*(fP(:)*ones(1,size(data_in,2)));
      fx=ifft(fX,'symmetric');
      %trim excess
      fx=fx(1:size(data_in,1),:);
      if p.Results.debug_plot
        m=numel(ff);
        X=fft(data_in(:,1),n);
        PX=X(1:m).*conj(X(1:m));
        PfX=fX(1:m,1).*conj(fX(1:m,1));
        figure
        subplot(2,1,1)
        title('frequency domain')
        loglog(ff,PX), hold on
        loglog(ff,PfX)
        loglog(ff,fP(1:m)*max([max(PX),max(PfX)]))
        legend('original','filtered','filter')

        subplot(2,1,2)
        title('time domain')
        plot(fx(:,1)), hold on
        plot(obj.y(:,1))
        legend('filtered','original')
        keyboard
      end
      %propagate
      obj=obj.assign(fx,'t',obj.t,'mask',obj.mask);
      %recompute PSD (sanitization done in psd_refresh)
      obj=obj.psd_refresh;
      %additional outputs
      if nargout>1
        filter_response.f=ff;
        filter_response.a=fP(1:numel(ff));
      end
    end
    function obj=butter_bandpass(obj,Wn,varargin)
      obj=psd_refresh_if_empty(obj);
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addRequired('Wn',                  @(i) isnumeric(i) && numel(i)==2);
      p.addParameter('gaps',      'interp',@(i) ischar(i));
      p.addParameter('detrend',   true,    @(i) islogical(i) && isscalar(i));
      p.addParameter('debug_plot',false,   @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(Wn,varargin{:});
      %handle gaps
      switch p.Results.gaps
      case 'trunc'
        %truncating bad data (not a good idea)
        data_in=obj.y_masked;
      case 'zeroed'
        %zeroing bad data
        data_in=obj.y;
        data_in(obj.mask,:)=0;
      case 'interp'
        %interpolating over bad data
        data_in=obj.resample.y;
      otherwise
        error([mfilename,': unknown gap handling mode ''',p.Results.gaps,'''.'])
      end
      
      %sanity
      if any(isnan(data_in(:)))
          error([mfilename,': found NaNs in the input data.'])
      end
      
      disp(['Butterworth filter: [',num2str(Wn),']'])

      %computational length
      n = 2^nextpow2(size(data_in,1));

      error([mfilename,': not yet implemented'])
      
      %sanitize
      obj.check_sf
    end
    function obj=bandpass(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addParameter('Wn',NaN,@(i) isnumeric(i) && numel(i)==2);
      p.addParameter('method',obj.bandpass_method,@(i) ischar(i));
      % parse it
      p.parse(varargin{:});
      %sanity
      if ~isfinite(p.Results.Wn)
        error([mfilename,': invalid value for input ''Wn'':',num2str(p.Results.Wn),'.'])
      end
      % branching
      switch p.Results.method
      case 'fft'
        obj=fft_bandpass(obj,p.Results.Wn,varargin{:});
      case 'butter'
        obj=fft_bandpass(obj,p.Results.Wn,varargin{:});
      otherwise
        error([mfilename,': unknown bandpass method ''',method,'''.'])
      end
    end
    %% time-domain operations, done at the level of the frequency domain
    function [despiked,spikes]=despike(obj,cutoff,varargin)
      obj=psd_refresh_if_empty(obj);
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addRequired('cutoff',                           @(i) isnumeric(i) && numel(i)==1);
      p.addParameter('outlier_sigma',obj.outlier_sigma, @(i) isnumeric(i));
      p.addParameter('debug_plot',false,                @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(cutoff,varargin{:});
      % clear varargin of some parameters
      varargin=cells.vararginclean(varargin,{'debug_plot'});
      % apply low-pass filter
      smoothed=obj.bandpass('Wn',[0,cutoff],varargin{:});
      % get high-frequency signal
      highfreq=obj-smoothed;
      % find outliers in high-frequency signal
      [despiked,spikes]=highfreq.outlier(p.Results.nSigma);
      spikes=spikes.psd_init;
      % restore smoothed signal
      despiked=despiked+smoothed;
      % debug plot
      if p.Results.debug_plot
        figure
        subplot(3,1,1)
             obj.plot('column',1,'line',{'-o'})
        despiked.plot('column',1,'line',{'+'})
        spikes.plus(smoothed).plot('column',1,'line',{'x'})
        smoothed.plot('column',1)
        legend('original','despiked','spikes','smooth')
        subplot(3,1,2)
        a=despiked+spikes-obj;
        a.plot('column',1)
        title('residual')
        subplot(3,1,3)
             obj.plot_psd('column',1,'line',{'-o'})
        despiked.plot_psd('column',1,'line',{'+'})
        spikes.plus(smoothed).plot_psd('column',1,'line',{'x'})
        smoothed.plot_psd('column',1)
        legend('original','despiked','spikes','smooth')
        keyboard
      end
      %sanitize
      obj.check_sf
    end
    %% overloading
    % uses a method from a superclass and resets PSD
    function obj=op(obj,operation,varargin)
      %operate
      obj=obj.(operation)(varargin{:});
      %reset PSD
      obj=obj.psd_init;
    end
    function obj=scale(obj,scale)
      obj=scale@simpledata(obj,scale);
      obj=obj.psd_init;
    end
    function obj=plus(obj,obj_new)
      obj=plus@simpledata(obj,obj_new);
      obj=obj.psd_init;
    end
    function obj=minus(obj,obj_new)
      obj=minus@simpledata(obj,obj_new);
      obj=obj.psd_init;
    end
    function obj=times(obj,obj_new)
      obj=times@simpledata(obj,obj_new);
      obj=obj.psd_init;
    end
    function obj=rdivide(obj,obj_new)
      obj=rdivide@simpledata(obj,obj_new);
      obj=obj.psd_init;
    end
    function obj=interp(obj,t_now,varargin)
      % call superclass
      obj=interp@simpletimeseries(obj,t_now,varargin{:});
      %initialize internal records
      obj.psdi=[];
      %save delta frequency
      obj.nyquist=2/obj.step_num;
    end
  end
end