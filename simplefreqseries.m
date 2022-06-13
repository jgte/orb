classdef simplefreqseries < simpletimeseries
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list={...
      'psd_method',     'periodogram',@ischar;...
      'bandpass_method','fft',        @ischar;...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={};
  end
  properties
    bandpass_method
    psd_method
  end
  %NOTE: edit this if you add a new parameter (if read only)
  properties(SetAccess=private)
    nyquist
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
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(simplefreqseries.parameter_list); end
      out=v.picker(varargin{:});
    end
    function out=transmute(in)
      if isa(in,'simplefreqseries')
        %trivial call
        out=in;
      else
        %transmute into this object
        if isprop(in,'t')
          out=simplefreqseries(in.t,in.y,in.varargin{:});
        elseif isprop(in,'x')
          out=simplefreqseries(in.x,in.y,in.varargin{:});
        else
          error('Cannot find ''t'' or ''x''. Cannot continue.')
        end
      end
    end
    %% frequency utilities
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
        error(['unknown mode ''',mode,'''.'])
      end
      if isempty(idx)
        error('empty output ''idx'', debug needed!')
      end
    end
    function out=f_domain(n,step,zero_f)
      if ~exist('zero_f','var') || isempty(zero_f)
        zero_f=true;
      end
      switch class(step)
        case 'duration'
          step=simpletimeseries.timescale(step);
        case 'double'
          %do nothing
        otherwise
          error(['Cannot handle input ''step'' of class ''',class(step),'''.'])
      end
      if zero_f
        %output frequency domain contains the zero frequency
        out=1/step/2*linspace(0,1,n);
      else
        %output frequency domain omits the zero frequency
        out=1/step/2*linspace(0,1,n+1);
        out(1)=[];
      end
    end
    function W=fourier2psd(X,n,dt)
      W=X.*conj(X)/n*dt;
    end
    %% constructors
    function obj=colored_noise(fn,Wn,tn,varargin)
      %Applies the weights defined in frequencies fn with values Wn to the time series defined in this object
      %NOTICE: this method operates on the time series (member 'y'' defined in the simpledata parent class)
      p=machinery.inputParser;
      % add stuff as needed
      p.addRequired( 'fn', @(i) (isnumeric(i) && isvector(i)) || ischar(i));
      p.addRequired( 'Wn', @(i) isnumeric(i) && isvector(i) && numel(i) == numel(fn) || isempty(i));
      p.addRequired( 'tn', @(i) simpletimeseries.valid_t(i));
      p.addParameter('seed','default',@(i)...   % does not support the form RNG(SD,GENERATOR)
        ischar(i) || ...                        % can be 'shuffle' or 'default'
        ( isnumeric(i) && isscalar(i) ) || ...  % can be the seed
        ( isstruct(i) && structs.iseq_field_list(i,struct('Type','','Seed',[],'State',[])) )... % can be a RNG structure
      );
      p.addParameter('columns',1,@(i)isnumeric(i) && isscalar(i));
      p.addParameter('plot_column',0,@(i)isnumeric(i) && isscalar(i));
      p.addParameter('Wn_interp','loglog',@ischar);
      % parse it
      p.parse(fn,Wn,tn,varargin{:});
      %init RNG
      rng(p.Results.seed);
      %branch on the type of data represented by 'fn'
      if ischar(fn)
        assert(file.exist(fn),['Cannot find file with spectrum: ',fn])
        %NOTICE: expecting data to be in two columns: fn, Wn; '#' are comments
        fnWn=file.textscan(fn,'%f %f',[],'#');
        %propagate data
        fn=fnWn(:,1);
        Wn=fnWn(:,2);
      end
      %trivial call: empty data, output is the same as input
      if all(Wn==0)
        yn=zeros(size(tn));
      else
        %resolve time step and save input time domain length
        dt=simpletimeseries.timestep(tn);
        %make some noise
        ti=tn(1):dt:tn(end);
        %make sure noise series ends after tn(end)
        if ti(end)<tn(end)
          ti(end+1)=ti(end)+dt;
        end
        ni=numel(ti);
        yi=randn(ni,p.Results.columns);
        %resolve frequency domain
        fi=simplefreqseries.f_domain(ni,dt,false);
        %interpolate fn and Wn to frequency domain of current timeseries
        switch p.Results.Wn_interp
        case 'loglog'
          Wi=10.^interp1(log10(fn),log10(Wn),log10(fi),'linear','extrap');
        case 'semilogy'
          Wi=10.^interp1(fn,log10(Wn),fi,'linear','extrap');
        case 'semilogx'
          Wi=interp1(log10(fn),Wn,log10(fi),'linear','extrap');
        case 'linear'
          Wi=interp1(fn,Wn,fi,'linear','extrap');
        otherwise
          error(['Cannot handle argumen ''Wn_interp'' with value ''',p.Results.Wn_interp,...
            ''', expecting one of ''loglog'', ''semilogy'', ''semilogx'' or ''linear''.'])
        end
        %determine fft computational length
        n=2^nextpow2(ni);
        %scaling to fourier coefficients
        Xi=zeros(n,p.Results.columns);
        Xi(1:ni,:)=complex(sqrt(Wi(:)))*ones(1,p.Results.columns);
        %compute Fourier coefficients of the noise
        X=fft(yi,n);
        %scale Fourier coefficients in the frequency domain
        Xo=X.*Xi;
        %convert back to the time domain
        yo=real(ifft(Xo,n,'symmetric'));
        %only need half of the time series
        no=n/2;
        yo=yo(1:no,:);
        to=linspace(tn(1),tn(end),no);
        %interpolate noise back to the requested time domain
        yn=interp1(to,yo,tn);
      end
      %build object
      obj=simplefreqseries(tn,yn,'descriptor','noise',varargin{:});
      %show plot, if requested
      if p.Results.plot_column>0
        c=p.Results.plot_column;
        plotting.figure;
        subplot(3,1,1)
        plot(fn,Wn*simpletimeseries.timescale(dt),'o-'), hold on
        plot(fi,Wi*simpletimeseries.timescale(dt))
        simplefreqseries(to,yo).plot_psd('columns',c,'method','fft','resample',false,'detrend',false)
        obj.plot_psd('columns',c,'method','fft','resample',false,'detrend',false)
        legend('Wn*dt','Wn*dt interp',['Noise ifft col.',num2str(c)],['Noise interp col.',num2str(c)],'location','west')
        subplot(3,1,2)
        plot(ti,yi(:,c)), hold on
        plot(to,yo(:,c),'-o')
        obj.plot('columns',c,'line',{'-+'})
        legend('white noise','colored noise ifft','colored noise interp');
        subplot(3,1,3)
        histogram(yi(:,c)); hold on
        histogram(yo(:,c));
        histogram(obj.y(:,c));
        legend('white noise','colored noise ifft','colored noise interp');
      end
    end
    %% general test for the current object
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
        out=[2e-7,1e-6];
      case 'time'
        %set the number of gaps to be 5% of the length of the time domain
        n_gaps=floor(l*0.1);
        %generate a time domain:
        % - <l> refers to the length of the data
        % - generate n_gaps more than needed, so that there's l number of data
        out=juliandate(datetime('now'),'modifiedjuliandate')+linspace(1,l,l+n_gaps);
        %insert implicit gaps
        while numel(out)>l
          out(round(2+rand*(numel(out)-3)))=[];
        end
      otherwise
        out=simpledata.test_parameters(field,l,w);
      end
    end
    function out=test(method,l,w)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      if ~exist('l','var') || isempty(l)
        l=1000;
      end
      if ~exist('w','var') || isempty(w)
        w=3;
      end
      %get common parameters
      args=simplefreqseries.test_parameters('args',l,w);
      %define time domain
      t=simplefreqseries.test_parameters('time',l,w);
      %init object
      y=simplefreqseries.test_parameters('y-sin',t,w);
      mask=simplefreqseries.test_parameters('mask',l,w);
      a=simplefreqseries(t,y,...
        'mask',mask,...
        args{:},...
        'format','modifiedjuliandate'...
      );
      switch method
        case 'all'
          for i={'init','despike','smooth','psd','band-pass','noises'}
            simplefreqseries.test(i{1},l);
          end
        case {'constructor','init'}
          out=a;
        case 'print'
          a.psd_refresh.print;
          out=a;
        case 'despike'
          [d,s]=a.despike(5e-7,'outlier_sigma',1);
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
          out=d;
        case 'smooth'
          figure
          a.psd_refresh('smooth',false).plot_psd('columns',1)
          a.psd_refresh('smooth',true ).plot_psd('columns',1)
          legend('original PSD','smoothed PSD')
          title('PSD smoothing')
          out=a;
        case 'psd'
          methods={'periodogram','fft'};
          figure
          for i=1:numel(methods)
            a=a.psd_refresh('method',methods{i}); hold on
            a.plot_psd('columns',1)
          end
          legend(methods)
          title('different PSD computation algorithms')
          out=a;
        case 'band-pass'
          Wn=simplefreqseries.test_parameters('Wn',t,w);
%           %TODO: fft_bandpass doesn't work if there are explicit gaps in the time domain
%           %(which is something uncommon anyway)
%           a=a.fill;
          b=a.fft_bandpass(Wn);
          figure
          subplot(2,1,1)
          a.plot('columns',1);
          b.plot('columns',1);
          legend('original','filtered')
          subplot(2,1,2)
          a.plot_psd('columns',1)
          b.plot_psd('columns',1)
          limits=axis;
          plot([Wn(1) Wn(1)],limits(3:4),'k:')
          plot([Wn(2) Wn(2)],limits(3:4),'k:')
          legend('original','filtered')
          title(method)
          out=b;
        case 'noises'
          for i={'step-down','step-up','pink','brown','blue','violet'}
            simplefreqseries.test(['noise-',i{1}],l);
          end
        case {'noise-step-down','noise-step-up','noise-pink','noise-brown','noise-blue','noise-violet'}
          if l<5000
            disp(['NOTICE: input ''l'' is ',num2str(l),', which is lower than the advised value of 5000'])
          end
          switch method
            case 'noise-step-down'
              fn=[1e-7,5.9e-7,6e-7,1e-5];
              Wn=[   1,     1,1e-3,1e-3];
            case 'noise-step-up'
              fn=[1e-7,5.9e-7,6e-7,1e-5];
              Wn=[1e-3,  1e-3,   1,   1];
            %https://en.wikipedia.org/wiki/Colors_of_noise
            case 'noise-pink'
              fn=[1e-7,1e-6,1e-5];
              Wn=[  10,   1, 0.1];
            case 'noise-brown'
              fn=[1e-7,1e-6,1e-5];
              Wn=[ 100,   1,0.01];
            case 'noise-blue'
              fn=[1e-7,1e-6,1e-5];
              Wn=[ 0.1,   1,  10];
            case 'noise-violet'
              fn=[1e-7,1e-6,1e-5];
              Wn=[0.01,   1, 100];
          end
          %shape that noise
          b=simplefreqseries.colored_noise(fn,Wn,a.t,'columns',w,'plot_column',randi(w,1,1));
          out=b;
      end
    end
  end
  methods
    %% constructor
    function obj=simplefreqseries(t,y,varargin)
      % input parsing
      p=machinery.inputParser;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y',     @(i) simpledata.valid_y(i));
      %create argument object, declare and parse parameters, save them to obj
      v=varargs.wrap('parser',p,'sources',{simplefreqseries.parameters('obj')},'mandatory',{t,y},varargin{:});
      %call superclass
      obj=obj@simpletimeseries(t,y,varargin{:});
      % save the arguments v into this object
      obj=v.save(obj,{'t','y'});
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
    function obj=copy_metadata(obj,obj_in,more_parameters,less_parameters)
      if ~exist('less_parameters','var')
        less_parameters={};
      end
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in,[simplefreqseries.parameters('list');more_parameters(:)],less_parameters);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      out=metadata@simpletimeseries(obj,[simplefreqseries.parameters('list');more_parameters(:)]);
    end
    %the varargin method can be called directly
    %% info methods
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      relevant_parameters={'bandpass_method','psd_method','nyquist'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print pdf, if available
      if isempty(obj.psdi)
        obj.disp_field('psd',tab,'empty')
      else
        obj.psdi.print(tab)
      end
      %print superclass
      print@simpletimeseries(obj,tab)
    end
    %% psd methods
    function out=get.f(obj)
      out=simplefreqseries.f_domain(obj.psd_nfft/2+1,obj.step_num);
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
      obj=psd_refresh_if_empty(obj,varargin{:});
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
    function obj=psd_refresh_if_empty(obj,varargin)
      if isempty(obj.psdi)
        obj=psd_refresh(obj,varargin{:});
      end
    end
    function n=psd_nfft(obj)
      n = 2^nextpow2(obj.length);
    end
    function obj=psd_refresh(obj,varargin)
      p=machinery.inputParser;
      % add stuff as needed
      p.addParameter('method',  obj.psd_method,@ischar);
      p.addParameter('resample',true, @(i)islogical(i) && isscalar(i));
      p.addParameter('detrend', true, @(i)islogical(i) && isscalar(i));
      p.addParameter('onesided',true, @(i)islogical(i) && isscalar(i));
      p.addParameter('smooth',  true, @(i)islogical(i) && isscalar(i));
      % parse it
      p.parse(varargin{:});
      %duplicate
      working=obj;
      %interpolate gaps for computing the PSD if requested
      if p.Results.resample
        working=working.resample;
      end
      %remove trend if requested
      if p.Results.detrend
        working=working.detrend;
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
          psd0=simplefreqseries.fourier2psd(X(1:m,:),obj.psd_nfft,obj.step_num);
          %one-sidedeness
          if p.Results.onesided
            psd0=2*psd0;
          end
        otherwise
          error(['unknown method ''',method,'''.'])
        end
      end
      if any(isnan(psd0(:)))
        error('detected NaNs in psd computation algorithm.')
      end
      units=cell(size(obj.units));
      for i=1:numel(obj.units)
        units{i}=[obj.units{i},'/(Hz)^{1/2}'];
      end
      obj.psdi=simpledata(f0,psd0,...
        'units',units,...
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
          error('the frequency domain must be positive.')
        end
        if any(obj.psdi.y<0)
          error('the power must be positive.')
        end
      end
    end
    %% PSD operations
    function obj=smooth(obj,varargin)
      obj=psd_refresh_if_empty(obj);
      dws=max([round(obj.length*1e-3),2]);
      p=machinery.inputParser;
      % add stuff as needed
      p.addParameter('window_size',    dws,@(i)num.isscalar(i) && round(i)==i);
      p.addParameter('iter',           2,  @(i)num.isscalar(i) && round(i)==i);
      p.addParameter('window_stretch', 20, @(i)num.isscalar(i) && round(i)==i);
      % parse it
      p.parse(varargin{:});
      % trivial call
      if p.Results.window_size < 2 || p.Results.iter < 1
          %nothing to do
          return
      end
      % add additional dependent parameters
      p.addParameter('min_window_size',max([1,round(p.Results.window_size/100)]),@(i)num.isscalar(i) && round(i)==i);
      % parse it again
      p.parse(varargin{:});
      %easier names
      n=obj.psd.length;
      %sanity
      if p.Results.window_size > n
          warning(['window size too big for input data (',...
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
      p=machinery.inputParser;
      % add stuff as needed
      p.addRequired('Wn',                  @(i) isnumeric(i) && numel(i)==2);
      p.addParameter('gaps',      'zeroed',@ischar);
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
        error(['unknown gap handling mode ''',p.Results.gaps,'''.'])
      end
      %sanity
      if any(isnan(data_in(:)))
          error('found NaNs in the input data.')
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
        idx=3;
        m=numel(ff);
        X=fft(data_in(:,idx),n);
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
        plot(fx(:,idx)), hold on
        plot(obj.y(:,idx))
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
      p=machinery.inputParser;
      % add stuff as needed
      p.addRequired('Wn',                  @(i) isnumeric(i) && numel(i)==2);
      p.addParameter('gaps',      'interp',@ischar);
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
        error(['unknown gap handling mode ''',p.Results.gaps,'''.'])
      end

      %sanity
      if any(isnan(data_in(:)))
          error('found NaNs in the input data.')
      end

      disp(['Butterworth filter: [',num2str(Wn),']'])

%       %computational length
%       n = 2^nextpow2(size(data_in,1));
%
%       error(['not yet implemented'])
%
%       %sanitize
%       obj.check_sf
    end
    function obj=bandpass(obj,varargin)
      p=machinery.inputParser;
      % add stuff as needed
      p.addParameter('Wn',NaN,@(i) isnumeric(i) && numel(i)==2);
      p.addParameter('method',obj.bandpass_method,@ischar);
      % parse it
      p.parse(varargin{:});
      %sanity
      if ~isfinite(p.Results.Wn)
        error(['invalid value for input ''Wn'':',num2str(p.Results.Wn),'.'])
      end
      % branching
      switch p.Results.method
      case 'fft'
        obj=fft_bandpass(obj,p.Results.Wn,varargin{:});
      case 'butter'
        obj=fft_bandpass(obj,p.Results.Wn,varargin{:});
      otherwise
        error(['unknown bandpass method ''',method,'''.'])
      end
    end
    %% time-domain operations, done at the level of the frequency domain
    function [despiked,spikes]=despike(obj,cutoff,varargin)
      obj=psd_refresh_if_empty(obj);
      p=machinery.inputParser;
      % add stuff as needed
      p.addRequired('cutoff',                           @(i) isnumeric(i) && numel(i)==1);
      p.addParameter('debug_plot',false,                @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(cutoff,varargin{:});
      % clear varargin of some parameters
      varargin=cells.vararginclean(varargin,{'debug_plot'});
      % apply low-pass filter
      smoothed=obj.bandpass('Wn',[0,cutoff],varargin{:});
      % get high-frequency signal
      highfreq=obj-smoothed;
      % find outliers in high-frequency signal (handles outlier_iter, outlier_sigma and detrend options)
      [despiked,spikes]=highfreq.outlier(varargin{:});
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