classdef gravity
  %static
  properties(Constant,GetAccess=private)
    %default value of some internal parameters
    default_list=struct(...
      'G',        6.67408e-11,...      % Gravitational constant [m3/kg/s2]
      'GM',       398600.4415e9,...    % Standard gravitational parameter [m^3 s^-2]
      'R',        6378136.460,...      % Earth's equatorial radius [m]
      'Rm',       6371000,...          % Earth's mean radius [m]
      'rho_earth',5514.3231,...        % average density of the Earth = (GM/G) / (4/3*pi*R_av^3) [kg/m3]
      'rho_water',1000,...             % water density
      'love',  [  0       0.000;...    % Love numbers
                  1       0.027;...
                  2      -0.303;...
                  3      -0.194;...
                  4      -0.132;...
                  5      -0.104;...
                  6      -0.089;...
                  7      -0.081;...
                  8      -0.076;...
                  9      -0.072;...
                  10     -0.069;...
                  12     -0.064;...
                  15     -0.058;...
                  20     -0.051;...
                  30     -0.040;...
                  40     -0.033;...
                  50     -0.027;...
                  70     -0.020;...
                  100    -0.014;...
                  150    -0.010;...
                  200    -0.007]...
    );
  end
  %read only
  properties(SetAccess=private)
    GM
    R
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    
  end
  %calculated only when asked for
  properties(Dependent)
    
  end
  methods(Static)
    function out=issquare(in)
      in_sqrt=sqrt(in);
      out=(in_sqrt-round(in_sqrt)==0);
    end
    function out=default
      out=gravity.default_list;
    end
    function out=y_valid(y)
      out=min(size(y))==1 && gravity.issquare(numel(y));
    end
    function out=mat_valid(mat)
      out=size(mat,1)==size(mat,2);
    end
    function out=cs_valid(in)
      out=isfield(in,'C') && isfield(in,'S');
      if ~out,return,end
      z=triu(in.C,1)+triu(in.S,0);
      out=all(size(in.C)==size(in.S)) && size(in.C,1) == size(in.C,2) && all(z(:)==0) ;
    end
    function out=mod_valid(mod)
      n=max(mod(:,1))+1;
      out=size(mod,2)==4 && size(mod,1)==n*(n+1)/2;
    end
    function mat=y2mat(y)
      mat=zeros(sqrt(numel(y)));
      mat(:)=y;
    end
    function y=mat2y(mat)
      y=mat(:);
    end
    function out=mat2cs(mat)
      out=struct(...
        'C',tril(mat,0),...
        'S',transpose(triu(mat,1))...
      );
    end
    function mat=cs2mat(in)
      mat=in.C+transpose(in.S);
    end
    function mod=cs2mod(in)
      %shortcuts
      n=size(in.C,1);
      %create lower triangle index matrix
      idxm=zeros(n);
      idxm(:)=1:n*n;
      idxm(idxm==triu(idxm,1))=NaN;
      %create index list
      [d,o]=ind2sub(n,idxm(:));
      %get location of NaNs
      i=isnan(d);
      %flatten coefficients
      C=in.C(:);S=in.S(:);
      %filter out upper diagonals
      d(i)=[];
      o(i)=[];
      C(i)=[];
      S(i)=[];
      %assemble
      mod=[d-1,o-1,C,S];
    end
    function out=mod2cs(mod)
      %shortcuts
      n=max(mod(:,1))+1;
      %create lower triangle index matrix
      idxm=zeros(n);
      idxm(:)=1:n*n;
      %make room
      out=struct(...
        'C',zeros(max(mod(:,1))+1),...
        'S',zeros(max(mod(:,1))+1)...
      );
      %propagate
      out.C(idxm==tril(idxm, 0)) = mod(:,3);
      out.S(idxm==tril(idxm, 0)) = mod(:,4);
      %assing outputs
    end
    %data type converter
    function out=dtc(from,to,in)
      %trivial call
      if strcmpi(from,to)
        out=in;
        return
      end
      %check input
      switch lower(from)
      case 'y';  c=gravity.y_valid(in);
      case 'mat';c=gravity.mat_valid(in);
      case 'cs'; c=gravity.cs_valid(in);
      case 'mod';c=gravity.mod_valid(in);
      otherwise
        error([mfilename,': unknown data type ''',from,'''.'])
      end
      if ~c
        error([mfilename,': invalid data of type ''',from,'''.'])
      end
      %convert to required types
      switch lower(from)
        case 'y'
          switch lower(to)
            case 'mat'; out=gravity.y2mat(in);
            case 'cs';  out=gravity.mat2cs(gravity.y2mat(in));
            case 'mod'; out=gravity.cs2mod(gravity.mat2cs(gravity.y2mat(in)));
          end
        case 'mat'
          switch lower(to)
            case 'y';   out=gravity.mat2y(in);
            case 'cs';  out=gravity.mat2cs(in);
            case 'mod'; out=gravity.cs2mod(gravity.mat2cs(in));
          end
        case 'cs'
          switch lower(to)
            case 'y';   out=gravity.mat2y(gravity.cs2mat(in));
            case 'mat'; out=gravity.cs2mat(in);
            case 'mod'; out=gravity.cs2mod(in);
          end
        case 'mod'
          switch lower(to)
            case 'y';   out=gravity.mat2y(gravity.cs2mat(gravity.mod2cs(in)));
            case 'mat'; out=gravity.cs2mat(gravity.mod2cs(in));
            case 'cs';  out=gravity.mod2cs(in);
          end
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
    function test(l)
      
      if ~exist('l','var') || isempty(l)
        l=20;
      end

      nr_coeff=(l+1)^2;
      
      dt={'y','mat','cs','mod'};
      y=randn(nr_coeff,1);
      dd={...
        y,...
        gravity.y2mat(y),...
        gravity.mat2cs(gravity.y2mat(y)),...
        gravity.cs2mod(gravity.mat2cs(gravity.y2mat(y)))...
      };
      for i=1:numel(dt)
        for j=1:numel(dt)
          out=gravity.dtc(dt{i},dt{j},dd{i});
          switch dt{j}
          case 'cs'
            c=any(any([out.C,out.S] ~= [dd{j}.C,dd{j}.S]));
          otherwise
            c=any(any(out~=dd{j}));
          end
          if c
            error([mfilename,': failed data type conversion between ''',dt{i},''' and ''',dt{j},'''.'])
          end
        end
      end
      
    end
  end
  methods
    %% constructor
    function obj=simplefreqseries(t,y,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't',     @(i)                 ~isscalar(i)); %this can be char, double or datetime
      p.addRequired( 'y',     @(i) simpledata.valid_y(i));
      % parse it
      p.parse(t,y,varargin{:});

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
    function obj=copy_metadata(obj,obj_in)
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in);
%TODO: fix this if there is relevant metadata in this object
%       parameters=fields(simpletimeseries.parameter_list);
%       for i=1:numel(parameters)
%         if isprop(obj,parameters{i}) && isprop(obj_in,parameters{i})
%           obj.(parameters{i})=obj_in.(parameters{i});
%         end
%       end
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
      p.addParameter('method',  simplefreqseries.default.psd_method,@(i) ischar(i));
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
        f=obj.f;
        psd=zeros(numel(f),obj.width);
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
          [psd,f]=periodogram(working.y_masked,[],2^nextpow2(obj.length),1/obj.step_num,freqrange);
        case 'fft'
          %build long filter domain
          f=obj.f;
          m=numel(f);
          %compute Fourier coefficients
          X=fft(working.y_masked,obj.psd_nfft);
          %compute power spectra
          psd=X(1:m,:).*conj(X(1:m,:))/obj.psd_nfft*obj.step_num;
          %one-sidedeness
          if p.Results.onesided
            psd=2*psd;
          end
        otherwise
          error([mfilename,': unknown method ''',method,'''.'])
        end
      end
      if any(isnan(psd(:)))
        error([mfilename,': detected NaNs in psd computation algorithm.'])
      end
      y_units=cell(size(obj.y_units));
      for i=1:numel(obj.y_units)
        y_units{i}=[obj.y_units{i},'/(Hz)^{1/2}'];
      end
      obj.psdi=simpledata(f,psd,...
        'y_units',y_units,...
        'x_units','Hz',...
        'labels' ,obj.labels,...
        'descriptor',obj.descriptor...
      );
      if p.Results.smooth
        obj=obj.smooth;
      end
      %sanitize
      obj.check
    end
    %% management methods
    function check(obj)
      %call superclass
      check@simpletimeseries(obj)
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
      obj.check
    end
    %% band-pass filtering
    function [obj,filter_response]=fft_bandpass(obj,Wn,varargin)
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
      obj=obj.assign(fx,obj.t,obj,mask);
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
      obj.check
    end
    function obj=bandpass(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      % add stuff as needed
      p.addParameter('Wn',NaN,@(i) isnumeric(i) && numel(i)==2);
      p.addParameter('method',simplefreqseries.default.bandpass_method,@(i) ischar(i));
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
      p.addRequired('cutoff',                          @(i) isnumeric(i) && numel(i)==1);
      p.addParameter('nSigma',simpledata.default.nSigma,@(i) isnumeric(i));
      p.addParameter('debug_plot',false,                @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(cutoff,varargin{:});
      % clear varargin of some parameters
      varargin=simpledata.vararginclean(varargin,{'debug_plot'});
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
      obj.check
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
