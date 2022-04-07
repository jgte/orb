classdef gravity < simpletimeseries
  %static
  properties(Constant)
    %this is used to define the epoch of static fields (in the gravity.load method)
   static_select_date=time.zero_date;
    static_start_date=time.zero_date;
     static_stop_date=time.inf_date;
%   Supported functionals are the following,
%       'nondim'    - non-dimensional Stoked coefficients.
%       'eqwh'      - equivalent water height [m]
%       'geoid'     - geoid height [m]
%       'potential' - [m2/s2]
%       'gravity'   - [m /s2], if input represents the disturbing potential,
%                     then the output represents the gravity disturbances.
%                     Otherwise, it represents the gravity accelerations
%                     (-dU/dr).
%       'anomalies' - If input represents the disturbing potential, then
%                     the output represents the gravity anomalies.
%                     Otherwise, it represents (-dU/dr - 2/r*U).
%       'vertgravgrad' - vertical gravity gradient.
    functional_details=struct(...
      'nondim',        struct('units',' ',         'name','Stokes Coeff.'),...
      'eqwh',          struct('units','m',         'name','Eq. Water Height'),...
      'geoid',         struct('units','m',         'name','Geoid Height'),...
      'potential',     struct('units','m^2.s^{-2}','name','Geopotential'),...
      'anomalies',     struct('units','m.s^{-2}',  'name','Gravity Anomalies'),...
      'vertgravgrad',  struct('units','s^{-2}',    'name','Vert. Gravity Gradient'),...
      'gravity',       struct('units','m.s^{-2}',  'name','Gravity Acc.')...
    );
  end
  properties(Constant,GetAccess=private)
    %TODO: the parameters below are a mix of physical constants and default values for properties; find a way to fix this
    %default value of some internal parameters
    parameter_list={...
      'GM',       398600.4415e9, @num.isscalar;...      % Standard gravitational parameter [m^3 s^-2]
      'R',        6378136.460,   @num.isscalar;...      % Earth's equatorial radius [m]
      'rho_earth',5514.32310829, @num.isscalar;...      % average density of the Earth = (GM/G) / (4/3*pi*R_av^3) [kg/m3]
      'rho_water',1000,          @num.isscalar;...      % water density [kg/m3]
      'G',        6.67408e-11,   @(i) (num.isscalar(i)) || isempty(i);...      % Gravitational constant [m3/kg/s2]
      'Rm',       6371000,       @(i) (num.isscalar(i)) || isempty(i);...      % Earth's mean radius [m]
      'love',  [  0       0.000;...               
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
                  200    -0.007], @(i) isnumeric(i) && size(i,2)==2;...     % Love numbers ?http://dx.doi.org/10.1029/98JB02844
        'origin',      'unknown', @ischar;...                        % (arbitrary string)
        'functional',   'nondim', @(i) ischar(i) && any(strcmp(i,gravity.functionals)); %see above
        'zf_love', 0.30190,       @num.isscalar;...      % zero frequency Love number: reported in IERS2003 Section 6.3 as "k20"
        'pt_factor',1.391413e-08, @num.isscalar;...      % permanent tide factor: reported in IERS2003 Section 6.3 as "A0*H0", (4.4228e-8)*(0.31460)
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'GM','R','functional','lmax'};
  end
  properties(SetAccess=public)
    GM
    R
    functional
    origin
    cdate
  end
  %calculated only when asked for
  properties(Dependent)
    lmax
    mat
    cs
    tri
    mod
    checksum
    funct %handles empty values of functional
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(gravity.parameter_list); end
      out=v.picker(varargin{:});
    end
    function out=issquare(in)
      in_sqrt=sqrt(in);
      out=(in_sqrt-round(in_sqrt)==0);
    end
    function out=functionals
      out=fieldnames(gravity.functional_details);
    end
    function out=functional_units(functional)
      assert(gravity.isfunctional(functional),['Ilegal functional ''',functional,'''.'])
      out=gravity.functional_details.(functional).units;
    end
    function out=functional_names(functional)
      assert(gravity.isfunctional(functional),['Ilegal functional ''',functional,'''.'])
      out=gravity.functional_details.(functional).name;
    end
    function out=functional_label(functional)
      out=[gravity.functional_names(functional),' [',gravity.functional_units(functional),']'];
    end
    function out=isfunctional(in)
      out=isfield(gravity.functional_details,in);
    end
    %% y representation
    function out=y_valid(y)
      out=size(y,1)==1 && gravity.issquare(numel(y));
    end
    function lmax=y_lmax(y)
      lmax=sqrt(numel(y))-1;
    end
    function s=y_length(lmax)
      s=(lmax+1).^2;
    end
    %% mat representation
    function out=mat_valid(mat)
      out=size(mat,1)==size(mat,2);
    end
    function lmax=mat_lmax(mat)
      lmax=gravity.y_lmax(mat(:));
    end
    function s=mat_length(lmax)
      s=[lmax+1,lmax+1];
    end
    %% cs representation
    function out=cs_valid(in)
      out=isfield(in,'C') && isfield(in,'S');
      if ~out,return,end
      z=triu(in.C,1)+triu([in.S(:,2:end),zeros(size(in.C,1),1)],0);
      out=all(size(in.C)==size(in.S)) && size(in.C,1) == size(in.C,2) && all(z(:)==0) ;
    end
    function lmax=cs_lmax(cs)
      lmax=gravity.y_lmax(cs(1).C(:));
      for i=2:numel(cs)
        assert(lmax==gravity.y_lmax(cs(i).C(:)),'found discrepancy in sizes in ''cs'' representation. Debug needed.');
      end
    end
    function s=cs_length(lmax)
      s=gravity.mat_length(lmax);
    end
    %% tri representation
    function out=tri_valid(tri)
      lmax=gravity.tri_lmax(tri);
      out=size(tri,1)==lmax+1 && round(lmax)==lmax && gravity.cs_valid(gravity.tri2cs(tri));
    end
    function lmax=tri_lmax(tri)
      lmax=(size(tri,2)+1)/2-1;
    end
    function s=tri_length(lmax)
      s1=gravity.mat_length(lmax);
      s=[s1(1),2*(lmax+1)-1];
    end
    %% mod representation
    function out=mod_valid(mod)
      n=gravity.mod_lmax(mod)+1;
      out=size(mod,2)==4 && size(mod,1)==n*(n+1)/2;
    end
    function lmax=mod_lmax(mod)
      lmax=max(mod(:,1));
    end
    function s=mod_length(lmax)
      s=[(lmax*(lmax+1))/2,2*(lmax+1)-1];
    end
    %% y and mat convertions
    function mat=y2mat(y)
      mat=zeros(sqrt(numel(y)));
      mat(:)=y;
    end
    function y=mat2y(mat)
      y=transpose(mat(:));
    end
    %% cs and mat convertions
    function cs=mat2cs(mat)
      S=transpose(triu(mat,1));
      cs=struct(...
        'C',tril(mat,0),...
        'S',[zeros(size(mat,1),1),S(:,1:end-1)]...
      );
    end
    function mat=cs2mat(cs)
      mat=cs.C+transpose([cs.S(:,2:end),zeros(size(cs.C,1),1)]);
    end
    %% cs and tri convertions
    function tri=cs2tri(cs)
      tri=[fliplr(cs.S(:,2:end)),cs.C];
    end
    function cs=tri2cs(tri)
      n=(size(tri,2)+1)/2;
      cs=struct(...
        'S',[zeros(n,1),fliplr(tri(:,1:n-1))],...
        'C',tri(:,n:end)...
      );
    end
    %% cs and mod convertions
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
      %make room
      out=struct(...
        'C',zeros(max(mod(:,1))+1),...
        'S',zeros(max(mod(:,1))+1)...
      );
      %propagate
      for i=1:size(mod,1)
        out.C(mod(i,1)+1,mod(i,2)+1) = mod(i,3);
        out.S(mod(i,1)+1,mod(i,2)+1) = mod(i,4);
      end          
    end
    %% agregator routines
    %data type converter
    function out=dtc(from,to,in)
      %trivial call
      if strcmpi(from,to)
        out=in;
        return
      end
      %check input
      if ~gravity.dtv(from,in)
        error([mfilename,': invalid data of type ''',from,'''.'])
      end
      %convert to required types
      switch lower(from)
        case 'y'
          switch lower(to)
            case 'mat'; out=gravity.y2mat(in);
            case 'cs';  out=gravity.mat2cs(gravity.y2mat(in));
            case 'tri'; out=gravity.cs2tri(gravity.mat2cs(gravity.y2mat(in)));
            case 'mod'; out=gravity.cs2mod(gravity.mat2cs(gravity.y2mat(in)));
          end
        case 'mat'
          switch lower(to)
            case 'y';   out=gravity.mat2y(in);
            case 'cs';  out=gravity.mat2cs(in);
            case 'tri'; out=gravity.cs2tri(gravity.mat2cs(in));
            case 'mod'; out=gravity.cs2mod(gravity.mat2cs(in));
          end
        case 'cs'
          switch lower(to)
            case 'y';   out=gravity.mat2y(gravity.cs2mat(in));
            case 'mat'; out=gravity.cs2mat(in);
            case 'tri'; out=gravity.cs2tri(in);
            case 'mod'; out=gravity.cs2mod(in);
          end
        case 'tri'
          switch lower(to)
            case 'y';   out=gravity.mat2y(gravity.cs2mat(gravity.tri2cs(in)));
            case 'mat'; out=gravity.cs2mat(gravity.tri2cs(in));
            case 'cs';  out=gravity.tri2cs(in);
            case 'mod'; out=gravity.cs2mod(gravity.tri2cs(in));
          end
        case 'mod'
          switch lower(to)
            case 'y';   out=gravity.mat2y(gravity.cs2mat(gravity.mod2cs(in)));
            case 'mat'; out=gravity.cs2mat(gravity.mod2cs(in));
            case 'cs';  out=gravity.mod2cs(in);
            case 'tri'; out=gravity.cs2tri(gravity.mod2cs(in));
          end
      end
    end
    %data type validity
    function c=dtv(type,in)
      switch lower(type)
      case 'y';  c=gravity.y_valid(in);
      case 'mat';c=gravity.mat_valid(in);
      case 'cs'; c=gravity.cs_valid(in);
      case 'tri';c=gravity.tri_valid(in);
      case 'mod';c=gravity.mod_valid(in);
      otherwise
        error([mfilename,': unknown data type ''',from,'''.'])
      end
    end
    %data type lmax
    function c=dtlmax(type,in)
      switch lower(type)
      case 'y';  c=gravity.y_lmax(in);
      case 'mat';c=gravity.mat_lmax(in);
      case 'cs'; c=gravity.cs_lmax(in);
      case 'tri';c=gravity.tri_lmax(in);
      case 'mod';c=gravity.mod_lmax(in);
      otherwise
        error([mfilename,': unknown data type ''',from,'''.'])
      end
    end
    %data type length
    function s=dtlength(type,in)
      switch lower(type)
      case 'y';  s=gravity.y_length(in);
      case 'mat';s=gravity.mat_length(in);
      case 'cs'; s=gravity.cs_length(in);
      case 'tri';s=gravity.tri_length(in);
      case 'mod';s=gravity.mod_length(in);
      otherwise
        error([mfilename,': unknown data type ''',from,'''.'])
      end
    end
    %% degree/order index mapping
    function out=mapping(lmax)
      %create triangular matrix for degrees
      d=(0:lmax)'*ones(1,2*lmax+1).*gravity.dtc('mat','tri',ones(lmax+1));
      %create triangular matrix for orders
      o=ones(lmax+1,1)*(-lmax:lmax).*gravity.dtc('mat','tri',ones(lmax+1));
      %convert to vector
      d=gravity.dtc('tri','y',d);
      o=gravity.dtc('tri','y',o);
      out=[d;o];
    end
    function [l,u]=labels(lmax,units_str)
      map=gravity.mapping(lmax);
      l=cell(1,size(map,1));
      for i=1:size(map,2)
        if map(2,i)<0
          l{i}=['S',num2str(map(1,i)),',',num2str(-map(2,i))];
        else
          l{i}=['C',num2str(map(1,i)),',',num2str( map(2,i))];
        end
      end
      if nargin>1 && nargout>1
        u=cell(size(l));
        u(:)={units_str};
      end
    end
    function out=colidx(d,o,lmax)
      if any(size(d)~=size(o))
        error([mfilename,': inputs ''d'' and ''o'' must have the same size(s).'])
      end
      m=gravity.mapping(lmax);
      out=false(1,size(m,2));
      for i=1:numel(d)
        out=out | (m(1,:)==d(i) & m(2,:)==o(i));
      end
      out=find(out);
    end
    %% constructors 
    function obj=unit(lmax,varargin)
      %create argument object, declare and parse parameters, save them to obj
      v=varargs.wrap('sources',{gravity.parameters('obj'),...
        {...
          'scale',           1,              @num.isscalar;...
          'scale_per_degree',ones(lmax+1,1), @(i) isvector(i) && lmax+1 == numel(i);...
          'scale_per_coeff', ones(lmax+1),   @(i) ismatrix(i) && all([lmax+1,lmax+1] == size(i));...
          't',               datetime('now'),@(i) isdatetime(i) || isvector(i);...
        }...
      },varargin{:});
      %create unitary triangular matrix
      u=gravity.dtc('mat','tri',ones(lmax+1));
      %scale along degrees (if needed)
      if any(v.scale_per_degree(:)~=1)
        u=v.scale_per_degree(:)*ones(1,size(u,2)).*u;
      end
      %scale per coefficient (if needed)
      if any(v.scale_per_coeff(:)~=1)
        u=gravity.dtc('mat','tri',v.scale_per_coeff).*u;
      end
      %replicate by the nr of elements of t
      u=ones(numel(v.t),1)*gravity.dtc('tri','y',u);
      %initialize
      obj=gravity(v.t,u,v.delete('t').varargin{:});
      % save the arguments v into this object
      obj=v.save(obj,{'t','lmax'});
      %call upstream scale method for global scale
      obj=obj.scale(v.scale);
    end
    % creates a unit model with per-degree amplitude equal to 1
    function obj=unit_amplitude(lmax,varargin)
      obj=gravity.unit(lmax,'scale_per_degree',1./sqrt(2*(0:lmax)+1),varargin{:});
    end
    % creates a unit model with per-degree RMS equal to 1
    function obj=unit_rms(lmax,varargin)
      obj=gravity.unit(lmax,'scale_per_degree',gravity.unit(lmax).drms,varargin{:});
    end
    % create a model with coefficients following Kaula's rule of thumb
    function obj=kaula(lmax,varargin)
      obj=gravity.unit(lmax,'scale_per_degree',[0,1e-5./(1:lmax).^2],varargin{:});
    end  
    % Creates a random model with mean 0 and std 1 (per degree)
    function obj=unit_randn(lmax,varargin)
      obj=gravity.unit(lmax,'scale_per_coeff',randn(lmax+1),varargin{:});
    end
    function obj=nan(lmax,varargin)
      obj=gravity.unit(lmax,'scale',nan,varargin{:});
    end
    function [m,e]=load(file_name,fmt,time,force,force_time)
      %default type
      if ~exist('fmt','var') || isempty(fmt) || strcmp(fmt,'auto')
        [~,fn,fmt]=fileparts(file_name);
        %get rid of the dot
        fmt=fmt(2:end);
        %check if this is CSR format
        if strcmp(fn,'GEO') || strcmp(fmt,'.GEO')
          fmt='csr';
        end
      end
      %default time
      if ~exist('time','var') || isempty(time)
        time=datetime('now');
      end
      %default force
      if ~exist('force','var') || isempty(force)
        force=false;
      end
      %default force_time
      if ~exist('force_time','var') || isempty(force_time)
        force_time=false;
      end
      %handle mat files
      [~,~,ext]=fileparts(file_name);
      if strcmp(ext,'.mat')
        mat_filename=file_name;
        file_name=strrep(file_name,'.mat','');
      else
        mat_filename=[file_name,'.mat'];
      end
      %check if mat file is already available
      if ~file.exist(mat_filename) || force
        switch lower(fmt)
        case 'gsm'
          [m,e]=load_gsm(file_name,time);
        case {'csr','geo'}
          [m,e]=load_csr(file_name,time);
        case {'icgem','gfc'}
          [m,e]=load_icgem(file_name,time);
        case 'mod'
          [m,e]=load_mod(file_name,time);
        case 'esamtm'
          error('BUG TRAP: The ''esamtm'' format is always storred in mat files')
        otherwise
          error([mfilename,': cannot handle models of type ''',fmt,'''.'])
        end
        try
          save(mat_filename,'m','e')
        catch
          disp(['Could not save ''',mat_filename,'''.'])
        end
      else
        switch lower(fmt)
        case 'esamtm'
          [m,e]=load_esamtm(mat_filename,time);
        otherwise
          %NOTICE: input argument 'time' is ignored here; only by coincidence (or design,
          %        e.g. if gravity.load_dir is used) will time be same as the one saved
          %        in the mat file.
          load(mat_filename,'m','e')
          %handle particular cases
          if ~exist('m','var')
            if exist('sol','var')
              m=sol.mod.(sol.names{1}).dat;
              e=[];
            else
              error(['Cannot handle mat file ',mat_filename,'; consider deleting so it is re-generated'])
            end
          end
        end
      end
      %enforce input 'time', if requested
      %NOTICE: this is used to be done automatically when loading the mat file
      %        (and there's no practical use for it at the moment)
      if force_time
        if m.t~=time;m.t=time;end 
        if ~isempty(e) && e.t~=time;e.t=time;end 
      end
      %update start/stop for static fields, i.e. those with time=gravity.static_start_date
      if time==gravity.static_select_date
        %enforece static start date
        m.t=gravity.static_start_date;
        %duplicate model and set it to static stop date
        m_stop=m;
        m_stop.t=gravity.static_stop_date;
        %append them together
        m=m.append(m_stop);
        %do the same to the error model if there is one
        if ~isempty(e)
          e.t=gravity.static_start_date;
          e_stop=e;
          e_stop.t=gravity.static_stop_date;
          e=e.append(e_stop);
        end
      end
    end
    function out=parse_epoch_grace(filename)
      [~,f]=fileparts(filename);
      %GSM-2_2015180-2015212_0027_JPLEM_0001_0005.gsm
      %123456789012345678901234567890
      start=dateshift(...
        time.doy2datetime(...
          str2double(f(7 :10)),...
          str2double(f(11:13))...
        ),...
      'start','day');
      stop =dateshift(...
        time.doy2datetime(...
          str2double(f(15:18)),...
          str2double(f(19:21))...
        ),...
      'end','day');
      out=mean([start,stop]);
    end
    function out=parse_epoch_aiub(filename)
      [~,f]=fileparts(filename);
      %AIUB_swarmABC_201312_70x70.gfc
      %123456789012345678901234567890
      start=dateshift(datetime([f(15:18),'-',f(19:20),'-01']),'start','month');
      stop =dateshift(datetime([f(15:18),'-',f(19:20),'-01']),'end',  'month');
      out=mean([start,stop]);
    end
    function out=parse_epoch_aiub_slr(filename)
      [~,f]=fileparts(filename);
      %AIUB-SLR_0301_1010.gfc
      %123456789012345678901234567890
      start=dateshift(datetime(['20',f(10:11),'-',f(12:13),'-01']),'start','month');
      stop =dateshift(datetime(['20',f(10:11),'-',f(12:13),'-01']),'end',  'month');
      out=mean([start,stop]);
    end
    function out=parse_epoch_asu(filename)
      [~,f]=fileparts(filename);
      %asu-swarm-2014-10-nmax40-orbits-aiub.gfc
      %123456789012345678901234567890
      start=dateshift(datetime([f(11:14),'-',f(16:17),'-01']),'start','month');
      stop =dateshift(datetime([f(11:14),'-',f(16:17),'-01']),'end',  'month');
      out=mean([start,stop]);
    end
    function out=parse_epoch_ifg(filename)
      [~,f]=fileparts(filename);
      %coeff-2015-03-SaSbSc-AIUB.gfc
      %123456789012345678901234567890
      start=dateshift(datetime([f(07:10),'-',f(12:13),'-01']),'start','month');
      stop =dateshift(datetime([f(07:10),'-',f(12:13),'-01']),'end',  'month');
      out=mean([start,stop]);
    end
    function out=parse_epoch_csr(filename)
      %2016-11.GEO.716424
      %12345678901234567890
      [~,file]=fileparts(filename);
      year =str2double(file(1:4));
      month=str2double(file(6:7));
      out=gravity.CSR_RL05_date(year,month);
    end
    function out=parse_epoch_gswarm(filename)
      [~,file]=fileparts(filename);
      persistent kv
      if isempty(kv)
      %GSWARM_GF_SABC_COMBINED_2014-11_01.gfc
      %0000000001111111111222222222233333333334
      %1234567890123456789012345678901234567890
        kv={...
          'COMBINED',25;...
          'AIUB',    21;...
          'ASU',     20;...
          'IFG',     20;...
          'OSU',     20;...
          'IGG',     20;...
        };
      end
      idx=0;
      for i=1:size(kv,1)
        if strcmp(file(16:16+length(kv{i,1})-1),kv{i,1})
          idx=kv{i,2}; break
        end
      end
      assert(idx>0,['Cannot find any of the keywords ''',strjoin(kv(:,1),''', '''),''' in the filename ''',filename,'''.'])
      year =str2double(file(idx  :idx+3));
      month=str2double(file(idx+5:idx+6));
      out=datetime(year,month,1);
      out=out+(dateshift(out,'end','month')-out)/2;
    end
    function out=parse_epoch_tn14(filename)
      %NOTICE: these data live in ~/data/SLR/TN-14/ (not set explicitly)
      %NOTICE: these data need to be converted with ~/data/SLR/TN-14/convert_TN-14.sh 
      %NOTICE: graceC20 handles these data directly from the published data file
      %TN-14_C30_C20_GSFC_SLR.55942.5.gfc
      %1234567890123456789012345678901234567890
      [~,file]=fileparts(filename);
      out=time.ToDateTime(str2double(file(24:30)),'mjd');
    end
    function out=parse_epoch_esamtm(filename)
      %mtmshc_A_19950104_06.180.mat
      [~,file]=fileparts(filename);
      file=strsplit(file,'_');
      out=time.ToDateTime([file{3},'T',file{4}(1:2),'0000'],'yyyyMMdd''T''HHmmss');
    end
    function [m,e]=load_dir(dirname,format,date_parser,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'dirname',     @ischar);
      p.addRequired( 'format',      @ischar);
      p.addRequired( 'date_parser', @(i) isa(i,'function_handle'));
      p.addParameter('wildcarded_filename',['*.',format], @ischar);
      p.addParameter('descriptor',        'unknown',     @ischar);
      %NOTICE: start/stop is only used to avoid loading models outside a certain time range
      p.addParameter('start', [], @(i) isempty(i) || (isdatetime(i)  &&  isscalar(i))); 
      p.addParameter('stop',  [], @(i) isempty(i) || (isdatetime(i)  &&  isscalar(i)));
      p.addParameter('overwrite_common_t',  false, @islogical);
      p.parse(dirname,format,date_parser,varargin{:})
      %retrieve all gsm files in the specified dir
      filelist=cells.scalar(file.unwrap(fullfile(dirname,p.Results.wildcarded_filename)),'set');
      assert(~isempty(filelist{1}),['Need valid dir, not ''',fileparts(filelist{1}),'''.'])
      %this counter is needed to report the duplicate models correctly
      c=0;init_flag=true;
      %loop over all models
      for i=1:numel(filelist)
        %skip png and yamls files (some files inheret the name of the models)
        [~,~,ext]=fileparts(filelist{i});
        if cells.isincluded({'.png','.yaml'},ext); c=c+1; continue; end
        %get time of the model in this file
        if strcmpi(func2str(p.Results.date_parser),'static')
          %if a static field is requested, there should be only one file
          if numel(filelist)~= 1
            error([mfilename,': when requested a static field, can only handle a single file, not ',...
              num2str(numel(filelist)),'.'])
          end
          %patch missing start epoch (static fields have no epoch)
          t=gravity.static_select_date;
        else
          t=p.Results.date_parser(filelist{i});
          %skip if this t is outside the required range (or invalid)
          if isempty(t) || ( time.isfinite(p.Results.start) && time.isfinite(p.Results.stop) && (t<p.Results.start || p.Results.stop<t) )
            c=c+1; continue
          end
        end
        %user feedback
        [~,f]=fileparts(filelist{i});
        disp(['Loading  ',f,' (to date ',datestr(t),')'])
        %load the data
        if init_flag
          %init output objects
          [m,e]=gravity.load(filelist{i},p.Results.format,t);
          %init no more
          init_flag=false;
        else
          %use temp container
          [m1,e1]=gravity.load(filelist{i},p.Results.format,t);
          %check if there are multiple models defined at the same epoch
          if any(m.istavail(m1.t))
            %find the model with the same epoch that has already been loaded
            [~,f_saved  ]=fileparts(filelist{m.idx(m1.t)+c});
            if p.Results.overwrite_common_t
              disp(['Replacing ',f_saved,' with ',f,' (same epoch).'])
            else
              disp(['Ignoring ',f,' because this epoch was already loaded from model ',f_saved,'.'])
              c=c+1; continue
            end              
          end
          %ensure R and GM are compatible append to output objects
          m1=m1.scale(m);
          m=m.append(m1.set_lmax(m.lmax),p.Results.overwrite_common_t);
          %same for error models, if there
          if ~isempty(e1)
            e1=e1.scale(e);
            e=e.append(e1.set_lmax(e.lmax),p.Results.overwrite_common_t);
          end
        end
      end
      %fix some parameters
      m.origin=dirname;
      e.origin=dirname;
      m.descriptor=p.Results.descriptor;
      e.descriptor=['error of ',p.Results.descriptor];
    end
    %% retrieves the Monthly estimates of C20 from 5 SLR satellites based on GRACE RL05/RL06 models
    function out=graceC20(varargin)
      %parse arguments that are required later
      v=varargs.wrap('sources',{...
        {...
          'time',        [], @(i) isdatetime(i) || isempty(i);...
          'mode',    'read', @ischar;...
          'descriptor',  '', @ischar;...
          'version','TN-14', @(i) ischar(i) || iscellstr(i);...
          'start',time.zero_date, @isdatetime;...
          'stop',  time.inf_date, @isdatetime;...
          'C20mean',-4.8416945732000E-04, @isscalar;...
        },...
      },varargin{:});
      %parse using a model or the original data
      if contains(v.version,'-model')
        new_mode=['model-',strrep(v.mode,'model-','')];
        new_version=strrep(v.version,'-model','');
        if strcmp(str.quote(v.mode),str.quote(new_mode))
          str.say('WARNING: input version has been updated from',str.quote(v.version),'to',str.quote(new_version),'.')
        else
          str.say('WARNING: over-writing input mode',str.quote(v.mode),'with',str.quote(new_mode),...
            'since input version is',str.quote(v.version),', now passed along as',str.quote(new_version),'.')
        end
        out=gravity.graceC20(varargin{:},'mode',new_mode,'version',new_version);
        return
      end
      switch v.mode
      case 'C20mean'
        out=v.C20mean; 
      case 'model-poly'
        out=2;
      case 'model-periods-datfile'
        [p,n]=fileparts(GetGRACEC20('mode','data_file','version',v.version));
        out=fullfile(p,[n,'_periods.mat']);
      case {'model-compute','model-periods'}
        %get data file
        f=gravity.graceC20(varargin{:},'mode','model-periods-datfile');
        %check if periods were already computed
        if ~file.exist(f)
          %loading necessary data
          c20=gravity.graceC20(varargin{:},'mode','read');
           np=gravity.graceC20(varargin{:},'mode','model-poly');
          %compute periods
          [~,pd]=c20.parametric_decomposition_search('np',np,'T',[365.2426,182.6213],'timescale','days');
          out=pd.T;
          %save periods
          save(f,'out');
        else
          disp(['Loading parametric decomposition from ',f])
          %load periods
          load(f,'out');
        end
%         %build strings to save to file
%         h=strsplit(c20.descriptor,newline);
%         s=strjoin([...
%             cellfun(@(i) ['%',i],h(:),'UniformOutput',false);...
%             {'out=[...'};...
%             arrayfun(@(i) ['  ',num2str(i,'%.12e'),';...'],out(:),'UniformOutput',false);...
%             {'];'};...
%           ],newline);
%         f=gravity.graceC20(varargin{:},'mode','model-periods-datfile');
%         b=file.strsave(f,s);
%         str.say('Written',b,'bytes of data to file:',newline,f,newline,'related to the periods of:',newline,c20.descriptor);
      case 'model-datfile'
        [p,n]=fileparts(GetGRACEC20(varargin{:},'mode','data_file'));
        out=fullfile(p,[n,'_pd.mat']);
      case 'model-md5file'
        [p,n]=fileparts(GetGRACEC20(varargin{:},'mode','data_file'));
        out=fullfile(p,[n,'.md5']);
      case 'model-md5'
        out=file.md5(GetGRACEC20(varargin{:},'mode','data_file'));
      case 'model-md5set'
        out=file.strsave(...
          gravity.graceC20(varargin{:},'mode','model-md5file'),...
          gravity.graceC20(varargin{:},'mode','model-md5')...
        );
      case 'model-md5get'
        md5file=gravity.graceC20(varargin{:},'mode','model-md5file');
        if ~file.exist(md5file)
          gravity.graceC20(varargin{:},'mode','model-md5set');
        end
        out=file.strload(md5file);
      case 'model-md5check'
        out=strcmp(...
          gravity.graceC20(varargin{:},'mode','model-md5get'),...
          gravity.graceC20(varargin{:},'mode','model-md5')...
        );
      case {'model','model-get','model-set','model-read','model-reload'}
        %loading necessary data
        c20=gravity.graceC20(varargin{:},'mode','read');
         np=gravity.graceC20(varargin{:},'mode','model-poly');
          T=gravity.graceC20(varargin{:},'mode','model-periods');
        %check if pdset is already available
        f_pdset=gravity.graceC20(varargin{:},'mode','model-list-datfile');
        f_pd   =gravity.graceC20(varargin{:},'mode','model-datfile');
        if ~file.exist(f_pdset) || ~file.exist(f_pd) || ~...
          gravity.graceC20(varargin{:},'mode','model-md5check')  
          %get the coefficients; NOTICE: always use c20.t so that f_pdset is not dependent on inputs
          [~,pd_set]=c20.parametric_decomposition('np',np,'T',T,...
            'timescale','days','time',c20.t_domain(days(7))); 
          %save them
          save(f_pdset,'pd_set')
          %update md5 of data
          gravity.graceC20(varargin{:},'mode','model-md5set')  
        else
          load(f_pdset,'pd_set')
        end
        %handle different time domains
        if isempty(v.time);      v.time=c20.t; end
        %implement default start/stop for modelling
        if time.iszero(v.start); v.start=v.time(1  ); end
        if time.isinf( v.stop );  v.stop=v.time(end); end
        %enforce start/stop
        v.time=v.time(v.time>=v.start & v.time<=v.stop);
        %evaluate model at requested time domain
        out=pardecomp.join( pd_set,'time',v.time);
      case 'model-list-datfile'
        [p,n]=fileparts(GetGRACEC20(varargin{:},'mode','data_file'));
        out=fullfile(p,[n,'_pdset.mat']);
      case {'model-list','model-list-tex'}
        %check if pdset is already available
        f=gravity.graceC20(varargin{:},'mode','model-list-datfile');
        if ~file.exist(f)
          %compute the model (saving is done inside)
          gravity.graceC20(varargin{:},'mode','model');
        end
        %load it
        load(f,'pd_set')
        %output it
        switch v.mode
          case 'model-list';     out=pardecomp.table(pd_set,'tablify',true);
          case 'model-list-tex'; out=pardecomp.table(pd_set,'tablify',false,'latex_table',true);
        end       
      case 'model-plot'
        %retrieve the orignal data
        c20o=gravity.graceC20(varargin{:},'mode','read');
        %resample to a finer time domain
        c20r=c20o.interp(c20o.t_domain(days(7)),...
          'interp_over_gaps_narrower_than',days(45),...
          'interp1_args',{'pchip'}...
        );
        c20m=gravity.graceC20(varargin{:},'mode','model','time',c20o.t_domain(days(7)));
        c20e=c20r-c20m;
        %compute means
        c20om=c20o.stats('mode','mean');
        c20mm=c20m.stats('mode','mean');
        c20em=c20e.stats('mode','mean');
        %compute stds
        c20os=c20o.stats('mode','std');
        c20ms=c20m.stats('mode','std');
        c20es=c20e.stats('mode','std');
        %remove the mean
        c20o=c20o-c20om(1);
        c20m=c20m-c20mm(1);
        c20e=c20e-c20em(1);
        %plot it
        plotting.figure;
        c20o.addgaps(days(35)).plot('columns',1);
        c20m.plot('columns',1);
        c20e.plot('columns',1);
        plotting.enforce(varargin{:},...
          'plot_legend',...
          {...
            ['Original \mu=',num2str(c20om(1),'%.3e'),' \sigma=',num2str(c20os(1),'%.3e')],...
            ['Modelled \mu=',num2str(c20mm(1),'%.3e'),' \sigma=',num2str(c20ms(1),'%.3e')],...
            ['Residual \mu=',num2str(c20em(1),'%.3e'),' \sigma=',num2str(c20es(1),'%.3e')],...
          },...
          'plot_xlabel','none',...
          'plot_title',c20o.descriptor...
        );
        out=c20e;
      case 'interp'
        assert(~isempty(v.time),['If mode is ''',v.mode,''', need argument ''time''.'])
        out=gravity.graceC20(varargin{:},'mode','read');
        out=out.interp(v.time);
      case 'plot-all'
        if iscellstr(v.version)
          version_list=v.version;
        else
          version_list={'GSFC-7day','GSFC','CSR-RL06','TN-14','TN-11','TN-07'};
        end
        dat_list=cell(size(version_list));
        plotting.figure(varargin{:});
        for i=1:numel(version_list)
          %pick keywords
          kw=strsplit(version_list{i},'-');
          %branch on them (NOTICE: always last)
          switch kw{end}
          case 'model'
            if i>1
              %NOTICE: to get the model with the most complete time domain, out it last in the 'version' input
              dat_list{i}=gravity.graceC20('mode','model',...
                'version',strrep(version_list{i},'-model',''),...
                'time',time.union(cellfun(@(j) j.t,dat_list(1:i-1),'UniformOutput',false),days(7))...
              );
            else
              dat_list{i}=gravity.graceC20('mode','model',...
                'version',strrep(version_list{i},'-model','')...
              );
            end              
            otherwise
            dat_list{i}=gravity.graceC20('mode','read','version',version_list{i}, v.delete({'version','mode'}).varargin{:});
            dat_list{i}=dat_list{i}.interp(dat_list{i}.t_domain(days(7)),...
              'interp_over_gaps_narrower_than',days(45)...
            );
          end
          dat_list{i}=dat_list{i}-v.C20mean;
          dat_list{i}.plot('columns',1);
        end
        plotting.enforce(varargin{:},...
          'plot_legend',version_list,...
          'plot_xlabel','none');
        title(['mean C_{20}: ',num2str(v.C20mean,'%e')]);
        out=dat_list;
      otherwise
        %call mother routine
        try
          [t,s,e,d]=GetGRACEC20(varargin{:},'mode',v.mode);
        catch ME
          switch ME.identifier
            case 'MATLAB:unassignedOutputs'
              t=GetGRACEC20(varargin{:},'mode',v.mode);
              s=[];
            otherwise
              rethrow(ME)
          end
        end
        if ~isempty(s)
          %create time series
          out=simpletimeseries(t,[s,e],...
            'labels',{'C20','error C20'},...
            'units',{'',''},...
            'timesystem','gps',...
            'descriptor',d...
          );
        else
          out=t;
        end
      end
    end
    %% vector operations to make models compatible
    function lmax=vlmax(model_list)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'model_list', @(i) iscell(i) && all(cellfun(isa(i,'gravity'))))
      p.parse(model_list)
      lmax=min(cell2mat(cellfun(@(i) i.lmax,model_list)));
    end
    function gm=vGM(model_list)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'model_list', @(i) iscell(i) && all(cellfun(isa(i,'gravity'))))
      p.parse(model_list)
      gm=min(cell2mat(cellfun(@(i) i.GM,model_list)));
    end
    function r=vR(model_list)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'model_list', @(i) iscell(i) && all(cellfun(isa(i,'gravity'))))
      p.parse(model_list)
      r=min(cell2mat(cellfun(@(i) i.R,model_list)));
    end
    %% common operations
    function mod=common_ops(mode,mod,varargin)
      %handle vector call on mode
      if iscellstr(mode)
        for i=1:numel(mode)
          mod=gravity.common_ops(mode{i},mod,varargin{:});
        end
      else
        % parse input arguments
        v=varargs.wrap('sources',{...
          {...
            'max_degree',            -1 , @(i) isnumeric(i) && isscalar(i);...
            'use_GRACE_C20',     'none' , @ischar; ...
            'use_GRACE_C20_plot', false , @islogical; ...
            'date_parser',       'none' , @ischar; ...
            'delete_C00',         false , @islogical; ...
            'delete_C20',         false , @islogical; ...
            'start',     time.zero_date , @isdatetime; ...
            'stop',       time.inf_date , @isdatetime; ...
            'static_model',      'none' , @ischar; ...
            'product_name',    'unknown', @ischar;...
            'model_type',       'signal', @ischar;...
            'show_msg',            true , @islogical;...
          },...
        },varargin{:});
        %sanity
        if exist('mod','var') && ~isa(mod,'gravity')
          str.say(['Ignoring model ',mod.descriptor,' since it is of class ''',class(mod),''' and not of class ''gravity''.'])
          return
        end
        %inits user feedback vars
        msg='';
        switch mode
        case 'mode_list'
          mod={'consistent_GM','consistent_R','max_degree','use_GRACE_C20','delete_C00','delete_C20','start','stop','static_model'};
          v.show_msg=false;
        case 'all'
          mod=gravity.common_ops(gravity.common_ops('mode_list'),mod,v.varargin{:},'show_msg',false);
          return
        case 'consistent_GM'
          mod=mod.setGM(gravity.parameters('GM'));
        case 'consistent_R'
          mod=mod.setR( gravity.parameters('R'));
        case 'max_degree'
          %set maximum degree (if requested)
          if v.max_degree>0
            mod.lmax=v.max_degree;
            msg=['set to ',num2str(v.max_degree)];
          end
        case 'use_GRACE_C20'
          %set C20 coefficient
          if ~str.none(v.use_GRACE_C20) && strcmp(v.model_type,'signal')
            %some sanity
            if strcmpi(v.date_parser,'static')
              error([mfilename,': there''s no point in replacing GRACE C20 coefficients in a static model.'])
            end
%             %legacy support
%             if islogical(v.use_GRACE_C20) && v.use_GRACE_C20
%               v.use_GRACE_C20='TN-07';
%             end
            %parse using a model or the original data
            if contains(v.use_GRACE_C20,'-model')
              c20mode='model';
              v.use_GRACE_C20=strrep(v.use_GRACE_C20,'-model','');
            else
              c20mode='interp';
            end
            %get C20 timeseries, interpolated to current time domain
            c20=gravity.graceC20('version',v.use_GRACE_C20,'mode',c20mode,'time',mod.t);
            if v.use_GRACE_C20_plot
              figure
              plot(c20.x_masked,c20.y_masked([],1),'x-','MarkerSize',10,'LineWidth',4), hold on
              plot(c20.x,spline(c20.x_masked,c20.y_masked([],1),c20.x),'o-','MarkerSize',10,'LineWidth',2)
              plot(c20.x,pchip(c20.x_masked,c20.y_masked([],1),c20.x),'+-','MarkerSize',10,'LineWidth',2)
              legend('data','spline','pchip')
            end
            mod=mod.setC(2,0,c20.y(:,1),mod.t);
            if v.use_GRACE_C20_plot
              figure
              plot(c20.t,c20.y(:,1),'o-'), hold on
              plot(mod.t,mod.C(2,0),'x-')
              legend('TN-11',mod.descriptor)
            end
            msg=['set to ',v.use_GRACE_C20];
          end
        case 'delete_C00'
          %remove C00 bias
          if v.delete_C00
            mod=mod.setC(0,0,0);
          end
        case 'delete_C20'
          %remove C20 bias
          if v.delete_C20
            mod=mod.setC(2,0,0);
          end
        case 'start' %NOTICE: this is only done when loading the raw data (afterwards the matlab data is read directly, bypassing this routine altogher)
          if v.start~=time.zero_date
            if v.start<mod.start
              %append extremeties
              mod=mod.append(gravity.nan(mod.lmax,'t',v.start,'R',mod.R,'GM',mod.GM));
            elseif v.start>mod.start
              %trim extremeties (this is redundant unless data is saved before the start metadata is increased)
              mod=mod.trim(v.start,mod.stop);
            end
            msg=['at ',datestr(v.start)];
          end
        case 'stop'  %NOTICE: this is only done when loading the raw data (afterwards the matlab data is read directly, bypassing this routine altogher)
          if v.stop~=time.inf_date
            if v.stop>mod.stop
              %append extremeties
              mod=mod.append(gravity.nan(mod.lmax,'t',v.stop,'R',mod.R,'GM',mod.GM));
            elseif v.stop<mod.stop
              %trim extremeties (this is redundant unless data is saved before the stop metadata is decreased)
              mod=mod.trim(mod.start,v.stop);
            end
            msg=['at ',datestr(v.stop)];
          end
        case 'static_model'
          %remove static field (if requested)
          if ~strcmpi(v.static_model,'none') && strcmp(v.model_type,'signal')
            %TODO: make this handle filenames as well, not only datastorage products
            static=datastorage().init(v.static_model,'show_msg',false).data_get_scalar(datanames(v.static_model).append_field_leaf('signal'));
            %adjust the start/stop
            switch static.length
            case 1
              static2=static;
              static.t=mod.start;
              static2.t=mod.stop;
              static=static.append(static2);
            case 2
              static.t=[mod.start;mod.stop];
            otherwise
              error(['Cannot handle static object of class gravity with length ', num2str(static.length),'.'])
            end
            %make sure max degree matches
            static.lmax=mod.lmax;
            %subtract it
            mod=mod-static.interp(mod.t);
            msg=[': subtracted ',v.static_model];
          end
        otherwise
          error(['Cannot handle operantion ''',mode,'''.'])
        end
        if v.show_msg;str.say(v.product_name,':',mode,msg);end
      end
    end
    %% model combination
    function out=combine(model_list,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'model_list', @(i) iscell(i) && all(cellfun(@(i) isa(i,'gravity'),model_list)))
      p.addParameter('mode','mean',@ischar);
      p.addParameter('type','signal',@ischar);
      p.addParameter('weights',ones(size(model_list))/numel(model_list),@(i) isnumeric(i) && all(size(model_list)==size(i)))
      p.parse(model_list,varargin{:})
      %trivial call
      if numel(model_list)==1
        out=model_list{1};
        return
      end
%       %sanity
%       for i=2:numel(model_list)
%         if ~model_list{1}.istequal(model_list{i})
%           error([mfilename,': time domain discrepancy between model ',...
%             model_list{1}.descriptor,' and model ',...
%             model_list{i}.descriptor,'.'])
%         end
%       end
      %set model compatibility
      model_list=simpledata.merge_multiple(model_list);
      %branch on mode
      switch lower(p.Results.mode)
      case {'mean','am','arithmeticmean'}
        switch lower(p.Results.type)
        case 'signal'
          out=simpledata.vmean(...
            model_list,...
          'weights',p.Results.weights);
        case 'error'
          out=simpledata.vmean(...
            simpledata.vtimes(model_list,model_list),...
          'weights',p.Results.weights).sqrt;
        otherwise
          error([mfilename,': unknown type ',p.Results.type,'.'])
        end
      case {'dtm','difftomean'}
        error([mfilename,': implementation needed'])
      otherwise
        error([mfilename,': unknown mode ''',p.Results.mode,'''.'])
      end
      %update the descriptor
      msg=cell(size(model_list));
      for i=1:numel(model_list)
        msg{i}=model_list{i}.descriptor;
      end
      out.descriptor=[p.Results.mode,' of ',strjoin(msg,', ')];
    end
    %% (half) wavelength to degree convertions
    function deg=wl2deg(wl)
      deg=2*pi*gravity.parameters('R')./wl;
    end
    function wl=deg2wl(deg)
      wl=gravity.wl2deg(deg);
    end
    function deg=hwl2deg(hwl)
      deg=gravity.wl2deg(hwl*2);
    end
    function hwl=deg2hwl(deg)
      hwl=gravity.deg2wl(deg)/2;
    end
    function rad=gauss_smoothing_radius(deg)
      rad=gravity.deg2hwl(deg)/2;
    end
    function deg=gauss_smoothing_degree(rad)
      deg=round(gravity.hwl2deg(rad*2));
    end
    function out=gauss_smoothing_type(in)
      % This is a weak criteria but will work as long as there's no need
      % to smooth at degrees larger than 10000 or radii smaller than 10km.
      if round(sum(in>10e3))
        out='radius';
      else
        out='deg';
      end      
    end
    function out=gauss_smoothing_degree_translate(in)
      switch gravity.gauss_smoothing_type(in)
      case 'radius'; out=gravity.gauss_smoothing_degree(in);
      case 'deg'   ; out=in;
      end
    end
    function out=gauss_smoothing_name(in)
      switch gravity.gauss_smoothing_type(in)
      case 'radius'
        %NOTICE: this test should be to see if in is zero but that is assumed to refer to degrees
        if isfinite(in)
          out=[num2str(in/1e3),'km'];
        else
          out='no';
        end
      case 'deg'
        %NOTICE: this test should only be to see if in is inf but zero is assumed to refer to degrees
        %        even if in practice is always refers to radius
        if isfinite(in) && in~=0
          out=[num2str(round(gravity.gauss_smoothing_radius(in)/1e4)*1e1),'km'];
        else
          out='no';
        end
      end
    end
    function out=smoothing_name(in)
      switch lower(in)
      case 'gauss';  out='Gaussian';
      case 'spline'; out='Spline';
      case 'trunc';  out='Truncation';
      otherwise; error(['Cannot handle smoothing method ''',in,'''.'])
      end
    end
    %% CSR specific stuff
    function t=CSR_RL05_date(year,month,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addParameter('EstimDirs_file',fullfile(getenv('HOME'),'data','csr','mascons','RL05','EstimDirs_RL05'),@ischar);
      p.parse(varargin{:})
      %sanity on year
      if year<1000
        year=year+2000;
      end
      persistent CSR_RL05_date_table 
      if isempty(CSR_RL05_date_table)
        fid=file.open(p.Results.EstimDirs_file);
        d=textscan(fid,'%d %f %s %s %s %s %s %s %f %f');
        fclose(fid);
        years=d{2};
        months=cellfun(@(i) time.month2int(i),d{3});
        start=time.ToDateTime(d{ 9},'modifiedjuliandate');
        stop= time.ToDateTime(d{10},'modifiedjuliandate');
        mid=time.ToDateTime((d{ 9}+d{10})/2,'modifiedjuliandate');
        CSR_RL05_date_table=table(years,months,start,stop,mid);
      end
      rows=CSR_RL05_date_table.years==year & CSR_RL05_date_table.months==month;
      t=CSR_RL05_date_table{rows,{'mid'}};
    end
    function [m,e]=CSR_RL05_SH(varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addParameter('datadir',fullfile(getenv('HOME'),'data','csr','RL05'),@(i) ischar(i) && exist(i,'dir')~=0)
      p.parse(varargin{:})
      [m,e]=gravity.load_dir(p.Results.datadir,'csr',@gravity.parse_epoch_csr,'wildcarded_filename','*.GEO.*',varargin{:});
    end
    function [m,e]=CSR_Mascons(varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addParameter('lmax',60,@(i) isscalar(i) && isnumeric(i))
      p.addParameter('RL','05',@(i) @ischar)
      p.parse(varargin{:})
      p.addParameter('datadir',fullfile(getenv('HOME'),'data','csr','mascons',['RL',p.Results.RL]),@(i) ischar(i) && exist(i,'dir')~=0)
      p.parse(varargin{:})
      datafile=fullfile(p.Results.datadir,['alldata-deg',num2str(p.Results.lmax),'.mat']);
      if file.exist(datafile)
        disp(['Loading all RL',p.Results.RL,' Mascons from ',datafile,'...'])
        load(datafile,'m','e')
      else
        [m,e]=gravity.load_dir(p.Results.datadir,'csr',@gravity.parse_epoch_csr,'wildcarded_filename','*.GEO',varargin{:});
        switch p.Results.RL
          case '05'
            s=gravity.static('GIF48');
          case '06'
            s=gravity.static('GGM05C');
          otherwise
            error(['Cannot handle input ''RL'' with value ''',p.Results.RL,'''.'])
        end
        m.lmax=p.Results.lmax;
        e.lmax=p.Results.lmax;
        s.lmax=p.Results.lmax;
        m=m-s;
        disp(['Saving all RL05 Mascons to ',datafile,'...'])
        save(datafile,'m','e')
      end
    end
    function [m,e]=static(model)
      switch upper(model)
        case 'GIF48'
          datafile=fullfile(file.orbdir('auxiliary'),'GIF48.2007.GEO');
          fmt='GEO';
        case 'GGM05C'
          datafile=fullfile(file.orbdir('auxiliary'),'GGM05C.gfc');
          fmt='gcf';
        case 'GGM05G'
          datafile=fullfile(file.orbdir('auxiliary'),'ggm05g.gfc');
          fmt='gcf';
        case 'WGS84'
          datafile=fullfile(file.orbdir('auxiliary'),'WGS84.GEO');
          fmt='GEO';
        otherwise
          error(['Cannot handle static model ''',model,'''.'])
      end
      [m,e]=gravity.load(datafile,fmt);
    end
    %% ESA's Earth System Model
    function out=ESA_MTM_dir(year,month,component)
      out=fullfile(getenv('HOME'),'media','data','data','ESA-MTM','data',num2str(year,'%04i'),component,num2str(month,'%02i'));
    end
    %year      : year (numeric)
    %month     : month (numeric)
    %component : A, AOHIS, H, I, O or S (char)
    function [m,e]=ESA_MTM(year,month,component,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addParameter('datadir',gravity.ESA_MTM_dir(year,month,component),@(i) ischar(i) && exist(i,'dir')~=0)
      p.parse(varargin{:})
      [m,e]=gravity.load_dir(p.Results.datadir,'esamtm',@gravity.parse_epoch_esamtm,'wildcarded_filename','mtmshc_*.180.mat',varargin{:});
    end
    function m=ESA_MTM_all(component,varargin)
      %init year-loop variables
      firstyear=true;
      %loop over the years
      for year=1995:2006
        %define yearly-aggregated data
        fileyear=fullfile(gravity.ESA_MTM_dir(year,[],component),['gravity_',component,'_',num2str(year,'%04i'),'.mat']);
        %load aggregated data is available
        if file.exist(fileyear)
          load(fileyear,'m_year');
        else
          %init month variables
          firstmonth=true;
          %loop over the months
          for month=1:12
            %define monthly-aggregated data
            filemonth=fullfile(gravity.ESA_MTM_dir(year,month,component),['gravity_',component,'_',num2str(year,'%04i'),'-',num2str(month,'%02i'),'.mat']);
            %load aggregated data is available
            if file.exist(filemonth)
              load(filemonth,'m_month');
            else
              %load every day
              m_month=gravity.ESA_MTM(year,month,component,varargin{:});
              %save this month for next time
              save(filemonth,'m_month');
            end
            %create/append
            if firstmonth
              firstmonth=false;
              m_year=m_month;
            else
              m_year=m_year.append(m_month);
            end
          end 
          %save this year for next time
          save(fileyear,'m_year');
        end
        %create/append
        if firstyear
          firstyear=false;
          m=m_year;
        else
          m=m.append(m_year);
        end
      end
    end
    %% permanent (solid earth) tide
    function C20=zero_tide(C20,tide_system)
    % https://geodesyworld.github.io/SOFTS/solid.htm
    %?http://www.springerlink.com/index/V1646106J6746210.pdf
    % https://www.ngs.noaa.gov/PUBS_LIB/EGM96_GEOID_PAPER/egm96_geoid_paper.html
    % https://www.iers.org/SharedDocs/Publikationen/EN/IERS/Publications/tn/TechnNote32/tn32_057.pdf, section 6.3
    %   The degree 2 zonal tide generating potential has a mean (time average)
    %   value that is nonzero. This time independent (nm) = (20) potential
    %   produces a permanent deformation and a consequent time independent
    %   contribution to the geopotential coefficient C20.
      debug=true;
      str.say(['Tide system is ',tide_system])
      switch tide_system
      case {'zero_tide'}
        % the zero-frequency value includes the indirect distortion, but not the direct distortion
        % NOTICE: do nothing, this is the default
        delta=0;
      case {'free_tide','tide_free','tide_tide','non_tidal'}
        % the tide-free value is the quantity from which all tidal effects have been removed
        delta=-gravity.parameters('zf_love')*gravity.parameters('pt_factor');
      case 'mean_tide'
        % the mean tide value includes both the direct and indirect permanent tidal distortions
        %http://mitgcm.org/~mlosch/geoidcookbook/node9.html
        delta=gravity.parameter('pt_factor');
      otherwise
        error(['Cannot handle the permanent tide system ''',tide_system,'''.'])
      end
      str.say('disp',debug,['Converted to zero tide by adding ',num2str(delta)])
      C20=C20+delta;
    end
    %% utilities
    function [degrees,orders]=resolve_degrees_orders(varargin)
      %[degrees,orders]=gravity.resolve_degrees_orders('degrees',[-1],'orders',[inf])
      % degrees =     []
      % orders  =     []
      %[degrees,orders]=gravity.resolve_degrees_orders('degrees',[2,3],'orders',[1,inf])
      % degrees = [     2     3     3     3     3     3     3     3 ]
      % orders  = [     1    -3    -2    -1     0     1     2     3 ]
      %[degrees,orders]=gravity.resolve_degrees_orders('degrees',[2,3,2],'orders',[1,2,1])
      % degrees = [     2     3     2 ]
      % orders  = [     1     2     1 ]
      v=varargs.wrap('sources',{{...
        'lmax',       inf,@isnumeric;...
        }},varargin{:});
      if isfinite(v.lmax)
        %set defaults according to lmax given
        v=varargs.wrap('sources',{{...
          'degrees',    0:v.lmax,       @isnumeric;...
          'orders',     inf(1,v.lmax+1),@(i) isnumeric(cell2mat(cells.num(i)));...
          }},varargin{:});
      else
        %set defaults according to degrees and orders given
        v=varargs.wrap('sources',{{...
          'degrees',    [2,3],     @isnumeric;...
          'orders',     [inf,inf], @(i) isnumeric(cell2mat(cells.num(i)));...
          }},varargin{:});
      end
      %convert degree-wise orders
      degrees_out=cell(size(v.degrees));
       orders_out=cell(size(v.degrees));
      %need to convert keywords in orders ('inf'), if there
      v.orders=cell2mat(cells.num(v.orders));
      %sanitize
      assert(numel(v.degrees)==numel(v.orders),...
        ['Numel of ''degrees'' (',num2str(numel(v.degrees)),...
        ') different than numel of ''orders'' (',num2str(numel(v.orders)),').'])
      for i=1:numel(v.degrees)
        if (v.degrees(i)<0)
          %do nothing
        elseif ~isfinite(v.orders(i))
           orders_out{i}=-v.degrees(i):v.degrees(i);
          degrees_out{i}=ones(size(orders_out{i}))*v.degrees(i);
        else
           orders_out{i}=v.orders(i);
          degrees_out{i}=v.degrees(i);
        end
      end
       orders=cell2mat(cells.rm_empty(orders_out));
      degrees=cell2mat(cells.rm_empty(degrees_out));
    end
    function lmax=width2lmax(width)
      lmax=sqrt(width)-1;
    end
    %% general test for the current object
    function out=test_parameters(field,l,w)
      switch field
      case 'something'
        %To be implemented
      otherwise
        out=simpledata.test_parameters(field,l,w);
      end
    end
    function test(method,l,t)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      if ~exist('l','var') || isempty(l)
        l=4;
      end
      if ~exist('t','var') || isempty(t)
        t=[datetime('now'),datetime('now')+seconds(1)];
      end
      nr_coeff=(l+1)^2;
      disp(['--- ',method,': l=',num2str(l)])
      switch lower(method)
      case 'all'
         for i={'reps','unit','unit rms','r','gm','minus','grid','ggm05g','stats','c','smoothing','deepoceanmaskplot','gracec20'}
           gravity.test(i{1},l);
         end
      case 'reps'
        dt={'y','mat','cs','tri','mod'};
        y=randn(1,nr_coeff);
        dd={...
          y,...
          gravity.y2mat(y),...
          gravity.mat2cs(gravity.y2mat(y)),...
          gravity.cs2tri(gravity.mat2cs(gravity.y2mat(y))),...
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
      case 'unit'
        a=gravity.unit_amplitude(l,'t',t);
        disp('- C')
        disp(a.cs(numel(t)).C)
        disp('- S')
        disp(a.cs(numel(t)).S)
        disp('- tri')
        disp(a.tri{numel(t)})
        disp('- mod')
        disp(a.mod{numel(t)})
        disp('- das')
        disp(a.das)
      case 'unit rms'
        a=gravity.unit_rms(l,'t',t);
        disp('- tri')
        disp(a.tri{numel(t)})
        disp('- drms')
        disp(a.at(t(numel(t))).drms)
      case 'r'
        a=gravity.unit_amplitude(l,'t',t);
        disp('- tri: start')
        disp(a.tri{numel(t)})
        disp('- tri: 2*R')
        disp(a.scale(a.R*2,'R').tri{numel(t)})
      case 'gm'
        a=gravity.unit_amplitude(l,'t',t);
        disp('- tri: start')
        disp(a.tri{numel(t)})
        disp('- tri: 2*GM')
        disp(a.scale(a.GM*2,'GM').tri{numel(t)})
      case 'minus'
        a=gravity.unit_amplitude(l,'t',t);
        disp('- tri: a=')
        disp(a.tri{numel(t)})
        disp('- tri: a-a.scale(2)=')
        b=a-a.scale(2);
        disp(b.tri{numel(t)})
      case 'grid'
        figure
        gravity.unit_randn(120,'t',t).grid.imagesc
      case 'ggm05g'
        m=gravity.ggm05g;
        disp('- print gravity')
        m.print
        disp('- print grid')
        m.grid.print
      case 'stats'
        m=gravity.ggm05g.setC(0,0,0).setC(2,0,0);
        for i={'dmean','cumdmean','drms','cumdrms','dstd','cumdstd','das','cumdas'}
          figure
          m.plot('method',i{1},'functional','geoid','title',[m.descriptor,' - ',i{1}]);
        end
      case 'c'
        if numel(t)<3
          now=juliandate(datetime('now'),'modifiedjuliandate');
          t=datetime(now,       'convertfrom','modifiedjuliandate'):...
            datetime(now+l*10-1,'convertfrom','modifiedjuliandate');
        end
        a=gravity.unit_randn(l,'t',t);
        disp('- tri: a=')
        disp(a.tri{numel(t)})
        d=round(rand*l);
        o=round(rand*2*d)-d;
        disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(t(numel(t))),')=',num2str(a.C(d,o,t(numel(t))))])
        disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(t(1)),')=',num2str(a.C(d,o,t(1)))])
        v=9.9999;
        a=a.setC(d,o,v,t(numel(t)));
        disp(['- tri: a.setC(',num2str(d),',',num2str(o),')=',num2str(v)])
        disp(a.tri{numel(t)})
        disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(t(numel(t))),')=',num2str(a.C(d,o,t(numel(t))))])
        disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(t(1)),')=',num2str(a.C(d,o,t(1)))])
      case 'smoothing'
        a=gravity.unit(l);
        if isdatetime(t)
          t=round(l/2);
        end
        methods={'gauss','spline','trunc'};
        for i=1:numel(methods)
          a.scale_plot(t,methods{i});
        end
        legend(methods);
      case 'deepoceanmaskplot'
        a=simplegrid.unit(l*90,l*45);
        a=a.spatial_mask('deep ocean');
        c=gray(l*16);
        c=flipud(c(l*4:l*12,:));
        plotting.figure;
        a.imagesc;
        colormap(c)
        colorbar off
        plotting.enforce('plot_legend_location','none');
      case 'gracec20'
        gravity.graceC20('mode','plot-all','version',{'GSFC-7DAY','GSFC','TN-14','TN-11'});
      case 'sh2grid'
        datadir='~/data/gswarm/aiub/gravity/';
        datafiles='GSWARM_GF_SABC_COMBINED_2020-*_09.gfc';
        functional='eqwh';
        %load the data
        g=gravity.load_dir(datadir,'gfc',@gravity.parse_epoch_gswarm,'wildcarded_filename',datafiles,'descriptor','Swarm');
        %remove static model
        g=g-gravity.static('ggm05c').set_lmax(g.lmax).scale(g);
        %convert to equivalent water height
        g=g.scale(750e3,'gauss').scale(functional,'functional');
        %convert to grid
        out=g.grid('spatial_step',2).center_resample;
        %show the grid
        out.imagesc
        %export to xyz format
        out.xyz(['./GSWARM_GF_SABC_COMBINED-',functional,'.xyz'])
      otherwise
        error(['Cannot handle test method ''',method,'''.'])
      end
    end
  end
  methods
    %% constructor
    function obj=gravity(t,y,varargin)
      % input parsing
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y', @(i) simpledata.valid_y(i));
      %create argument object, declare and parse parameters, save them to obj
      [v,p]=varargs.wrap('parser',p,'sources',{gravity.parameters('obj')},'mandatory',{t,y},varargin{:});
      % get some parameters
      lmax=gravity.y_lmax(y(1,:));
      [labels,units]=gravity.labels(lmax,gravity.functional_units(p.Results.functional));
      % call superclass
      obj=obj@simpletimeseries(t,y,varargin{:},...
        'labels',labels,...
        'units',units...
      );
      % save the arguments v into this object
      obj=v.save(obj,{'t','y'});
    end
    function obj=assign(obj,y,varargin)
      %pass it upstream
      obj=assign@simpletimeseries(obj,y,varargin{:});
      %update labels and units
      obj=obj.setlabels;
    end
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in,[gravity.parameters('list');more_parameters(:)]);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      out=metadata@simpletimeseries(obj,[gravity.parameters('list');more_parameters(:)]);
    end
    function out=simpletimeseries(obj)
      % call superclass
      out=simpletimeseries(obj.t,obj.y,...
        'labels',obj.labels,...
        'units',obj.functional_unit(true),...
        'timesystem',obj.timesystem...
      );
    end
    %% labels and units
    function [obj,reset_flag]=setlabels(obj)
      [l,u]=gravity.labels(obj.lmax,obj.functional_unit);
      reset_flag=false;
      if numel(obj.labels)~=obj.width || ~cells.isequalstr(obj.labels,l)
        obj.labels=l;
        obj.units=u;
        reset_flag=true;
      end
    end
    %% functional
    function out=get.funct(obj)
      if isempty(obj.functional)
        obj.functional=varargs(gravity.parameter_list).functional;
      end
      out=obj.functional;
    end
    function obj=set.funct(obj,in)
      assert(gravity.isfunctional(in),['Ilegal functional ''',in,'''.'])
      obj.functional=in;
    end
    function out=functional_name(obj)
      out=gravity.functional_names(obj.funct);
    end
    function out=functional_unit(obj,cellstr_flag)
      if ~exist('cellstr_flag','var') || isempty(cellstr_flag)
        cellstr_flag=false;
      end
      out=gravity.functional_units(obj.funct);
      if cellstr_flag
        tmp=out;
        out=cell(1,obj.width);
        out(:)={tmp};
      end
    end
    %% lmax
    function obj=set.lmax(obj,l)
      %trivial call
      if obj.lmax==l
        return
      end
      %make room for new matrix representation
      mat_new=cell(obj.length,1);
      %get existing mat representation
      mat_old=obj.mat;
      %branch between extending and truncating
      if obj.lmax<l
        for i=1:obj.length
          %make room for this epoch
          mat_new{i}=zeros(l+1);
          %propagate existing coefficients
          mat_new{i}(1:obj.lmax+1,1:obj.lmax+1)=mat_old{i};
        end
      else
        for i=1:obj.length
          %truncate
          mat_new{i}=mat_old{i}(1:l+1,1:l+1);
        end
      end
      %assign result
      obj.mat=mat_new;
    end
    function obj=set_lmax(obj,l)
      obj.lmax=l;
    end
    function out=get.lmax(obj)
      out=gravity.width2lmax(obj.width);
    end
    %% representations
    %returns a cell array with matrix representation
    function out=get.mat(obj)
      out=cell(obj.length,1);
      for i=1:obj.length
        out{i}=gravity.dtc('y','mat',obj.y(i,:));
      end
    end
    function obj=set.mat(obj,in)
      %sanity
      if ~iscell(in)
        error([mfilename,': input <in> must be a cell array of matrices.'])
      end
      if (numel(in)~=obj.length)
        error([mfilename,': cannot handle input <in> if it does not have the same number of elements as obj.length.'])
      end
      %make room for data
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('mat',in{1})));
      %build data
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(in{i});
      end
      %assign data
      obj=obj.assign(y_now,'reset_width',obj.width~=size(y_now,2));
    end
    %returns a structure array with C and S representation
    function out=get.cs(obj)
      out(obj.length)=struct('C',[],'S',[]);
      for i=1:obj.length
        out(i)=gravity.dtc('y','cs',obj.y(i,:));
      end
    end
    function obj=set.cs(obj,in)
      %sanity
      if ~isstruct(in)
        error([mfilename,': input <in> must be a structure array.'])
      end
      if (numel(in)~=obj.length)
        error([mfilename,': cannot handle input <in> if it does not have the same number of elements as obj.length.'])
      end
      %make room for data
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('cs',in)));
      %build data
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(gravity.cs2mat(in(i)));
      end
      %assign data
      obj=obj.assign(y_now);
    end
    %return a cell array with triangular matrix representation
    function out=get.tri(obj)
      out=cell(obj.length,1);
      for i=1:obj.length
        out{i}=gravity.dtc('y','tri',obj.y(i,:));
      end
    end
    function obj=set.tri(obj,in)
      %sanity
      if ~iscell(in) && ~ismatrix(in{1})
        error([mfilename,': input <in> must be a cell array of matrices.'])
      end
      if (numel(in)~=obj.length)
        error([mfilename,': cannot handle input <in> if it does not have the same number of elements as obj.length.'])
      end
      %make room for data
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('tri',in{1})));
      %build data
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(gravity.cs2mat(gravity.tri2cs(in{i})));
      end
      %assign data
      obj=obj.assign(y_now);
    end
    %return a cell array with the mod representation
    function out=get.mod(obj)
      out=cell(obj.length,1);
      for i=1:obj.length
        out{i}=gravity.dtc('y','mod',obj.y(i,:));
      end
    end
    function obj=set.mod(obj,in)
      %sanity
      if ~iscell(in) && ~ismatrix(in{1})
        error([mfilename,': input <in> must be a cell array of matrices.'])
      end
      if (numel(in)~=obj.length)
        error([mfilename,': cannot handle input <in> if it does not have the same number of elements as obj.length.'])
      end
      %make room for data
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('mod',in{1})));
      %build data
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(gravity.cs2mat(gravity.mod2cs(in{i})));
      end
      %assign data
      obj=obj.assign(y_now);
    end
    function out=obj2mod(obj)
      mod_list=obj.mod;
      out=cell(size(mod_list));
      for i=1:obj.length
      	out{i}=struct('mod',mod_list{i},'GM',obj.GM,'R',obj.R,'t',obj.t(i));
      end
    end
    %% checksum of the models
    function out=get.checksum(obj)
      out=obj.norm-sqrt(obj.C(0,0).^2);
    end
    %% info methods
    function print(obj,tab,lprint)
      if ~exist('lprint','var') || isempty(lprint)
        lprint=3;
      end
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      tri_now=obj.tri;
      n=min(obj.lmax,lprint)+1;
      disp(['Coefficients up to degree ',num2str(n-1),' of the first model, at ',datestr(obj.t(1))])
      format shortg
      disp(tri_now{1}(1:n,obj.lmax-n+2:obj.lmax+n))
      format short
      %parameters
      relevant_parameters={'GM','R','functional','origin','checksum'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print superclass
      print@simpletimeseries(obj,tab)
    end
    %% coefficient access
    function out=C(obj,d,o,time)
      if ~exist('time','var') || isempty(time)
        time=obj.t;
      end
      if any(size(d)~=size(o))
        error([mfilename,': inputs ''d'' and ''o'' must have the same size(s).'])
      end
      %get indexes of the cosine-coefficients
      C_idx=find(o>=0);
      S_idx=find(o<0);
      %make room for output
      out=cell(size(time));
      for i=1:numel(out)
        out{i}=zeros(size(d));
        for j=1:numel(C_idx)
          %retrive cosine coefficients
          out{i}( C_idx(j))=obj.mat{obj.idx(time(i))}( d( C_idx(j))+1,o( C_idx(j))+1);
        end
        for j=1:numel(S_idx)
          %skip zonals
          if -o(S_idx(j))==0; continue; end
          %retrieve sine coefficients
          out{i}( S_idx(j))=obj.mat{obj.idx(time(i))}(-o( S_idx(j))  ,d( S_idx(j))+1);
        end
      end
      %handle scalar quantities gracefully
      if isscalar(d)
        out_tmp=out;
        out=zeros(size(time));
        out(:)=[out_tmp{:}];
      end
    end
    function ts=ts_C(obj,d,o)
      idx=gravity.colidx(d,o,obj.lmax);
      if isempty(idx)
        ts=[];
        return
      end
      [labels,units]=gravity.labels(obj.lmax,obj.functional_unit);
      ts=simpletimeseries(...
        obj.t,...
        obj.y(:,idx),...
        'format','datetime',...
        'labels',labels(idx),...
        'timesystem',obj.timesystem,...
        'units',units(idx),...
        'descriptor',obj.descriptor...
      );    
    end
    function obj=setC(obj,d,o,values,time)
      if ~exist('time','var') || isempty(time)
        time=obj.t;
      end
      if any(size(d)~=size(o))
        error([mfilename,': inputs ''d'', ''o'' must have the same size'])
      end
      if numel(d)~=size(values,2) && numel(values)>1
        error([mfilename,': inputs ''d'' and ''o'' must have the same number of elements as the number of columns of ''values'''])
      end
      if numel(time)~=size(values,1) && numel(values)>1
        error([mfilename,': input ''time'' must have the same number of elements as the number of rows of ''values'''])
      end
      if numel(values)==1
        values=values*ones(numel(time),numel(d));
      end
      %interpolate over time if needed
      if numel(time)>1 && ~obj.istequal(time)
        values=interp1(time,values,obj.t_masked);
      end
      %retrieve matrix form
      mat_now=obj.mat;
      %get indexes of the cosine-coefficients
      for ti=1:numel(time)
        for i=1:numel(d)
          if o(i)>=0
            %set cosine coefficients
            mat_now{obj.idx(time(ti))}( d( i)+1,o( i)+1)=values(ti,i);
          else
            %retrieve sine coefficients
            mat_now{obj.idx(time(ti))}(-o( i)  ,d( i)+1)=values(ti,i);
          end
        end
      end
      %save new values
      obj.mat=mat_now;
    end
    function obj=setdegree(obj,d,values,time)
      if ~exist('time','var') || isempty(time)
        time=obj.t;
      end
      if any(size(d)~=size(values))
        error([mfilename,': inputs ''d'' and ''values'' must have the same size'])
      end
      %retrieve triangular form
      tri_now=obj.tri;
      %get indexes of the cosine-coefficients
      for ti=1:numel(time)
        for i=1:numel(d)
          %set cosine coefficients
          tri_now{obj.idx(time(ti))}(d(i)+1,:)=values( i);
        end
      end
      %save new values
      obj.tri=tri_now;
    end
    function obj=setorder(obj,o,values,time)
      if ~exist('time','var') || isempty(time)
        time=obj.t;
      end
      if any(size(o)~=size(values))
        error([mfilename,': inputs ''d'' and ''values'' must have the same size'])
      end
      %retrieve triangular form
      tri_now=obj.tri;
      %get indexes of the cosine-coefficients
      for ti=1:numel(time)
        for i=1:numel(o)
          %set cosine coefficients
          tri_now{obj.idx(time(ti))}(:,o(i)+1+obj.lmax)=values( i);
        end
      end
      %save new values
      obj.tri=tri_now;
    end
    %% scaling functions
    % GM scaling
    function s=scale_GM(obj,gm)
      s=gm/obj.GM;
    end
    % radius scaling
    function s=scale_R(obj,r)
      if obj.R==r
        s=1; %quicker downstream
      else
        s=(obj.R/r).^((0:obj.lmax)+1);
      end
    end
    % functional scaling
    function s=scale_functional(obj,functional_new)
      %   Supported functionals are the following,
      %       'non-dim'   - non-dimensional Stoked coefficients.
      %       'eqwh'      - equivalent water height [m]
      %       'geoid'     - geoid height [m]
      %       'potential' - [m2/s2]
      %       'gravity'   - [m /s2], if input represents the disturbing potential,
      %                     then the output represents the gravity disturbances.
      %                     Otherwise, it represents the gravity accelerations
      %                     (-dU/dr).
      %       'anomalies' - If input represents the disturbing potential, then
      %                     the output represents the gravity anomalies.
      %                     Otherwise, it represents (-dU/dr - 2/r*U).
      %       'vert-grav-grad' - vertical gravity gradient.
      %get nr of degrees
      N=obj.lmax+1;
      %get pre-scaling
      switch lower(obj.funct)
        case 'nondim'
          %no scaling
          pre_scale=ones(N);
        otherwise
          %need to bring these coefficients down to 'non-dim' scale
          functional_old=obj.funct;
          obj.funct='nondim';
          pre_scale = obj.scale_functional(functional_old);
      end
      %get pos-scaling
      switch lower(functional_new)
        case 'nondim'
          %no scaling
          pos_scale=ones(N);
        case 'eqwh' %[m]
          love=gravity.parameters('love');
          pos_scale=zeros(N);
          rho_earth=gravity.parameters('rho_earth');
          rho_water=gravity.parameters('rho_water');
          %converting Stokes coefficients from non-dimensional to equivalent water layer thickness
          for i=1:N
            deg=i-1;
            lv=interp1(love(:,1),love(:,2),deg,'linear','extrap');
            pos_scale(i,:)=obj.R*rho_earth/rho_water/3*(2*deg+1)/(1+lv);
          end
        case 'geoid' %[m]
          pos_scale=ones(N)*obj.R;
        case 'potential' %[m2/s2]
          pos_scale=ones(N)*obj.GM/obj.R;
        case 'gravity' %[m/s2]
          %If input represents the disturbing potential, then the output
          %represents the gravity disturbances. Otherwise, it represents the
          %gravity accelerations (-dU/dr).
          deg=(0:N-1)'*ones(1,N);
          pos_scale=-obj.GM/obj.R^2*(deg+1);
        case 'anomalies'
          %If input represents the disturbing potential, then the output
          %represents the gravity anomalies. Otherwise, it represents:
          %-dU/dr - 2/r*U
          deg=(0:N-1)'*ones(1,N);
          pos_scale=obj.GM/obj.R^2*max(deg-1,ones(size(deg)));
        case 'vertgravgrad'
          deg=(0:N-1)'*ones(1,N);
          pos_scale=obj.GM/obj.R^3*(deg+1).*(deg+2);
        otherwise
          error([mfilename,': unknown scale ',functional_new])
      end
      %outputs
      s=pos_scale./pre_scale;
    end
    % Gaussan smoothing scaling
    function s=scale_gauss(obj,fwhm_degree)
      % translate smoothing radius to degree (criteria inside) 
      fwhm_degree=gravity.gauss_smoothing_degree_translate(fwhm_degree);
      %NOTICE: gravity.gauss_smoothing_degree_translate(0) will assume the input is in degrees
      %        which should mean infinite smoothing, but generally the 0 will actually refer
      %        to a zero smoothing radius, i.e. no smoothing. This is the assumed convention.
      if fwhm_degree==0
        s=ones(1,obj.lmax+1);
      else
        %https://en.wikipedia.org/wiki/Gaussian_function#Properties
        s=exp(-4*log(2)*((0:obj.lmax)/fwhm_degree/2).^2);
      end
%       find(abs(s-0.5)==min(abs(s-0.5)))
%       %NOTICE: the tail of this function is unstable
%       %http://dx.doi.org/10.1029/98JB02844
%       b=log(2)/(1-cos(radius/obj.R));
%       c=exp(-2*b);
%       s=zeros(1,obj.lmax+1);
%       s(1)=1;                         %degree 0
%       s(2)=s(1)*((1+c)/(1-c) - 1/b);	%degree 1
%       for l=2:obj.lmax
%         s(l+1)=-(2*l+1)/b*s(l)+s(l-1);
%       end
    end
    function s=scale_spline(obj,fwhm_degree,width_degree)
      if ~exist('width_degree','var') || isempty(width_degree)
        width_degree=5;
      end
      % translate smoothing radius to degree (criteria inside) 
      fwhm_degree=gravity.gauss_smoothing_degree_translate(fwhm_degree);
      %NOTICE: usually smoothing up to degree zero would mean zeroing the data
      %        but the convention here is that is means no smoothing.
      if fwhm_degree==0
        s=ones(1,obj.lmax+1);
        return
      end
      if fwhm_degree<=width_degree
        width_degree=fwhm_degree-1;
      end
      s=nan(1,obj.lmax+1);
      w0=max([1,fwhm_degree-width_degree]);
      w1=min([obj.lmax,fwhm_degree+width_degree]);
      if w0>=1
        s(1:w0)=1;
      end
      if width_degree<=0
        s_transition=0.5;  
      else
        s_transition=spline([1 (2*width_degree+1)],[0 1 0 0],1:(2*width_degree+1));
      end
      s(w0+1:w1+1)=s_transition(1:w1-w0+1);
      if w1+2<=obj.lmax+1
        s(w1+2:obj.lmax+1)=0;
      end
    	%sanity
      assert(~any(isnan(s)),...
        [mfilename,': found NaNs in the output. Debug needed!'])
    end
    function s=scale_trunc(obj,fwhm_degree)
      % translate smoothing radius to degree (criteria inside) 
      fwhm_degree=gravity.gauss_smoothing_degree_translate(fwhm_degree);
      %https://en.wikipedia.org/wiki/Gaussian_function#Properties
      s=[ones(1,min([fwhm_degree+1,obj.lmax+1])),zeros(1,obj.lmax-fwhm_degree)];
%       find(abs(s-0.5)==min(abs(s-0.5)))
%       %NOTICE: the tail of this function is unstable
%       %http://dx.doi.org/10.1029/98JB02844
%       b=log(2)/(1-cos(radius/obj.R));
%       c=exp(-2*b);
%       s=zeros(1,obj.lmax+1);
%       s(1)=1;                         %degree 0
%       s(2)=s(1)*((1+c)/(1-c) - 1/b);	%degree 1
%       for l=2:obj.lmax
%         s(l+1)=-(2*l+1)/b*s(l)+s(l-1);
%       end
    end
    % scale according to number of orders in each degree
    function s=scale_nopd(obj,~)
      s=1./obj.nopd;
    end
    % trap for no scaling
    function s=scale_none(obj,~)
      s=ones(1,obj.lmax+1);
    end
    % scale operation agregator
    function out=scale_factor(obj,s,method)
      out=obj.(['scale_',method])(s);
    end
    function obj=scale(obj,s,method)
      if ~exist('method','var') || isempty(method)
        %handle object-to-object operations
        if isa(s,'gravity')
          %scale to object 's' values of R and GM
          obj=obj.setR(s.R);
          obj=obj.setGM(s.GM);
          %done!
          return
        end
        %trivial call
        if all(s==1)
          return
        end
        switch numel(s)
        case obj.lmax+1
          %get unit model, scaled per degree and get y-representation
          s=gravity.unit(obj.lmax,'scale_per_degree',s).y;
        case 1
          %global scaling: already handled downstream
        case obj.width
          %get CS representation
          c=obj.cs;
          %per-coefficients scaling: expand over all time domain
          for i=1:numel(c)
            c(i).C= c(i).C.*s;
            c(i).S= c(i).S.*s;
          end
          obj.cs=c;
          %we're done, bail now
          return
        otherwise
          error([mfilename,': cannot handle scaling factors with number of elements equal to ',...
            num2str(numel(s)),'; either 1, max degree+1 (',num2str(obj.lmax+1),') or nr of coeffs (',...
            num2str(obj.width),').'])
        end
        %call mother routine
        obj=scale@simpledata(obj,s);
      else
        % input 's' assumes different meanings, dependending on the method; invoke as:
        % obj.scale('geoid','functional')
        % obj.scale(300e3,'gauss')
        obj=obj.scale(obj.scale_factor(s,method));
        %need to update metadata in some cases
        switch lower(method)
          case 'gm'
            obj.GM=s;
          case 'r'
            obj.R=s;
          case 'functional'
            obj.funct=s;
            obj.units(:)={gravity.functional_units(s)};
        end
      end
    end
    % scale plotting
    function out=scale_plot(obj,s,method)
      out.axishandle=plot(0:obj.lmax,obj.scale_factor(s,method),'LineWidth',2);
      hold on
      grid on
      xlabel('SH degree')
      ylabel('scale factor [ ]')
      out.legend_str=method;
    end
    %% GM and R operations
    function out=getGM(obj)
      out=obj.GM;
    end
    function obj=setGM(obj,gm)
      obj=obj.scale(gm,'GM');
    end
    function out=getR(obj)
      out=obj.R;
    end
    function obj=setR(obj,r)
      obj=obj.scale(r,'R');
    end
    %% derived quantities
    % number of orders in each degree
    function out=nopd(obj)
      out=2*(1:obj.lmax+1)-1;
    end
    % aux function to create time series from derived quantities
    function ts=makets(obj,d,l,title,label)
      %if data is only a column, then there is no degree dependence
      if size(d,2)==1
        ts=simpletimeseries(...
          obj.t,...
          d,...
          'format','datetime',...
          'labels',{label},...
          'timesystem',obj.timesystem,...
          'units',{obj.functional_unit},...
          'descriptor',[title,' @ degree ',num2str(l),' of ',obj.descriptor]...
        );    
      else
        assert(size(d,2)==l+1,['Input ''d'' must have ',num2str(l+1),' columns, not ',num2str(size(d,2))])
        ts=simpletimeseries(...
          obj.t,...
          d,...
          'format','datetime',...
          'labels',cellfun(@(i) ['deg. ',i],strsplit(num2str(0:l)),'UniformOutput',false),...
          'timesystem',obj.timesystem,...
          'units',repmat({obj.functional_unit},1,l+1),...
          'descriptor',[title,' of ',obj.descriptor]...
        );
      end
    end
    %degree mean
    function [d,ts,title,label]=dmean(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      d=zeros(obj.length,obj.lmax+1);
      tri_now=obj.tri;
      m=transpose(obj.nopd);
      for i=1:obj.length
        %compute mean over each degree (don't use 'mean' here, there's lots of zeros in tri_now for the low degrees)
        d(i,:) = sum(tri_now{i},2)./m;
      end
      if nargout>1
        title='Degree mean';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    %degree STD
    function [d,ts,title,label]=dstd(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      d=sqrt(obj.drms(l).^2-obj.dmean(l).^2);
      if nargout>1
        title='Degree STD';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    %degree RMS
    function [d,ts,title,label]=drms(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      das=obj.das(l);
      d=zeros(size(das));
      m=obj.nopd;
      for i=1:obj.length
        d(i,:)=das(i,:)./sqrt(m);
      end
      if nargout>1
        title='Degree RMS';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    % returns degree amplitude spectrum for each row of obj.y. The output
    % matrix has in each row the epochs of obj.y (corresponding to the epochs
    % of the models) and for each column i the degree i-1 (this is implicit).
    function [d,ts,title,label]=das(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      d=zeros(obj.length,l+1);
      tri_now=obj.tri;
      for i=1:obj.length
        %compute Degree amplitude
        d(i,:) = sqrt(sum(tri_now{i}(1:l+1,:).^2,2));
      end
      if nargout>1
        title='Degree amplitude';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    %cumulative degree mean
    function [d,ts]=cumdmean(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      d=sqrt(cumsum(obj.dmean(l).^2,2));
      if nargout>1
        title='Cum. degree mean';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    %cumulative degree STD
    function [d,ts,title,label]=cumdstd(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      d=sqrt(cumsum(obj.dstd(l).^2,2));
      if nargout>1
        title='Cum. degree STD';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    %cumulative degree RMS
    function [d,ts]=cumdrms(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      d=sqrt(cumsum(obj.drms(l).^2,2));
      if nargout>1
        title='Cum. degree RMS';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    % cumulative degree amplitude spectrum
    function [d,ts,title,label]=cumdas(obj,l)
      if ~exist('l','var') || isempty(l)
        l=obj.lmax;
      end
      %NOTICE: this is squared-sum-rooted, not simply summed!
      d=sqrt(cumsum(obj.das.^2,2));
      if nargout>1
        title='Cum. degree amplitude';
        label=[title,' @ degree ',num2str(l)];
        ts=obj.makets(d(:,end),l,title,label);
      end
    end
    % created a timeseries object with the derived quantities above, *for all degrees* (unlike the above which is only for one degrees)
    function ts=derived(obj,quantity)
      %NOTICE: there is no 'l' argument in this call so that obj.lmax is used
      [d,~,title,label]=obj.(quantity);
      ts=obj.makets(d,obj.lmax,title,label);
    end
    % gridded wrapers for statistics (the mean and std are not possible with SH)
    function [d,ts,title,label]=gridmean(obj)
      [d,~,ts]=obj.grid.mean('total');
      title='Grid mean';
      label=title;
    end
    function [d,ts,title,label]=gridrms(obj)
      [d,~,ts]=obj.grid.rms('total');
      title='Grid RMS';
      label=title;
    end
    function [d,ts,title,label]=gridstd(obj)
      [d,~,ts]=obj.grid.std('total');
      title='Grid STD';
      label=title;
    end
    %% convert to grid
    function out=grid(obj,varargin)
      v=varargs.wrap('sources',{{...
        'Nlat', obj.lmax+1    , @num.isscalar;...
        'Nlon',(obj.lmax+1)*2 , @num.isscalar;...
        'Nfactor',          1 , @num.isscalar;...
        'spatial_step',   inf , @num.isscalar;...
      }},varargin{:});
      %easier names for the resolution
      Nlat=max([v.Nlat,simplegrid.lat_stepped_length(v.spatial_step)])*v.Nfactor;
      Nlon=max([v.Nlon,simplegrid.lon_stepped_length(v.spatial_step)])*v.Nfactor;
      %retrieve CS-representation
      cs_now=obj.cs;
      %initiate grid 3d-matrix
      map=nan(Nlat,Nlon,obj.length);
      %get default latitude domain
      lat=simplegrid.lat_default(Nlat);
      %make synthesis on each set of coefficients
      for i=1:obj.length
        [lon,~,map(:,:,i)]=mod_sh_synth(cs_now(i).C,cs_now(i).S,deg2rad(lat),Nlon);
      end
      %create grid structure in the form of a list
      sl=simplegrid.dti(obj.t,map,transpose(rad2deg(lon(:))),lat(:),'list');
      %create container for units
      units=cell(1,size(sl.map,2));units(:)={obj.functional_unit};
      %initialize grid object
      out=simplegrid(sl.t,sl.map,...
        'lon',sl.lon,...
        'lat',sl.lat,...
        'descriptor',obj.descriptor,...
        'units',units...
      );
    end
    function out=point(obj,lon,lat,h)
      assert(all(size(lat)==size(lon)) && all(size(lon)==size(h)),'all inputs must have the same length')
      assert(isvector(lat),'all inputs must be vectors')
      %wrap inputs
      lon=wrapTo360(lon);
      %make room for outputs
      out=zeros(obj.length,numel(lat));
      %determine number of longitude points
      if numel(lon)>1 && mean(diff(lon))>0
        Nlon=360/mean(diff(lon));
      else
        Nlon=360;
      end
      %retrieve CS-representation
      cs_now=obj.cs;
      %make synthesis on each set of coefficients
      for i=1:obj.length
        for j=1:numel(lat)
          %scale to requested height
          obj=obj.scale(obj.R+h(j),'R');
          %perform synthesis
          [lon_out,lat_out,tmp]=mod_sh_synth(cs_now(i).C,cs_now(i).S,deg2rad(lat(j)),Nlon);
          %sanity (probably unnecessary)
          assert(abs(lat(j)-rad2deg(lat_out))<1e-12,'discrepancy in input/output lat')
          %retrieve longitude
          lon_out=rad2deg(lon_out);
          if any(lon_out==lon(j))
            out(i,j)=tmp(lon_out==lon(j));
          else
            out(i,j)=interp1(lon_out,tmp,lon(j));
          end
        end
      end
    end
    %% spatial operations
    function obj=spatial_mask(obj,mode,varargin)
      obj=obj.grid('Nfactor',3,varargin{:}).spatial_mask(mode,varargin{:}).sh(obj.lmax,'functional',obj.functional);
    end
    %% overloading
    function [obj,stats]=component_split(obj,varargin)
      v=varargs.wrap('sources',{{...
        't_extrapolate',  obj.stop,   @isdatetime;...
      }},varargin{:});
      %get collection of degrees/orders
      [degrees,orders]=gravity.resolve_degrees_orders('lmax',obj.lmax,varargin{:});
      %make room for optional outputs
      stats=struct(...
        'degrees',degrees,...
        'orders',orders,...
        'idx',gravity.colidx(degrees,orders,obj.lmax),...
        'dat',{cell(1,obj.width)}...
      );
      %prepare for parloop
      ts=cell(size(degrees));
      for i=1:numel(degrees)
        %get time series of this coefficient
        ts{i}=obj.ts_C(degrees(i),orders(i));
      end
      dat=cell(size(degrees));
      out=cell(size(degrees));
      t_extrapolate=v.t_extrapolate;
      %operate over each coefficient individually
      parfor i=1:numel(degrees)
        %inform
        str.say('Looking at degree ',degrees(i),' order ',orders(i))
        %skip constant coefficients (usually zero)
        if all(ts{i}.y_masked(1)==ts{i}.y_masked)
          %extend time domain
          ts{i}=ts{i}.append_epochs( (ts{i}.stop+ts{i}.step):ts{i}.step:t_extrapolate , ts{i}.y_masked(1) );
        else
          %call mother routine, save stats to index consistent with columns of y
          [ts{i},dat{i}]=ts{i}.component_split(varargin{:}); %#ok<PFBNS>
        end
      end
      %init outputs
      stats.dat=cell(size(dat));
      %collect
      for i=1:numel(degrees)
        d=degrees(i);
        o=orders(i);
        %initialize
        if i==1
          out=gravity.unit(obj.lmax,'scale',0,'t',ts{i}.t);
          out=out.copy_metadata(obj);
        end
        %accumulate
        out=out.setC(d,o,ts{i}.y);
        %get column index
        colidx=gravity.colidx(d,o,obj.lmax);
        %propagate stats
        stats.dat(colidx)=dat(i);
      end
      %propagate
      obj=out;
    end
    %% multiple operands
    function compatible(obj1,obj2,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('skip_par_check',{''},@iscellstr)
      p.parse(varargin{:});
      %call mother routine
      compatible@simpletimeseries(obj1,obj2,varargin{:});
      %shorter names
      par=gravity.compatible_parameter_list;
      for i=1:numel(par)
        % if a parameter is empty, no need to check it
        if ( iscell(obj1.(par{i})) && isempty([obj1.(par{i}){:}]) ) || ...
           ( ischar(obj1.(par{i})) && isempty( obj1.(par{i})    ) ) || ...
           ( iscell(obj2.(par{i})) && isempty([obj2.(par{i}){:}]) ) || ...
           ( ischar(obj2.(par{i})) && isempty( obj2.(par{i})    ) )
          continue
        end
        if ~cells.isincluded(p.Results.skip_par_check,par{i}) && ~isequal(obj1.(par{i}),obj2.(par{i}))
          error([mfilename,': discrepancy in parameter ',par{i},'.'])
        end
      end
    end
    function [obj1,obj2,idx1,idx2]=merge(obj1,obj2,varargin)
      %match the minimum degree (truncate)
      if obj1.lmax<obj2.lmax
        obj2.lmax=obj1.lmax;
      else
        obj1.lmax=obj2.lmax;
      end
      %match the compatible parameters
      parameters=gravity.compatible_parameter_list;
      for i=1:numel(parameters)
        p=(parameters{i});
        if ~isequal(obj1.(p),obj2.(p))
          obj2=obj2.scale(obj1.(p),p);
        end
      end
      %extend the time-domain of both objects to be in agreement with the each other.
      [obj1,obj2,idx1,idx2]=merge@simpletimeseries(obj1,obj2);
    end
    function obj1=glue(obj1,obj2,varargin)
      %objects need to have the same time domain
      assert(obj1.istequal(obj2),'Input objects do not share the same time domain.')
      %get maximum degree
      lmax_now=max([obj1.lmax,obj2.lmax]);
      %match lmax
      obj1.lmax=lmax_now; obj2.lmax=lmax_now;
      %make sure objects are compatible
      compatible(obj1,obj2,varargin{:})
      %get mapping
      map=gravity.mapping(lmax_now);
      %augment the data, labels and units, coefficient-wise
      for i=1:size(map,2)
        y1=obj1.y_masked([],i);
        y2=obj2.y_masked([],i);
        if     any(y1~=0) && all(y2==0)
          %do nothing, obj1 has the correct data at this column
        elseif all(y1==0) && any(y2~=0)
          %propagate this column from obj2 to obj1
          obj1=obj1.set_cols(i,obj2.y(:,i));
        elseif all(y1==0) && all(y2==0)
          %do nothing, obj1 has already zeros at this column
        else
          %ensure data is consistent
          assert(all(y1==y2),['Cannot glue objects at degree ',num2str(map(1,i)),' and  order ',num2str(map(2,i)),...
            'because data is not consistent.'])
        end
      end
    end
    %% plot functions
    function out=plot(obj,varargin)
      v=varargs.wrap('sources',{{...
        'method',   'drms',    @ischar;...
        'showlegend',false,    @islogical;...
        'line',     '-',       @ischar;...
        'title',    '',        @ischar;...
        'functional',obj.funct,@ischar;...
        'time',      [],       @(i) simpletimeseries.valid_t(i) || isempty(i);...
        'colormap',  '',       @(i) ischar(i) || ismatrix(i);...
        'degrees',[2,3],       @isnumeric;...
        'orders', [0,0],       @isnumeric;...
      }},varargin{:}); 
      % enforce requested functional
      if ~strcmpi(obj.funct,v.functional)
        obj=obj.scale(v.functional,'functional');
      end
      %consider only requested epoch(s)
      if ~isempty(v.time)
        obj=obj.at(v.time);
      end
      %if data columns are specified, then fallback to timeseries plotting
      if cells.isincluded(varargin,'columns')
        v.method='timeseries';
      end
      %branch on method
      switch lower(v.method)
      case 'timeseries'
        %call superclass
        out=plot@simpletimeseries(obj,'columns',gravity.colidx(v.degrees,v.orders,obj.lmax),varargin{:});
      case {'cumdmean-timeseries','cumdrms-timeseries','cumdstd-timeseries','cumdas-timeseries'}
        %compute cumdrms for all epochs
        out=obj.(strrep(v.method,'-timeseries','')).plot(varargin{:});
      case {'triang','trianglog10'}
        %get triangular plots (don't plot invalid entries)
        obj_masked=obj.masked;
        if isempty(obj_masked); out=[];return; end
        tri_now=obj_masked.tri;
        for i=1:numel(tri_now)
          bad_idx=(tri_now{i}==0);
          tri_now{i}(bad_idx)=NaN;
          if i>1;figure;end
          if strcmpi(v.method,'trianglog10')
            out.image=imagesc([-obj.lmax,obj.lmax],[0,obj.lmax],log10(abs(tri_now{i})));
          else
            out.image=imagesc([-obj.lmax,obj.lmax],[0,obj.lmax],tri_now{i});
          end
          out.ylabel='SH degree';
          out.xlabel='SH order';
          ylabel(out.ylabel)
          xlabel(out.xlabel)
          %TODO: merge this with a to-be-implemented plotting.colormap
          if ~isempty(v.colormap)
            colormap(gca,v.colormap)
          end
          cb.nan;
          out.cb=colorbar;
          if strcmpi(v.method,'trianglog10')
            cb.label(['log_{10}(',obj.functional_name,') [',obj.functional_unit,']']);
          else
            cb.label([obj.functional_name,' [',obj.functional_unit,']']);
          end
          if isempty(v.title)
            out.title=[obj.descriptor,' - ',datestr(obj.t(i)),', \mu=',num2str(mean(tri_now{i}(~bad_idx)))];
          else
            out.title=v.title;
          end
          title(out.title)
          out.axis_handle=gca;
          out.fig_handle=gcf;
        end
      case {'dmean','cumdmean','drms','cumdrms','dstd','cumdstd','das','cumdas'}
        d=transpose(obj.(v.method));
        switch lower(v.method)
        case {'dmean','cumdmean'}
          out.line=semilogy(0:obj.lmax,abs(d),v.line);hold on
        otherwise
          out.line=semilogy(0:obj.lmax,d,v.line);hold on
        end
        grid on
        xlabel('SH degree')
        ylabel([obj.functional_name,' [',obj.functional_unit,']'])
        if v.showlegend
          out.legend=legend(datestr(obj.t));
        end
        %title: append functional if no title specified
        if isempty(v.title)
          switch lower(v.method)
            case 'dmean',    title_str='degree mean';
            case 'cumdmean', title_str='cumul. degree mean';
            case 'drms',     title_str='degree RMS';
            case 'cumdrms',  title_str='cumul. degree RMS';
            case 'dstd',     title_str='degree STD';
            case 'cumdstd',  title_str='cumul. degree STD';
            case 'das',      title_str='degree amplit.';
            case 'cumdas',   title_str='cumul. degree amplit.';
          end
          out.title=[obj.descriptor,' - ',title_str];
        else
          out.title=v.title;
        end
        title(str.clean(out.title,'title'))
      otherwise
        error([mfilename,': unknonw method ''',v.method,'''.'])
      end
    end
    %% export functions
    function filelist=icgem(obj,varargin)
      % Parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      % optional arguments
      p.addParameter('prefix',   '',  @ischar);
      p.addParameter('suffix',   '',  @ischar);
      p.addParameter('delim',    '.', @ischar);
      p.addParameter('path',     '.', @ischar);
      p.addParameter('timefmt',  'yyyymmddTHHMMSS', @ischar);
      p.addParameter('error_obj','',  @(i)isa(i,'gravity') || isempty(i));
      p.addParameter('modelname','unknown',         @ischar);
      p.addParameter('errors',   'formal',          @ischar);
      p.addParameter('norm',     'fully_normalized',@ischar);
      % parse it
      p.parse(varargin{:});
      % get only valid entries
      obj_now=obj.masked;
      % init outputs
      filelist=cell(obj.length,1);
      % build filenames
      for i=1:obj_now.length
        f={p.Results.prefix,datestr(obj_now.t(i),p.Results.timefmt),p.Results.suffix};
        f=strjoin(cells.rm_empty(f),p.Results.delim);
        filelist{i}=fullfile(p.Results.path,[f,'.gfc']);
      end
      header={...
       ''... %place holder for model epoch
       ''... %place holder for time of exporting
       ['timesystem                  ',obj_now.timesystem],...
       ['descriptor                  ',obj_now.descriptor],...
       ['geopotential_functional     ',obj_now.functional],...
        'begin_of_head ===========================================================================',...
        'product_type                gravity_field',...
       ['modelname                   ',p.Results.modelname],...
       ['earth_gravity_constant      ',num2str(obj_now.GM,'%.15g')],...
       ['radius                      ',num2str(obj_now.R, '%.15g')],...
        'min_degree                  0',...
       ['max_degree                  ',num2str(obj_now.lmax,'%d')],...
        'tide_system                 zero_tide',...
       ['modelname                   ',p.Results.errors],...
       ['norm                        ',p.Results.norm],...
        'key   n   m             C                      S                 sigma C       sigma S',...
        'end_of_head ============================================================================='...
      };
      %get the data in a suitable representations
      sig=obj_now.mod;
      % handle empty error models
      if isempty(p.Results.error_obj)
        err=cell(size(sig));
        err(:)=gravity.unit(obj.lmax,'scale',0).mod;
      else
        err=p.Results.error_obj.mod;
      end
      % write data
      for i=1:obj_now.length
        %check if file is already available
        if ~exist(filelist{i},'file')
          %open file
          fid=file.open(filelist{i},'w');
          %update the header
          header{1}=['model_epoch                 ',datestr(obj_now.t(i),   p.Results.timefmt)];
          header{2}=['exported_at                 ',datestr(datetime('now'),p.Results.timefmt)];
          %write the header
          fprintf(fid,'%s\n',strjoin(header,'\n'));
          %paranoid sanity
          assert(all(sig{i}(:,1)==err{i}(:,1)) && all(sig{i}(:,2)==err{i}(:,2)),...
            'degree/order domain discrepancy between signal and error')
          %write the data
          for j=1:numel(sig{i}(:,1))
            fprintf(fid,'gfc %3d %3d %23.15e %23.15e %14.8g %14.8g\n',...
              [sig{i}(j,1),sig{i}(j,2),sig{i}(j,3),sig{i}(j,4),err{i}(j,3),err{i}(j,4)]...
            );
          end
          %close the file
          fclose(fid);
          disp(['Finished writing:  ',filelist{i}])
        else
          disp(['Already available: ',filelist{i}])
        end
      end
    end
    function filelist=csr(obj,varargin)
      % Parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      % optional arguments
      p.addParameter('prefix',   '',  @ischar);
      p.addParameter('suffix',   '',  @ischar);
      p.addParameter('delim',    '.', @ischar);
      p.addParameter('path',     '.', @ischar);
      p.addParameter('timefmt',  'yyyymmddTHHMMSS', @ischar);
      p.addParameter('error_obj','',  @(i)isa(i,'gravity') || isempty(i));
      % parse it
      p.parse(varargin{:});
      % get only valid entries
      obj_now=obj.masked;
      % init outputs
      filelist=cell(obj.length,1);
      % build filenames
      for i=1:obj_now.length
        f={p.Results.prefix,datestr(obj_now.t(i),p.Results.timefmt),p.Results.suffix};
        f=strjoin(cells.rm_empty(f),p.Results.delim);
        filelist{i}=fullfile(p.Results.path,[f,'.GEO']);
      end
      %build header
      header={...
        '(2a10,2e20.13)';...
        [' SOLUTION  FIELD     ',num2str(obj.GM,'%20.13E'),' ',num2str(obj.R,'%20.13E')];...
        '(A6,2I3,2D21.14,2D11.4,F4.1)';...
      };...
      %get the data in a suitable representations
      sig=obj_now.mod;
      % handle empty error models
      if isempty(p.Results.error_obj)
        err=cell(size(sig));
        err(:)=gravity.unit(obj.lmax,'scale',0).mod;
      else
        err=p.Results.error_obj.mod;
      end
      % write data
      for i=1:obj_now.length
        %check if file is already available
        if ~exist(filelist{i},'file')
          %open file
          fid=file.open(filelist{i},'w');
          %write the header
          fprintf(fid,'%s\n',strjoin(header,'\n'));
          %paranoid sanity
          assert(all(sig{i}(:,1)==err{i}(:,1)) && all(sig{i}(:,2)==err{i}(:,2)),...
            'degree/order domain discrepancy between signal and error')
          %write the data
          for j=1:numel(sig{i}(:,1))
            fprintf(fid,'%6s%3d%3d%21.14e%21.14e%11.4e%11.4e%4.1f\n',...
              'RECOEF',sig{i}(j,1),sig{i}(j,2),sig{i}(j,3),sig{i}(j,4),err{i}(j,3),err{i}(j,4),-1 ...
            );
          end
          %close the file
          fclose(fid);
          disp(['Finished writing:  ',filelist{i}])
        else
          disp(['Already available: ',filelist{i}])
        end
      end
      %
    end
  end
end

%% load interfaces
function [m,e]=load_gsm(filename,time)
  %default time
  if ~exist('time','var') || isempty(time)
    time=datetime('now');
  end
  %open file
  fid=file.open(filename);
  modelname=''; GM=0; radius=0; Lmax=0; %Mmax=0;
  % Read header
  s=fgets(fid); c=0;
  while(strncmp(s, 'GRCOF2', 6) == 0)
     if (keyword_search(s, 'FIRST'))
       modelname =strtrim(s(7:end));
     end
     if (keyword_search(s, 'EARTH'))
%        x=str2num(s(7:end));
       x=str.num(s(7:end));
       GM=x(1);
       radius=x(2);
     end
     if keyword_search(s, 'SHM') && ~keyword_search(s, 'SHM*')
%        x=str2num(s(4:16));
       x=str.num(s(4:16));
       Lmax=x(1);
       %Mmax=x(2);
       permanent_tide=strfind(s,'inclusive permanent tide');
     end
     %handle yaml headers
     if str.contains(s,'End of YAML header')
       %need to read YAML (load packages dir, inside which the +yaml package sits)
       addpath(fullfile(file.orbdir('packages')));
       %rewind
       frewind(fid)
       %build header strings
       yaml_header=cell(c,1);
       for i=1:c
         yaml_header{i}=fgets(fid);
       end
       %define header filename
       [~,header_filename]=fileparts(filename);
       header_filename=fullfile('tmp',[header_filename,'.yaml']);
       %write to file
       file.strsave(header_filename,strjoin(yaml_header,'\n'));
       %read yaml from file (that's how it works...)
       d=yaml.ReadYaml(header_filename);
       %delete file
       delete(header_filename)
       %translate header info
       modelname=d.header.global_attributes.title;
       GM=d.header.non0x2Dstandard_attributes.earth_gravity_param.value;
       radius=d.header.non0x2Dstandard_attributes.mean_equator_radius.value;
       permanent_tide=strfind(d.header.non0x2Dstandard_attributes.permanent_tide_flag,'inclusive');
       %read end of header string
       fgets(fid);
     end
     %increment counter
     c=c+1;
     %get next line
     s=fgets(fid);
  end
  %sanity
  if sum(s)<0
     error([mfilename,': Problem with reading the GRCOF2 file ''',filename,'''.'])
  end
  %make room for coefficients
  mi.C=zeros(Lmax+1);
  mi.S=zeros(Lmax+1);
  ei.C=zeros(Lmax+1);
  ei.S=zeros(Lmax+1);
  % read coefficients
  while (s>=0)
    %skip empty lines
    if numel(s)<5
      s=fgets(fid);
      continue
    end
    %get numeric data
%     x=str2num(s(7:76)); %#ok<*ST2NM>
    x=str.num(s(7:76));
    %save degree and order
    d=x(1)+1;
    o=x(2)+1;
    %check if this is a valid line
    if strcmp(s(1:6),'GRCOF2')
      mi.C(d,o)=x(3);
      mi.S(d,o)=x(4);
      if (numel(x)>=6)
         ei.C(d,o)=x(5);
         ei.S(d,o)=x(6);
      end
    else
      error([mfilename,': unexpected tag in line: ''',s,''.'])
    end
    % read next line
    s=fgets(fid);
  end
  %fix permanent_tide
  if ~permanent_tide
    mi.C(3,1)=gravity.zero_tide(mi.C(3,1),'tide_free');
  end
  %initializing data object
  m=gravity(...
    time,...
    gravity.dtc('cs','y',mi),...
    'GM',GM,...
    'R',radius,...
    'descriptor',modelname,...
    'cdate',datetime('now'),...
    'origin',filename...
  );
  %initializing error object
  e=gravity(...
    time,...
    gravity.dtc('cs','y',ei),...
    'GM',GM,...
    'R',radius,...
    'descriptor',['error of ',modelname],...
    'cdate',datetime('now'),...
    'origin',filename...
  );
end
function [m,e]=load_csr(filename,time)
  %default time
  if ~exist('time','var') || isempty(time)
    time=datetime('now');
  end
  %open file
  fid=file.open(filename);
  modelname=''; Lmax=360; %Mmax=0;
  % Read header
  % skip first line (with a fortran format specifier)
  fgets(fid);
  %get GM and radius from the second line
  s=fgets(fid);
      GM=str.num(s(21:40));
  radius=str.num(s(41:60));
  %sanity
  if GM==0 || radius==0
    error([mfilename,': Problem with reading the CSR file ''',filename,''', could not find GM and R constants.'])
  end
  if sum(s)<0
    error([mfilename,': Problem with reading the CSR file ''',filename,'''.'])
  end
  % skip third line
  fgets(fid);
  %make room for coefficients
  mi.C=zeros(Lmax+1);
  mi.S=zeros(Lmax+1);
  ei.C=zeros(Lmax+1);
  ei.S=zeros(Lmax+1);
  Lmax=0;
  % get first line with coefficients
  s=fgets(fid);
  % read coefficients
  while (s>=0)
%     %skip irrelevant lines
%     if strncmp(s, 'RECOEF', 6) == 0 && strncmp(s, 'SUMGEO', 6) == 0
%       s=fgets(fid);
%       continue
%     end
    if length(s)==81
      %save degree and order
      d=str.num(s(7:9))+1;
      if d>Lmax;Lmax=d;end
      o=str.num(s(10:12))+1;
      mi.C(d,o)=str.num(s(13:33));
      mi.S(d,o)=str.num(s(34:54));
      ei.C(d,o)=str.num(s(55:65));
      try 
        ei.S(d,o)=str.num(s(66:76));
      catch
        ei.S(d,o)=str.num(s(66:74));
      end
    else
      disp(['Skipped line: ''',s,'''.'])
    end
    % read next line
    s=fgets(fid);
  end
  %crop containers
  mi.C=mi.C(1:Lmax,1:Lmax);
  mi.S=mi.S(1:Lmax,1:Lmax);
  ei.C=ei.C(1:Lmax,1:Lmax);
  ei.S=ei.S(1:Lmax,1:Lmax);
%   %fix permanent_tide
%   if permanent_tide
%     mi.C(3,1)=mi.C(3,1)-4.173e-9;
%   end
  %initializing data object
  m=gravity(...
    time,...
    gravity.dtc('cs','y',mi),...
    'GM',GM,...
    'R',radius,...
    'descriptor',modelname,...
    'cdate',datetime('now'),...
    'origin',filename...
  );
  %initializing error object
  e=gravity(...
    time,...
    gravity.dtc('cs','y',ei),...
    'GM',GM,...
    'R',radius,...
    'descriptor',['error of ',modelname],...
    'cdate',datetime('now'),...
    'origin',filename...
  );
end
function [m,e,trnd,acos,asin]=load_icgem(filename,time)
%This function is an adaptation of icgem2mat.m from rotating_3d_globe, by
%Ales Bezdek, which can be found at:
%
%http://www.asu.cas.cz/~bezdek/vyzkum/rotating_3d_globe/
%
%The original header is transcribed below.
%
%J.Encarnacao (j.g.deteixeiradaencarnacao@tudelft.nl) 11/2013
%
% ICGEM2MAT   Reads geopotential coefficients from an ICGEM file and saves them in a mat file.
%
% Usage:
%
%       icgem2mat
%
% finds all the ICGEM files (*.gfc) in the current directory,
% reads the geopotential coefficients, transforms them into Matlab variables:
%       header...structure with Icgem header information
%       mi.C(n+1,m+1), mi.S(n+1,m+1)...harmonic coefficients C(n,m), S(n,m)
%
% The new mat file with the same name is moved into 'data_icgem' subdirectory;
% the original gfc file is moved into 'data_icgem/gfc/' subdirectory.
%
% Add the 'data_icgem' folder into your Matlab path.
% The model coefficients are then loaded by typing, e.g.:
%
%          load egm2008
%
% To display the C(2,0) zonal term type
%
%          mi.C(3,1)
%
%
% See also compute_geopot_grids
%
% Ales Bezdek, bezdek@asu.cas.cz, 11/2012
%
% clear
% NMAX=360;
% NMAX=1e100;  %it is possible to limit the maximum degree read from the gfc file
% adr_data='./';
% adr_kam='./data_icgem/';
%
% seznam_soub=dir(adr_data);
% soub={seznam_soub.name};   %cell with filenames
% for i=1:length(soub)
%    jm=soub{i};
%    if length(jm)>4 && strcmpi(jm(end-3:end),'.gfc')
%       soub1=jm(1:end-4);
%       fprintf('Gfc file processed: %s\n',soub1);
%       filename=[adr_data soub1 '.gfc'];

  %default time
  if ~exist('time','var') || isempty(time)
    time=datetime('now');
  end
  %open file
  fid=file.open(filename);
  % init header
  header=struct(...
      'product_type',           '',...
      'modelname',              '',...
      'model_content',          '',...
      'earth_gravity_constant', [],...
      'radius',                 [],...
      'max_degree',             [],...
      'errors',                 '',...
      'norm',                   '',...
      'tide_system',            '',...
      'filename',               filename...
  );
  % Read header
  s=fgets(fid);
  fn=fieldnames(header);
  while(strncmp(s, 'end_of_head', 11) == 0 && sum(s)>=0)
    for j=1:numel(fn)
      f=fn{j};
      if (keyword_search(s,f))
        valuestr=strtrim(s(length(f)+1:end));
        switch class(header.(f))
          case 'double'
            header.(f)=str2double(strrep(valuestr,'D','e'));
          case 'char'
            header.(f)=valuestr;
          otherwise
          error([mfilename,': cannot handle class ',class(header.(f)),'.'])
        end
      end
    end
    s=fgets(fid);
  end
  if s<0
    error([mfilename,'Problem reading the gfc file.'])
  end
  % sanity on max degree
  if isempty(header.max_degree)
    error([mfilename,': could not determine maximum degree of model ''',filename,'''.'])
  end
  % make room for coefficients
  mi=struct('C',zeros(header.max_degree+1),'S',zeros(header.max_degree+1),'t0',[]);
  ei=struct('C',zeros(header.max_degree+1),'S',zeros(header.max_degree+1));
  trnd=struct('C',[],'S',[]);
  acos=struct('C',[],'S',[]);
  asin=struct('C',[],'S',[]);
  %iterators
  i_t0=0;
  i_trnd=0; %pocet clenu s trendem
  i_acos=0; %pocet clenu
  i_asin=0; %pocet clenu
  i_gfc=0;
  %read data
  s=fgets(fid);
  while (s>=0)
    %skip empty lines
    if numel(s)<5
      s=fgets(fid);
      continue
    end
    %retrieve keywords, degree and order
%     x=str2num(s(5:end));
    x=str.num(s(5:end));
    n=x(1)+1;
    m=x(2)+1;
    if strcmp(s(1:4),'gfct')
      i_t0=i_t0+1;
      mi.C(n,m)=x(3);
      mi.S(n,m)=x(4);
      [yr,mn,dy]=ymd2cal(x(end)/1e4);
      yrd=jd2yr(cal2jd(yr,mn,dy));
      if isempty(mi.t0)
        mi.t0=zeros(grep_nr_occurences(filename,'gfct'));
        if numel(mi.t0)==0; error([mfilename,'Problem with gfct']); end
      end
      mi.t0(i_t0,:)=[n m yrd];
      if (strcmp(header.errors,'formal') || ...
          strcmp(header.errors,'calibrated') || ...
          strcmp(header.errors,'calibrated_and_formal'))
        ei.C(n,m)=x(5);
        ei.S(n,m)=x(6);
      end
    elseif strcmp(s(1:3),'gfc')
      mi.C(n,m)=x(3);
      mi.S(n,m)=x(4);
      if (strcmp(header.errors,'formal') || ...
          strcmp(header.errors,'calibrated') || ...
          strcmp(header.errors,'calibrated_and_formal'))
        ei.C(n,m)=x(5);
        ei.S(n,m)=x(6);
      end
      i_gfc=i_gfc+1;
    elseif strcmp(s(1:4),'trnd') || strcmp(s(1:3),'dot')
      if isempty(trnd.C)
        trnd.C=zeros(grep_nr_occurences(filename,'trnd')+grep_nr_occurences(filename,'dot'),3); trnd.S=trnd.C;
        if numel(trnd.C)==0; error([mfilename,'Problem with trnd']); end
      end
      i_trnd=i_trnd+1;
      trnd.C(i_trnd,:)=[n m x(3)];
      trnd.S(i_trnd,:)=[n m x(4)];
    elseif strcmp(s(1:4),'acos')
      if isempty(acos.C)
        acos.C=zeros(grep_nr_occurences(filename,'acos'),4); acos.S=acos.C;
        if numel(asin.C)==0; error([mfilename,'Problem with acos']); end
      end
      i_acos=i_acos+1;
      acos.C(i_acos,:)=[n m x(3) x(end)];
      acos.S(i_acos,:)=[n m x(4) x(end)];
    elseif strcmp(s(1:4),'asin')
      if isempty(asin.C)
        asin.C=zeros(grep_nr_occurences(filename,'asin'),4); asin.S=asin.C;
        if numel(asin.C)==0; error([mfilename,'Problem with asin']); end
      end
      i_asin=i_asin+1;
      asin.C(i_asin,:)=[n m x(3) x(end)];
      asin.S(i_asin,:)=[n m x(4) x(end)];
    else
      error([mfilename,'A problem occured in gfc data.']);
    end
    s=fgets(fid);
  end
  fclose(fid);
  % the tide system is not documented in some models
  if strcmp(header.modelname,'GROOPS')
    header.tide_system='zero_tide';
  end
  % handle the permanent tide deformation system
  mi.C(3,1)=gravity.zero_tide(mi.C(3,1),header.tide_system);
  %initializing data object
  m=gravity(...
    time,...
    gravity.dtc('cs','y',mi),...
    'GM',header.earth_gravity_constant,...
    'R',header.radius,...
    'descriptor',header.modelname,...
    'cdate',datetime('now'),...
    'origin',header.filename...
  );
  if any(ei.C(:)~=0) || any(ei.S(:)~=0)
    %initializing error object
    e=gravity(...
      time,...
      gravity.dtc('cs','y',ei),...
      'GM',header.earth_gravity_constant,...
      'R',header.radius,...
      'descriptor',['error of ',header.modelname],...
      'cdate',datetime('now'),...
      'origin',header.filename...
    );
  else
    e=gravity.unit(m.lmax,...
      'scale',0,...
      't',time,...
      'GM',header.earth_gravity_constant,...
      'R',header.radius,...
      'descriptor',['error of ',header.modelname],...
      'cdate',datetime('now'),...
      'origin',header.filename...
    );
  end
end
function [m,e]=load_mod(filename,time)
  %default time
  if ~exist('time','var') || isempty(time)
    time=datetime('now');
  end
  %loading data
  [mi,headerstr]=file.textscan(filename);
  %init constants
  header=struct(...
    'GM',gravity.parameters('GM'),...
    'R',gravity.parameters('R'),...
    'name','unknown'...
  );
  %going through all header lines
  for i=1:numel(headerstr)
    %checking if this is one of the ID strings
    for j=1:fieldnames(header)
      fieldname=j{1};
      if keyword_search(headerstr{i},[fieldname,':'])
        valuestr=strrep(headerstr{i},[fieldname,':'],'');
        switch class(header.(fieldname))
        case 'double'
          header.(fieldname)=str2double(valuestr);
        case 'char'
          header.(fieldname)=valuestr;
        otherwise
          error([mfilename,': cannot handle class ',class(header.(fieldname)),'.'])
        end
      end
    end
  end
  %initializing data object
  m=gravity(...
    time,...
    gravity.dtc('mod','y',mi),...
    'GM',header.GM,...
    'R',header.R,...
    'descriptor',header.name,...
    'cdate',datetime('now'),...
    'origin',filename...
  );
  %no error info on mod format
  e=[];
end
function [m,e]=load_esamtm(filename,time)
  %default time
  if ~exist('time','var') || isempty(time)
    time=datetime('now');
  end
  %loading data
  load(filename,'mod')
  %initializing data object
  m=gravity(...
    time,...
    gravity.dtc('mod','y',mod.mod),...
    'GM',mod.GM,...
    'R',mod.R,...
    'descriptor','ESA Mass Transport Model',...
    'cdate',datetime('now'),...
    'origin',filename...
  );
  %no error info on mod format
  e=[];
end
%% Aux functions
function out=keyword_search(line,keyword)
    out=strncmp(strtrim(line),       keyword,         length(keyword)) || ...
        strncmp(strtrim(line),strrep(keyword,' ','_'),length(keyword));
end
function out=grep_nr_occurences(filename,pattern)
   [~, result] =system(['grep -c ',pattern,' ',filename]);
   out=str2double(result);
end

%% Spherical Harmonic synthesis
function [long,lat,grid_out]=mod_sh_synth(c,s,lat,NLon)
% [LONG,LAT,GRID]=MOD_SH_SYNTH(C,S) is the low-level spherical harmonic
% synthesis routine. It should be a self-contained script. Inputs C and S
% are the cosine and sine coefficients, organized in lower triangle
% matrices, with constant degree in each row and constant order in each
% column, sectorial coefficients are in the diagonals, zonal coefficients
% are in the first column of matrix C (first column of S matrix is zeros).
%
%   MOD_SH_SYNTH(C,S,LAT,NLON) with the optional inputs LAT, NLON causes
%   the synthesis to be done on the latitudes specified with vector LAT and
%   along the #NLON equidistant number of points along the longitudes
%   domain.
%
%   Output LONG is a horizontal vector [0:2*pi].
%   Output LAT is a vertical vector [-pi/2:pi/2].
%   Output GRID is a matrix with size [length(lat) NLon]
%
%   All angular quantities are in radians.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>
% List of changes:
%
%   P.Inacio <p.m.g.inacio@tudelft.nl>, 10/2011, Added check for repeated
%       meridians at 0 and 2*pi. If they both exist and have the same Also
%       added a more robust check that grids are regular taking into
%       account a numerical error threshold.
%   P.Inacio <p.m.g.inacio@tudelft.nl>, 03/2012, The algorightm upon which
%       this function relies only uses the number of points along
%       longitude. Therefore, to avoid misinterpretations the input LONG
%       has been replaces with the number of points along the longitude
%       domain. Further interpolation to arbitrary sets of longitudes,
%       regular or irregular is implemented outside this routine, in
%       mod_synthesis.

  % Defaults
  %check dimensions
  if any(size(c) ~= size(s))
       error([mfilename,': inputs <c> and <s> must have the same size.'])
  end
  if size(c,1) ~= size(c,2)
      error([mfilename,': inputs <c> and <s> must be square matrices'])
  end

  %calculating max resolution
  L = size(c,1)-1;
  N = max([1 L]);
  % %%% OLD CODE
  % N = max([1,min([L,360])]);
  % %%%

  %building latitude domain (if not given)
  if ~exist('lat','var') || isempty(lat)
      lat=linspace(-pi/2,pi/2,N)';
  end
  %setting default number of points along latitude
  if ~exist('NLon','var') || isempty(NLon)
      NLon = 2*N;
  elseif ~isscalar(NLon) || ~isnumeric(NLon)
      error([mfilename,': input <Nlon> must be a numeric scalar.'])
  end

  %create internal longitude domain - only used for output.
  long=linspace(0,2*pi,NLon+1);
  % removing duplicate zero longitude
  long(end)=[];

  %checking (possible) inputs
  % NOTE: Requiring latitude to be vertical is a bit stupid, but it general
  %       it allows it to be distinguished from the longitude avoiding
  %       possible usage error.
  if size(lat,2) ~= 1
      error([mfilename,': input <lat> must be a vertical vector.'])
  end

  %check latitude domain
  if max(lat) > pi/2 || min(lat) < -pi/2
      error([mfilename,': input <lat> does not seem to be in radians or outside legal domain [-pi/2,pi/2].'])
  end

  %need co-latitude
  clat=pi/2-lat;

  % calculating Legendre polynomials, P{latitude}
  P=legendre_latitude(L,clat);

  %building the grids
  grid_out = zeros(length(clat),NLon);
  for i=1:length(clat)
      %calculating Fourier coefficients
      a=zeros(L+1,1);
      b=zeros(L+1,1);
      for m=0:L
          n=m:L;
          a(m+1) = sum(c(n+1,m+1).*P{i}(n+1,m+1));
          b(m+1) = sum(s(n+1,m+1).*P{i}(n+1,m+1));
      end

      %FFT
      grid_out(i,:) = real(fft(a+1i*b,NLon))';
  end
end
function out = legendre_latitude(L,lat)
  %getting legendre coefficient, per degree
  P=legendre_degree(L,lat);
  %building legendre coefficients, per latitude
  out = cell(length(lat),1);
  for i=1:length(lat)
      out{i} = zeros(L+1,L+1);
      for n=0:L
          out{i}(n+1,1:n+1) = P{n+1}(:,i)';
      end
  end
end
function out = legendre_degree(L,lat)
  if ~isvector(lat)
      error([mfilename,': input <lat> must be a vector.'])
  end
  %getting legendre coefficients, per degree
  out=cell(1,L+1);
  for n=0:L
      out{n+1}=legendre(n,cos(lat'),'norm')*2;
      out{n+1}(1,:)=out{n+1}(1,:)/sqrt(2);
  %     % P.Inacio - testing - manual normalization
  %     m = (0:n)';
  %     out{n+1}=legendre(n,cos(lat'));
  %     out{n+1}=repmat(sqrt((2*n+1)*factorial(n-m)./factorial(n+m)),[1 length(lat)]).*out{n+1};
  %     % P.Inacio - testing - no normalization
  %     out{n+1}=legendre(n,cos(lat'),'norm')*2;
  end
end

%% Auxiliarly data
function [t,s,e,d]=GetGRACEC20(varargin)
  %parse arguments that are required later
  v=varargs.wrap('sources',{...
    {...
      'version',     'TN-14', @ischar;...
      'mode',         'read', @ischar;...
      'start',time.zero_date, @isdatetime;...
      'stop',  time.inf_date, @isdatetime;...
      'data_dir', file.orbdir('auxiliary'),  @(i) ischar(i) && exist(i,'dir')~=0;
    },...
  },varargin{:});
  %some default parameters
  t_idx=1;
  s_idx=3;
  e_idx=5;
  CommentStyle='*';
  datfmt='%7.1f%10.4f%22.13f%8.4f%8.4f';
  %upper-case version name 
  v.version=upper(v.version);
  %parse dependent arguments (can be over-written)
  %(NOTICE: upper is redundant but the preprocessor shows non-capitalized cases)
  switch upper(v.version)
    case 'TN-07'
      datfil=[v.version,'_C20_SLR.txt'];
      daturl=['ftp://podaac.jpl.nasa.gov/allData/grace/docs/',datfil]; %NOTICE: this no longer works
    case 'TN-11'
      datfil=[v.version,'_C20_SLR.txt'];
%       daturl=['https://podaac-w10n.jpl.nasa.gov/w10n/allData/grace/docs/',datfil];
      daturl=['https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/',datfil];
      datfmt='%7.1f%10.4f%22.13f%8.4f%8.4f%8.1f%10.4f';
      t_idx=[1 6];
    case 'CSR-RL06'
      datfil='C20_RL06.txt';
      daturl=['http://download.csr.utexas.edu/pub/slr/degree_2/',datfil];
      datfmt='%10.4f%19.10f%8.4f%8.4f%8.4f%16.4f%16.4f';
      CommentStyle='#';
      s_idx=2;
      e_idx=4;
    case 'TN-14'
      datfil=[v.version,'_C30_C20_GSFC_SLR.txt'];
      daturl=['https://podaac-tools.jpl.nasa.gov/drive/files/allData/grace/docs/',datfil];
      datfmt='%7.1f%10.4f%22.13f%8.4f%8.4f%22.13f%8.4f%8.4f%8.1f%10.4f';
    case 'GSFC'
      datfil='GSFC_SLR_C20_GSM_replacement.txt';
      daturl=['https://earth.gsfc.nasa.gov/sites/default/files/neptune/images_db/',datfil];
    case 'GSFC-7DAY'
      datfil='GSFC_SLR_C20_7day.txt';
      daturl='personal communication: email from bryant.d.loomis@nasa.gov';
%       datfmt='%10.4f%23.13f'; %NOTICE: this was used for the format of the data sent by email
      datfmt='%12.6f%21.13f'; %NOTICE: this is used for the format of the data resulting from the online file
      t_idx=1;
      s_idx=2;
      e_idx=0;
    otherwise
      error(['Cannot handle version ''',v.version,'''.'])
  end
  %append v.file and v.url
  v=varargs.wrap('sources',{v,{...
        'file',fullfile(v.data_dir,datfil),  @ischar;...
        'url',daturl,                         @ischar;...
      }},varargin{:});
  switch lower(v.mode)
  case 'data_file'
    t=v.file;
  case 'read' %'set' if not already, then 'get'
    if ~file.exist(v.file)
      GetGRACEC20('mode','set','version',v.version,'data_dir',v.data_dir);
    end
    [t,s,e,d]=GetGRACEC20('mode','get','version',v.version,'file',v.file,'start',v.start,'stop',v.stop);
  case 'reload' %'set' then 'get'
    GetGRACEC20('mode','set','version',v.version,'data_dir',v.data_dir);
    [t,s,e,d]=GetGRACEC20('mode','get','version',v.version,'file',v.file,'start',v.start,'stop',v.stop);
  case 'plot'
    %NOTICE: this is a low-level plot, without much features
    [t,s,e,d]=GetGRACEC20(varargin{:},'mode','read');
    plotting.figure;
    plot(t ,s ,'o-','MarkerSize',4), hold on
    plotting.enforce('plot_title',d,'plot_legend_location','none');
  case 'get' %read the downloaded data
    %open the file
    fid=file.open(v.file);
    %sanity
    if fid<=0
      error([mfilename,': cannot open the data file ''',v.file,'''.'])
    end
    %default descriptor
    d=['C20 time series, version ',v.version];
    %get header info (NOTICE: upper is redundant but the preprocessor shows non-capitalized cases)
    switch upper(v.version)
    case {'GSFC-7DAY','CSR-RL06'}
      found=true;
    otherwise
      found=false;
    end
    c=0;
    while ~found
      c=c+1;
      line=fgetl(fid);
      if line<0
        error(['Could not find keyword PRODUCT in data file ',v.file])
      end
      found=any(str.contains(line,{'PRODUCT:','Product:'}));
      %these strings should come in the same order as here
      if any(str.contains(line,{'TITLE:','Title:'}))
        d=strtrim(str.clean(line,{'TITLE:','Title:'}));
      end
      if str.contains(line,'UPDATE HISTORY:')
        d=strjoin({d,strtrim(strrep(line,'UPDATE HISTORY:',''))},newline);
      end
      if str.contains(line,'Data span:')
        d=strjoin({d,strtrim(strrep(line,'Data span:',''))},newline);
      end
    end
%     str.say('read through ',c,'header lines')
    %read the data
    dat=textscan(fid,datfmt,'CommentStyle',CommentStyle);
    %close the file
    fclose(fid);
    % outputs (NOTICE: upper is redundant but the preprocessor shows non-capitalized cases)
    switch upper(v.version)
    case {'GSFC-7DAY','CSR-RL06'}
      t=years(mean([dat{t_idx}],2))+datetime('0000-01-01 00:00:00');
    otherwise
      t=datetime(mean([dat{t_idx}],2),'ConvertFrom','modifiedjuliandate');
    end
    % enforce start/stop
    idx=t>=v.start & t<=v.stop;
    t=t(idx);
    %inform
    str.say('start/stop of',d,':',t(1),t(end))
    s=dat{s_idx}(idx);
    if e_idx>0
      e=dat{e_idx}(idx)*1e-10;
    else
      e=zeros(size(s));
    end
  case 'set' %download the data
    if url.is(v.url)
      fid=file.open(v.file,'w+');
      t=webread(v.url);
      fprintf(fid,'%s',t);
      fclose(fid);
    else
      t=[];
    end
    s=[];e=[];d=[];
  otherwise
    error([mfilename,': unknown mode ''',v.mode,'''.'])
  end
end
% function out=GRACEC20_periods(version,descriptor)
%   switch version
%   case 'TN-07'
%     error('implementation needed!')
%   case 'TN-11'
%     if str.contains(descriptor,'May 10, 2019')
%       out=[...
%         3.652426000000e+02;...
%         1.826213000000e+02;...
%         3.840000000000e+03;...
%         1.920000000000e+03;...
%         8.533333333333e+02;...
%         1.097142857143e+03;...
%         5.485714285714e+02;...
%         1.536000000000e+03;...
%         2.648275862069e+02;...
%         2.327272727273e+02;...
%       ];
%       return
%     else
%       %this was relevant to the version previous to the one above
%       out=[365.2426 182.6213 7936 992 566.85714 1984 91.218391 881.77778 721.45455 ...
%         214.48649 1587.2 529.06667 377.90476 122.09231 661.33333 273.65517 330.66667];
%     end
%   otherwise
%     error(['Cannot handle version ''',version,'''.'])
%   end
% end