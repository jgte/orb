classdef gravity < simpletimeseries
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
                  200    -0.007],...
        'units',  'non-dim'...
    );
%   Supported units are the following,
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
    valid_units={...
      'non-dim',...
      'eqwh',...
      'geoid',...
      'potential',...
      'anomalies',...
      'vert-grav-grad',...
      'gravity'...
    };
    parameter_list=struct(...
      'GM',    struct('default',gravity.default_list.GM,    'validation',@(i) isdouble(i)),...
      'R',     struct('default',gravity.default_list.R,     'validation',@(i) isdouble(i)),...
      'units', struct('default',gravity.default_list.units, 'validation',@(i) ischar(i))...
    );

  end
  %read only
  properties(SetAccess=private)
    GM
    R
    units
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    
  end
  %calculated only when asked for
  properties(Dependent)
    lmax
    mat
    cs
    tri
    mod
  end
  methods(Static)
    function out=parameters
      out=fieldnames(gravity.parameter_list);
    end
    function out=issquare(in)
      in_sqrt=sqrt(in);
      out=(in_sqrt-round(in_sqrt)==0);
    end
    function out=default
      out=gravity.default_list;
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
      lmax=y_lmax(mat(:));
    end
    function [s1,s2]=mat_length(lmax)
      s1=lmax+1;
      s2=s1;
    end
    %% cs representation
    function out=cs_valid(in)
      out=isfield(in,'C') && isfield(in,'S');
      if ~out,return,end
      z=triu(in.C,1)+triu([in.S(:,2:end),zeros(size(in.C,1),1)],0);
      out=all(size(in.C)==size(in.S)) && size(in.C,1) == size(in.C,2) && all(z(:)==0) ;
    end
    function lmax=cs_lmax(cs)
      lmax=gravity.y_lmax(cs.C(:));
    end
    function [s1,s2]=cs_length(lmax)
      [s1,s2]=gravity.mat_length(lmax);
    end
    %% tri representation
    function out=tri_valid(tri)
      lmax=gravity.tri_lmax(tri);
      out=size(tri,1)==lmax+1 && round(lmax)==lmax && gravity.cs_valid(gravity.tri2cs(tri));
    end
    function lmax=tri_lmax(tri)
      lmax=(size(tri,2)+1)/2-1;
    end
    function [s1,s2]=tri_length(lmax)
      s1=mat_length(lmax);
      s2=2*(lmax+1)-1;
    end
    %% mod representation
    function out=mod_valid(mod)
      n=gravity.mod_lmax(mod)+1;
      out=size(mod,2)==4 && size(mod,1)==n*(n+1)/2;
    end
    function lmax=mod_lmax(mod)
      lmax=max(mod(:,1));
    end
    function [s1,s2]=mod_length(lmax)
      s1=(lmax*(lmax+1))/2;
      s2=2*(lmax+1)-1;
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
      %shortcuts
      n=gravity.mod_lmax(mod)+1;
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
    function [s1,s2]=dtlength(type,in)
      switch lower(type)
      case 'y';  s2=0;s1=gravity.y_length(in);
      case 'mat';[s1,s2]=gravity.mat_length(in);
      case 'cs'; [s1,s2]=gravity.cs_length(in);
      case 'tri';[s1,s2]=gravity.tri_length(in);
      case 'mod';[s1,s2]=gravity.mod_length(in);
      otherwise
        error([mfilename,': unknown data type ''',from,'''.'])
      end
    end
    %% constructors
    function obj=unit(lmax,scale)
      if ~exist('scale','var') || isempty(scale)
          scale=ones(lmax+1,1);
      end
      %sanity
      if min(size(scale)) ~= 1 || lmax+1 ~= numel(scale)
        error([mfilename,': input ''scale'' has to be a vector with length lmax+1'])
      end
      %create unitary triangular matrix
      t=gravity.dtc('mat','tri',ones(lmax+1));
      %scale along rows
      t=scale(:)*ones(1,size(t,2)).*t;
      %initialize
      obj=gravity(...
        datetime('now'),...
        gravity.dtc('tri','y',t)...
      );
    end
    function obj=unit_amplitude(lmax)
      obj=gravity.unit(lmax,1./sqrt(2*(0:lmax)+1));
    end
    function obj=unit_rms(lmax)
      obj=gravity.unit(lmax,gravity.unit(lmax).drms);
    end
    %general test for the current object
    function out=test_parameters(field,l,w)
      switch field
      case 'something'
        %To be implemented
      otherwise
        out=simpledata.test_parameters(field,l,w);
      end
    end
    function test(l)
      
      if ~exist('l','var') || isempty(l)
        l=4;
      end

      nr_coeff=(l+1)^2;
      
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
      
      disp('--- unit amplitude')
      a=gravity.unit_amplitude(l);
      disp('- C')
      disp(a.cs(1).C)
      disp('- S')
      disp(a.cs(1).S)
      disp('- tri')
      disp(a.tri{1})
      disp('- mod')
      disp(a.mod{1})
      disp('- das')
      disp(a.das)
      
      disp('--- unit rms')
      a=gravity.unit_rms(l);
      disp('- tri')
      disp(a.tri{1})
      disp('- drms')
      disp(a.drms)

      disp('--- change R')
      b=a.change_R(a.R*2);
      disp('- tri')
      disp(b.tri{1})
      
    end
  end
  methods
    %% constructor
    function obj=gravity(t,y,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't'); %this can be char, double or datetime
      p.addRequired( 'y',     @(i) simpledata.valid_y(i));
      %declare parameters
      for j=1:numel(gravity.parameters)
        %shorter names
        pn=gravity.parameters{j};
        %declare parameters
        p.addParameter(pn,gravity.parameter_list.(pn).default,gravity.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(t,y,varargin{:});
      % call superclass
      obj=obj@simpletimeseries(p.Results.t,p.Results.y,varargin{:});
      % save parameters
      for i=1:numel(gravity.parameters)
        %shorter names
        pn=gravity.parameters{i};
%         %parameter 'units' has already been handled when calling simpledata
%         %constructor, so skip it
%         if strcmp(pn,'units')
%           continue
%         end
        if ~isscalar(p.Results.(pn))
          %vectors are always lines (easier to handle strings)
          obj.(pn)=transpose(p.Results.(pn)(:));
        else
          obj.(pn)=p.Results.(pn);
        end
      end
    end
    function obj=assign(obj,y,varargin)
      %pass it upstream
      obj=assign@simpletimeseries(obj,y,varargin{:});
    end
    function obj=copy_metadata(obj,obj_in)
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in);
      %propagate parameters of this object
      parameters=gravity.parameters;
      for i=1:numel(parameters)
        if isprop(obj,parameters{i}) && isprop(obj_in,parameters{i})
          obj.(parameters{i})=obj_in.(parameters{i});
        end
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
          mat_new{i}=zeros(l,l);
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
    function out=get.lmax(obj)
      out=sqrt(obj.width)-1;
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
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('mat',in)));
      %build data 
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(in{i});
      end
      %assign data
      obj=obj.assign(y_now);
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
      if ~isstruct(in)
        error([mfilename,': input <in> must be a cell array of matrices.'])
      end
      if (numel(in)~=obj.length)
        error([mfilename,': cannot handle input <in> if it does not have the same number of elements as obj.length.'])
      end
      %make room for data
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('tri',in)));
      %build data 
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(gravity.cs2mat(gravity.tri2cs(in)));
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
      if ~isstruct(in)
        error([mfilename,': input <in> must be a cell array of matrices.'])
      end
      if (numel(in)~=obj.length)
        error([mfilename,': cannot handle input <in> if it does not have the same number of elements as obj.length.'])
      end
      %make room for data
      y_now=zeros(obj.length,gravity.dtlength('y',gravity.dtlmax('mod',in)));
      %build data 
      for i=1:obj.length
        y_now(i,:)=gravity.mat2y(gravity.cs2mat(gravity.mod2cs(in)));
      end
      %assign data
      obj=obj.assign(y_now);
    end
    %% radius
    function obj=change_R(obj,r)
      deg=(0:obj.lmax)'*ones(1,obj.lmax+1);
      scale=(obj.R/r).^(deg+1);
      cs_now=obj.cs;
      for i=1:obj.length
        cs_now(i).C=cs_now(i).C.*scale;
        cs_now(i).S=cs_now(i).S.*scale;
      end
      obj.cs=cs_now;
    end
    %% derived quantities
    function out=dorders(obj)
      % number of orders in each degree
      out=2*(1:obj.lmax+1)-1;
    end
    %degree RMS
    function out=drms(obj)
      das=obj.das;
      out=zeros(size(das));
      l=obj.dorders;
      for i=1:obj.length
        out(i,:)=das(i,:)./sqrt(l);
      end
    end
    %cumulative degree RMS
    function out=cumdrms(obj)
      out=sqrt(cumsum(obj.drms.^2,2));
    end
    % returns degree amplitude spectrum for each row of obj.y. The output
    % matrix has in each row the epochs of obj.y (corresponding to the epochs
    % of the models) and for each column i the degree i-1 (this is implicit).
    function out=das(obj)
      out=zeros(obj.length,obj.lmax+1);
      tri_now=obj.tri;
      for i=1:obj.length
        %compute DAS
        out(i,:) = sqrt(sum(tri_now{i}.^2,2));
      end
    end
    % returns the cumulative degree amplitude spectrum for each row of obj.y.
    % the output matrix is arranged as <das>.
    function out=cumdas(obj)
      out=cumsum(obj.das,2);
    end
    
  end
end
