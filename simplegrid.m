classdef simplegrid < simpletimeseries
  %static
  properties(Constant)
    %default value of some internal parameters
    parameter_list={...
      'sp_tol',   1e-8,  @num.isscalar;...
      'sp_units', 'deg',@ischar;...  %TODO: need to convert non-deg units at input/output (internally deg are expected)
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'sp_units'};
    %catchment list: name, lat, lon
    catchment_list={...
      'Alaska',        [  56   65],[-151 -129];...
      'Columbia',      [  38   50],[-125 -110];...
      'Mississippi',   [  29   44],[-101  -80];...
      'Amazon',        [ -17,   3],[ -76  -47];...
      'Orinoco'        [   3,  12],[ -72  -59];...
      'La Plata',      [ -34, -19],[ -65  -50];...
      'Greenland',     [  60   85],[ -60  -37];...
      'W Sub-Sahara',  [   5   15],[ -15   -1];...
      'E Sub-Sahara',  [   1   13],[   8   35];...
      'Congo Zambezi', [ -23   -3],[  14   38];...
      'W Antarctica',  [ -80  -70],[-140  -85];...
      'E Antarctica',  [ -80  -68],[  80  130];...
      'N Australia',   [ -24  -10],[ 124  145];...
      'Ganges-Brahm',  [  15   30],[  72   89];...
      'Indochina',     [  12   29],[  93  105];...
      'Siberia',       [  57   72],[  68  109];...
      'Volga',         [  53   61],[  34   56];...
      'Danube',        [  43   48],[  13   28];...
};
%       'Great Lakes',   [  50  67],[-100  -70];...
%       'Severny Island',[  66  78],[  38   67];...

  end
  %read only
  properties(SetAccess=private)
    loni
    lati
  end
  properties(GetAccess=public,SetAccess=public)
    sp_tol
    sp_units
  end
  %calculated only when asked for
  properties(Dependent)
    lat
    lon,lon180,lon360
    map
    vecmat
    list
    flatlist
    matrix
    latSpacing
    lonSpacing
    latLimit
    lonLimit
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(simplegrid.parameter_list); end
      out=v.picker(varargin{:});
    end
    %% spacing
    function out=spacing(in)
      out= min(diff(unique(in(:),'sorted')));
      if isempty(out); out=1; end
    end
    function out=getSpacing(lon,lat)
      out = [ simplegrid.spacing(lon) simplegrid.spacing(lat) ];
    end
    function out=spacing_size(lon,lat,radius)
      if ~exist('radius','var') || isempty(radius)
        radius=gravity.parameters('Rm');
      end
      out=deg2rad(simplegrid.getSpacing(lon,lat))*radius;
    end
    function out=sh_lmax(lon,lat)
      %NOTICE: this function is the same as:
      % floor(gravity.hwl2deg(simplegrid.spacing_size(lon,lat)))
      out = min([numel(lat)-1,floor((numel(lon)-1)/2)]);
    end
    function out=lat_spacing_valid(spacing)
      out= isnumeric(spacing) && isscalar(spacing) && mod(180,spacing)==0;
    end
    function out=lon_spacing_valid(spacing)
      out= isnumeric(spacing) && isscalar(spacing) && mod(360,spacing)==0;
    end
    %% lon and lat handlers
    function out=getLimits(lon,lat)
      lon=simplegrid.setlon360(lon);
      out = [ min(lon(:)) max(lon(:)) min(lat(:)) max(lat(:)) ];
      %some sanity
      assert(~(any(lon <   0) || any(lon > 360)),'Ilegal longitude domain. Debug needed!')
      assert(~(any(lat < -90) || any(lat >  90)),'Ilegal latitude domain. Debug needed!')
    end
    function out=delta(l,u,n)
      out=(u-l)/(n-1);
    end
    function out=lon_default(n)
      %NOTICE: this used to be simply linspace(0,360,n) but that may introduce NaN at 360  
      %NOTICE: the use of delta may look redundant but linspace is not very accurate and introduces small variations
      out=0:simplegrid.delta(0,360,n):360;
      %NOTICE: there is particular requirement for having (or not) longitude at degree 360 (degree 0 should be there, if relevant);
      %        wrap the code under obj.{add|del}360lon for those cases when either is needed
%       out=out(1:end-1);
    end
    function out=lat_default(n)
      %NOTICE: this looks redundant but linspace is not very accurate and introduces small variations
      out=transpose(-90:simplegrid.delta(-90,90,n):90);
    end
    function out=lat_stepped(delta_deg)
      %NOTICE: this is the same as: out=simplegrid.lat_default((180+delta_deg)/delta_deg);
      out=-90:delta_deg:90;
      assert(out(end)==90,'Input delta_deg is ilegal because it does not produce a complete latitude domain')
    end
    function out=lon_stepped(delta_deg)
      %NOTICE: this is the same as: out=simplegrid.lon_default((360+delta_deg)/delta_deg);
      out=0:delta_deg:360;
      assert(out(end)==360,'Input delta_deg is ilegal because it does not produce a complete latitude domain')
    end
    function out=lat_stepped_length(delta_deg)
      if isfinite(delta_deg)
        if delta_deg==0
          out=inf;
        else
          out=(180+delta_deg)/delta_deg;
        end
      else
        out=0;
      end
    end
    function out=lon_stepped_length(delta_deg)
      if isfinite(delta_deg)
        if delta_deg==0
          out=inf;
        else
          out=(360+delta_deg)/delta_deg;
        end
      else
        out=0;
      end
    end
    function out=setlon360(lon)
      out=wrapTo360(lon);
    end
    function out=islon360complete(lon)
      lon360=simplegrid.setlon360(lon);
      sp_tol=simplegrid.parameters('sp_tol')^2;
      out = ...
        (lon360(end)-360                          )^2<sp_tol && ... %check if lon 360 is already there
        (lon360(1)                                )^2<sp_tol && ... %check if current spatial domain includes longitude 0
        (lon360(end-1)+simplegrid.spacing(lon)-360)^2<sp_tol;       %check if current spatial domain should have longitude 360
    end
    function [idx_a,idx_b,out]=common(a,b)
      %a(idx_a)=b(idx_b)=out and is the common entries between a and b
      [ag,bg]=meshgrid(a,b);
      [idx_b,idx_a]=find((ag-bg).^2<simplegrid.parameters('sp_tol')^2);
      if nargout>2; out=a(idx_a); end
    end
    %% vector/matrix representation (lat, lon and t are vectors, map is a 3D matrix with time changing along the right-most dim)
    function     out=vecmat_valid(vecmat)
      try
        simplegrid.vecmat_check(vecmat);
        out=true;
      catch
        out=false;
      end
    end
    function     out=vecmat_check(vecmat)
      
      %easier names
      type='vectmat'; 
      %need structure
      if ~isstruct( vecmat);       error([mfilename,': invalid ',type,': not a structure.']); end
      %need to check lon/lat fields existence now, before determining if they are uniform
      if ~isfield(    vecmat,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end 
      if ~isfield(    vecmat,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(vecmat.lon) && ~isempty(vecmat.lat)
        % Calculate the spacing and limits
        GridLimits  = simplegrid.getLimits( vecmat.lon(:),vecmat.lat(:));
        GridSpacing = simplegrid.getSpacing(vecmat.lon(:),vecmat.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon',          GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',transpose(GridLimits(3):GridSpacing(2):GridLimits(4))...
        );
        % sanity
        if numel(out.lon) ~= numel(vecmat.lon) || any( (out.lon(:)-vecmat.lon(:)).^2>simplegrid.parameters('sp_tol')^2 )
          error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.'])
        end
        if numel(out.lat) ~= numel(vecmat.lat) || any( (out.lat(:)-vecmat.lat(:)).^2>simplegrid.parameters('sp_tol')^2 )
          error([mfilename,': invalid ',type,': field ''lat'' is not an uniform domain.'])
        end
      else
        out=struct('lon',[],'lat',[]);
      end
      %map
      if ~isfield(    vecmat,'map'); error([mfilename,': invalid ',type,': field ''map'' is missing.']);end
      if ~isempty(    vecmat.map) && ...
         ~isnumeric(  vecmat.map);   error([mfilename,': invalid ',type,': field ''map'' is not numeric.']);end
      %lon
      if ~isempty(    vecmat.lon)
        if ~isrow(    vecmat.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not a row vector.']);end
        if ~isnumeric(vecmat.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not numeric.']);end
        if  ~isempty( vecmat.map) && ...
            numel(    vecmat.lon) ~= ...
            size(     vecmat.map,2); error([mfilename,': invalid ',type,': field ''lon'' has different nr of elements than rows of ''map''.']);end %#ok<ALIGN>
      end
      %lat
      if ~isempty(    vecmat.lat)
        if ~iscolumn( vecmat.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not a column vector.']);end
        if ~isnumeric(vecmat.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not numeric.']);end
        if ~isempty(  vecmat.map) && ...
            numel(    vecmat.lat) ~= ...
            size(     vecmat.map,1); error([mfilename,': invalid ',type,': field ''lat'' has different nr of elements than columns of ''map''.']);end %#ok<ALIGN>
      end
      %t
      if ~isfield(   vecmat,'t');   error([mfilename,': invalid ',type,': field ''t'' is missing.']);end
      if ~isvector(  vecmat.t  );   error([mfilename,': invalid ',type,': field ''t'' is not a vector.']);end
      if ~isdatetime(vecmat.t  ) && ...
         ~isnumeric( vecmat.t  );   error([mfilename,': invalid ',type,': field ''t'' is not datetime or numeric.']);end
      if ~isempty(   vecmat.map) && ...
          numel(     vecmat.t  )~=...
          size(      vecmat.map,3); error([mfilename,': invalid ',type,': field ''t'' has different nr of elements than pages of ''map''.']);end %#ok<ALIGN>
    end
    function     out=vecmat_size(vecmat)
      out=size(vecmat.map);
    end
    function vectmat=vecmat_init(t,map,lon,lat)
      vectmat=struct('t',t,'lon',simplegrid.setlon360(lon),'lat',lat,'map',map);
    end
    %% list representation (lat, lon and t are vectors, map is a 2D matrix, with time changing along the rows)
    function       out=list_valid(list)
      try
        simplegrid.list_check(list);
        out=true;
      catch
        out=false;
      end
    end
    function [out,idx]=list_check(list)
      %easier names
      type='list';
      %need structure
      if ~isstruct(  list);      error([mfilename,': invalid ',type,': not a structure.']); end
      %need to check lon/lat fields existence now, before determining if they are uniform
      if ~isfield(   list,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end
      if ~isfield(   list,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(list.lon) && ~isempty(list.lat)
        % Calculate the spacing and limits
        GridLimits  = simplegrid.getLimits( list.lon(:),list.lat(:));
        GridSpacing = simplegrid.getSpacing(list.lon(:),list.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon', GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',(GridLimits(3):GridSpacing(2):GridLimits(4))'...
        );
        % get indexes
        %NOTICE: this is slow
        idx.lon=nan(1,numel(list.lon));
        for i=1:numel(out.lon)
          idx.lon( (out.lon(i)-list.lon).^2<simplegrid.parameters('sp_tol')^2 )=i;
        end
        idx.lat=nan(1,numel(list.lat));
        for i=1:numel(out.lat)
          idx.lat( (out.lat(i)-list.lat).^2<simplegrid.parameters('sp_tol')^2 )=i;
        end
        if any(isnan(idx.lat)); error([mfilename,': invalid ',type,': field ''lat'' is not an uniform domain.']);end
        if any(isnan(idx.lon)); error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.']);end
      else
        idx=struct('lon',[],'lat',[]);
        out=struct('lon',[],'lat',[]);
      end  
      %map
      if ~isfield(    list,'map'); error([mfilename,': invalid ',type,': field ''map'' is missing.']);end
      if     ~isempty(list.map)
        if ~isnumeric(list.map); error([mfilename,': invalid ',type,': field ''map'' is not numeric.']);end
        if ~ismatrix( list.map); error([mfilename,': invalid ',type,': field ''map'' is not a matrix.']);end
      end
      %lon
      if ~isempty(    list.lon)
        if ~isrow(    list.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not a row vector.']);end
        if ~isnumeric(list.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not numeric.']);end
        if ~isempty(  list.map) && ...
            numel(    list.lon) ~= ...
            size(     list.map,2); error([mfilename,': invalid ',type,': field ''lon'' has different elements than columns of ''map''.']);end %#ok<ALIGN>
      end
      %lat
      if ~isempty(    list.lat)
        if ~isrow(    list.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not a row vector.']);end
        if ~isnumeric(list.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not numeric.']);end
        if ~isempty(  list.map) && ...
            numel(    list.lat) ~= ...
            size(     list.map,2); error([mfilename,': invalid ',type,': field ''lat'' has different sizes than columns of ''map''.']);end %#ok<ALIGN>
      end
      %t
      if ~isfield(   list,'t');   error([mfilename,': invalid ',type,': field ''t'' is missing.']);end
      if ~isvector(  list.t);     error([mfilename,': invalid ',type,': field ''t'' is not a vector.']);end
      if ~isdatetime(list.t  ) && ...
         ~isnumeric( list.t  );   error([mfilename,': invalid ',type,': field ''t'' is not datetime or numeric.']);end
      if ~isempty(   list.map) && ...
          numel(     list.t)   ~= ...
          size(      list.map,1); error([mfilename,': invalid ',type,': field ''t'' has different nr of elements than rows of ''map''.']);end %#ok<ALIGN>
    end
    function       out=list_size(list)
      out = [ numel(unique(list.lon(:))) numel(unique(list.lat(:))) numel(list.t) ];
    end
    function      list=list_init(t,map,lon,lat)
      list=struct('t',t,'lon',simplegrid.setlon360(lon),'lat',lat,'map',map);
    end
    %% flatlist representation (lat, lon, t and map are vectors, all have the same length)
    function       out=flatlist_valid(flatlist)
      try
        simplegrid.flatlist_check(flatlist);
        out=true;
      catch
        out=false;
      end
    end
    function [out,idx]=flatlist_check(flatlist) 
      %easier names
      type='flatlist';
      %need structure
      if ~isstruct(  flatlist);      error([mfilename,': invalid ',type,': not a structure.']); end
      %need to check lon/lat fields existence now, before determining if they are uniform
      if ~isfield(   flatlist,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end
      if ~isfield(   flatlist,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(flatlist.lon) && ~isempty(flatlist.lat)
        % Calculate the spacing and limits
        GridLimits  = simplegrid.getLimits( flatlist.lon(:),flatlist.lat(:));
        GridSpacing = simplegrid.getSpacing(flatlist.lon(:),flatlist.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon', GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',(GridLimits(3):GridSpacing(2):GridLimits(4))',...
            't', unique(flatlist.t,'sorted')...
        );
        % get indexes
        %NOTICE: this is slow
        idx.lon=nan(1,numel(flatlist.lon));
        for i=1:numel(out.lon)
          idx.lon( (out.lon(i)-flatlist.lon).^2<simplegrid.parameters('sp_tol')^2 )=i;
        end
        if any(isnan(idx.lon)); error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.']);end
        idx.lat=nan(1,numel(flatlist.lat));
        for i=1:numel(out.lat)
          idx.lat( (out.lat(i)-flatlist.lat).^2<simplegrid.parameters('sp_tol')^2 )=i;
        end
        if any(isnan(idx.lat)); error([mfilename,': invalid ',type,': field ''lat'' is not an uniform domain.']);end
        idx.t=nan(1,numel(flatlist.t));
        fltdn=datenum(flatlist.t);
        for i=1:numel(out.t)
          idx.t( (datenum(out.t(i))-fltdn).^2<simplegrid.parameters('sp_tol')^2 )=i;
        end
        if any(isnan(idx.t)); error([mfilename,': invalid ',type,': field ''t'' contains NaNs, debug needed.']);end
      else
        idx=struct('lon',[],'lat',[],'t',[]);
        out=struct('lon',[],'lat',[],'t',[]);
      end  
      %map
      if ~isfield(    flatlist,'map'); error([mfilename,': invalid ',type,': field ''map'' is missing.']);end
      if     ~isempty(flatlist.map)
        if ~isnumeric(flatlist.map); error([mfilename,': invalid ',type,': field ''map'' is not numeric.']);end
        if ~iscolumn( flatlist.map); error([mfilename,': invalid ',type,': field ''map'' is not a column vector.']);end
      end
      %lon
      if ~isempty(    flatlist.lon)
        if ~iscolumn( flatlist.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not a column vector.']);end
        if ~isnumeric(flatlist.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not numeric.']);end
        if ~isempty(  flatlist.map) && ...
            numel(    flatlist.lon) ~= ...
            numel(    flatlist.map);   error([mfilename,': invalid ',type,': field ''lon'' has different length than ''map''.']);end %#ok<ALIGN>
      end
      %lat
      if ~isempty(    flatlist.lat)
        if ~iscolumn( flatlist.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not a iscolumn vector.']);end
        if ~isnumeric(flatlist.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not numeric.']);end
        if ~isempty(  flatlist.map) && ...
            numel(    flatlist.lat) ~= ...
            numel(    flatlist.map);   error([mfilename,': invalid ',type,': field ''lat'' has different length than ''map''.']);end %#ok<ALIGN>
      end
      %t
      if ~isfield(   flatlist,'t');   error([mfilename,': invalid ',type,': field ''t'' is missing.']);end
      if ~isvector(  flatlist.t);     error([mfilename,': invalid ',type,': field ''t'' is not a vector.']);end
      if ~isdatetime(flatlist.t  ) && ...
         ~isnumeric( flatlist.t  );   error([mfilename,': invalid ',type,': field ''t'' is not datetime or numeric.']);end
      if ~isempty(   flatlist.map) && ...
          numel(     flatlist.t)   ~= ...
          numel(     flatlist.map);   error([mfilename,': invalid ',type,': field ''t'' has different length than ''map''.']);end %#ok<ALIGN>
    end
    function       out=flatlist_size(list)
      out = [ numel(unique(list.lon(:))) numel(unique(list.lat(:))) numel(unique(list.t)) ];
    end
    function      list=flatlist_init(t,map,lon,lat)
      list=struct('t',t,'lon',simplegrid.setlon360(lon),'lat',lat,'map',map);
    end
    %% matrix representation (lat and lon are 2D matrices, map is a 3D matrix with time changing along the right-most dim)
    function    out=matrix_valid(matrix)
      try
        simplegrid.matrix_check(matrix);
        out=true;
      catch
        out=false;
      end
    end
    function    out=matrix_check(matrix)
      %easier names
      type='matrix';
      %need structure
      if ~isstruct( matrix);       error([mfilename,': invalid ',type,': not a structure.']); end
      %need to check lon/lat fields existence now, before determining if they are uniform
      if ~isfield(    matrix,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end 
      if ~isfield(    matrix,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(matrix.lon) && ~isempty(matrix.lat)
        % Calculate the spacing and limits
        GridLimits  = simplegrid.getLimits( matrix.lon(:),matrix.lat(:));
        GridSpacing = simplegrid.getSpacing(matrix.lon(:),matrix.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon',          GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',transpose(GridLimits(3):GridSpacing(2):GridLimits(4))...
        );
        %build meshed lon and lat
        [lon_m,lat_m]=meshgrid(out.lon,out.lat);
        % sanity
        if numel(lon_m) ~= numel(matrix.lon) || any( (lon_m(:)-matrix.lon(:)).^2>simplegrid.parameters('sp_tol')^2 )
          error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.'])
        end
        if numel(lat_m) ~= numel(matrix.lat) || any( (lat_m(:)-matrix.lat(:)).^2>simplegrid.parameters('sp_tol')^2 )
          error([mfilename,': invalid ',type,': field ''lat'' is not an uniform domain.'])
        end
      else
        out=struct('lon',[],'lat',[]);
      end
      %map
      if ~isfield(   matrix,'map'); error([mfilename,': invalid ',type,': field ''map'' is missing.']);end
      if ~isempty(   matrix.map) && ...
          ~isnumeric(matrix.map);   error([mfilename,': invalid ',type,': field ''map'' is not numeric.']);end
      %lon
      if ~isfield(    matrix,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end 
      if ~isempty(    matrix.lon)
        if ~ismatrix( matrix.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not a matrix.']);end
        if ~isnumeric(matrix.lon);   error([mfilename,': invalid ',type,': field ''lon'' is not numeric.']);end
        if  ~isempty( matrix.map) && ...
            any(size( matrix.lon) ~= ...
            size(     matrix.map(:,:,1))) %#ok<ALIGN>
                                     error([mfilename,': invalid ',type,': field ''lon'' has different sizes than rows and columns of ''map''.']);end
      end
      %lat
      if ~isfield(    matrix,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(    matrix.lat)
        if ~ismatrix( matrix.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not a matrix.']);end
        if ~isnumeric(matrix.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not numeric.']);end
        if ~isempty(  matrix.map) && ... 
            any(size( matrix.lat) ~= ...
            size(     matrix.map(:,:,1))) %#ok<ALIGN>
                                     error([mfilename,': invalid ',type,': field ''lat'' has different sizes than rows and columns of ''map''.']);end
      end
      %t
      if ~isfield(   matrix,'t'  ); error([mfilename,': invalid ',type,': field ''t'' is missing.']);end
      if ~isvector(  matrix.t  );   error([mfilename,': invalid ',type,': field ''t'' is not a vector.']);end
      if ~isdatetime(matrix.t  ) && ...
         ~isnumeric( matrix.t  );   error([mfilename,': invalid ',type,': field ''t'' is not datetime or numeric.']);end
      if  ~isempty(  matrix.map) && ... 
          numel(     matrix.t  ) ~= ...
          size(      matrix.map,3) %#ok<ALIGN>
                                    error([mfilename,': invalid ',type,': field ''t'' has different nr of elements than pages of ''map''.']);end
    end
    function    out=matrix_size(matrix)
      out=size(matrix.map);
    end
    function matrix=matrix_init(t,map,lon,lat)
      matrix=struct('t',t,'lon',simplegrid.setlon360(lon),'lat',lat,'map',map);
    end
    %% vecmat and list convertions
    function [vecmat,idx]=list2vecmat(list)
      % inputs are lists of points in rows, changing with times along the columns
      [vecmat,idx]=simplegrid.list_check(list);
      % save indexes
      idx.map = sub2ind([numel(vecmat.lat),numel(vecmat.lon)],idx.lat,idx.lon);
      % save t
      vecmat.t=list.t;
      % assign map
      vecmat.map =  NaN(numel(vecmat.lat),numel(vecmat.lon),numel(list.t));
      if ~isempty(list.map)
        for i=1:numel(vecmat.t)
          tmp=NaN(numel(vecmat.lat),numel(vecmat.lon));
          tmp(idx.map)=list.map(i,:);
          vecmat.map(:,:,i) = tmp;
        end
      end
    end
    function list=vecmat2list(vecmat)
      %check input and get uniform domain
      e=simplegrid.vecmat_check(vecmat);
      % build meshed lon and lat
      [lon_m,lat_m]=meshgrid(e.lon,e.lat);
      % outputs
      list=struct(...
        'lon',transpose(lon_m(:)),...
        'lat',transpose(lat_m(:)),...
        't'  ,vecmat.t,...
        'map',NaN(numel(vecmat.t),numel(lon_m))...
      );
      %skip converting map if input is empty
      if ~isempty(vecmat.map)
        for i=1:numel(list.t)
          list.map(i,:)=reshape(vecmat.map(:,:,i),numel(list.lon),1);
        end
      end
    end
    %% list and flatlist convertions
    function flatlist=vecmat2flatlist(vecmat)
      % need list (vecmat is checked inside vecmat2list)
      list=simplegrid.vecmat2list(vecmat);
      % build 3d meshgrids     
      [lat,t,lon]=meshgrid(vecmat.lat(:),datenum(list.t(:)),vecmat.lon(:));
      %assign outputs: collapse everything into vectors
      flatlist=struct(...
        'lon',lon(:),...
        'lat',lat(:),...
        't',datetime(t(:),'ConvertFrom','datenum'),...
        'map',list.map(:)...
      );
    end
    function vecmat=flatlist2vecmat(flatlist)
      % inputs are lists of points in rows, including lat, lon and time dependencies
      [vecmat,idx]=simplegrid.flatlist_check(flatlist);
      % save indexes
      idx.map = sub2ind([numel(vecmat.lat),numel(flatlist.t),numel(vecmat.lon)],idx.lat,idx.t,idx.lon);
      % assign map
      vecmat.map =  NaN(numel(vecmat.lat),numel(vecmat.lon),numel(vecmat.t));
      if ~isempty(flatlist.map)
        for i=1:numel(flatlist.map)
          vecmat.map(idx.lat(i),idx.lon(i),idx.t(i))=flatlist.map(i);
        end
      end
    end
    %% vecmat and matrix convertions
    function vecmat=matrix2vecmat(matrix)
      % Check if grid is regular
      if all(diff(matrix.lon(1,:),1,2) == 0)
          % Transpose everything
          matrix.lon = transpose(matrix.lon);
          matrix.lat = transpose(matrix.lat);
          matrix.map = permute(  matrix.map,[2,1,3]);
      end
      % inputs are matrices of coordinates and values
      vecmat=simplegrid.matrix_check(matrix);
      % copy the rest
      vecmat.map = matrix.map;
      vecmat.t   = matrix.t;
    end
    function matrix=vecmat2matrix(vecmat)
      e=simplegrid.vecmat_check(vecmat);
      % build meshed lon and lat
      [lon_m,lat_m]=meshgrid(e.lon,e.lat);
      % outputs
      matrix=struct(...
        'lon',lon_m,...
        'lat',lat_m,...
        't'  ,vecmat.t,...
        'map',vecmat.map...
      );
    end
    %% agregator routines
    %data type converter
    function out=dtc(from,to,in)
      %trivial call
      if strcmpi(from,to)
        out=in;
        return
      end
      %convert to required types
      switch lower(from)
        case 'vecmat'
          switch lower(to)
            case 'list';     out=simplegrid.vecmat2list(in);
            case 'flatlist'; out=simplegrid.vecmat2flatlist(in);
            case 'matrix';   out=simplegrid.vecmat2matrix(in);
          end
        case 'list'
          switch lower(to)
            case 'vecmat';   out=                           simplegrid.list2vecmat(in);
            case 'flatlist'; out=simplegrid.vecmat2flatlist(simplegrid.list2vecmat(in));
            case 'matrix';   out=simplegrid.vecmat2matrix(  simplegrid.list2vecmat(in));
          end
        case 'flatlist'
          switch lower(to)
            case 'vecmat'; out=                         simplegrid.flatlist2vecmat(in);
            case 'list';   out=simplegrid.vecmat2list(  simplegrid.flatlist2vecmat(in));
            case 'matrix'; out=simplegrid.vecmat2matrix(simplegrid.flatlist2vecmat(in));
          end
        case 'matrix'
          switch lower(to)
            case 'vecmat';   out=                           simplegrid.matrix2vecmat(in);
            case 'flatlist'; out=simplegrid.vecmat2flatlist(simplegrid.matrix2vecmat(in));
            case 'list';     out=simplegrid.vecmat2list(    simplegrid.matrix2vecmat(in));
          end
      end
    end
    %data type validity
    function c=dtv(type,in)
      switch lower(type)
      case 'vecmat';   c=simplegrid.vecmat_valid(in);
      case 'list';     c=simplegrid.list_valid(in);
      case 'flatlist'; c=simplegrid.flatlist_valid(in);
      case 'matrix';   c=simplegrid.matrix_valid(in);
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    %data type check
    function dtcheck(type,in)
      switch lower(type)
      case 'vecmat';   simplegrid.vecmat_check(in);
      case 'list';     simplegrid.list_check(in);
      case 'flatlist'; simplegrid.flatlist_check(in);
      case 'matrix';   simplegrid.matrix_check(in);
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    %data type length
    function l=dtl(type,in)
      switch lower(type)
      case 'vecmat';   l=simplegrid.vecmat_size(in);
      case 'list';     l=simplegrid.list_size(in);
      case 'flatlist'; l=simplegrid.flatlist_size(in);
      case 'matrix';   l=simplegrid.mat_size(in);
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    %data type init
    function s=dti(t,map,lon,lat,type)
      %default type is vecmat
      if ~exist('type','var') || isempty(type)
        type='vecmat';
      end
      %cell arrays of maps
      if iscell(map) && numel(map) == numel(t)
        %collapse all matrices in the cell array into a 3d matrix
        map=reshape(cat(2,map{:}),[size(map{1}),numel(t)]);
      %3d array with maps on each page
      elseif isnumeric(map) && size(map,3) == numel(t)
        %do nothing, already in expected format
      %2d array with maps on each row
      elseif isnumeric(map) && size(map,1) == numel(t)
        %do nothing, this will be handled later
      else
        error([mfilename,': error in sizes of inputs ''t'' and ''map''.'])
      end
      % build structure with inputs
      si=struct(...
        't',t,...
        'lon',simplegrid.setlon360(lon),...
        'lat',lat,...
        'map',map...
      );
      % determine format of inputs
      if simplegrid.list_valid(si)
        fmt='list';
      elseif simplegrid.vecmat_valid(si)
        fmt='vecmat';
      elseif simplegrid.flatlist_valid(si)
        fmt='flatlist';
      elseif simplegrid.matrix_valid(si)
        fmt='matrix';
      else
        for i={'list','vecmat','flatlist','matrix'}
          try
            simplegrid.([i{1},'_check'])(si)
          catch ME
            disp(ME.message)
          end
        end
        error([mfilename,': cannot handle the format of inputs'])
      end
      % convert to requested type
      s=simplegrid.dtc(fmt,type,si);
    end
    %% constructors
    function obj=unit(n_lon,n_lat,varargin)
      p=inputParser;
      p.addRequired( 'n_lon',  @isnumeric);
      p.addRequired( 'n_lat',  @isnumeric);
      p.addParameter('scale',1,@num.isscalar);
      p.addParameter('bias', 0,@num.isscalar);
      p.addParameter('t',datetime('now'),@isdatetime);
      p.parse(n_lon,n_lat,varargin{:});
      %initialize with unitary grid
      if isscalar(p.Results.n_lat) && isscalar(p.Results.n_lon)
        obj=simplegrid(...
          p.Results.t,...
          ones(p.Results.n_lat,p.Results.n_lon,numel(p.Results.t))*p.Results.scale+p.Results.bias,...
          varargin{:}...
        );
      else
        obj=simplegrid(...
          p.Results.t,...
          ones(numel(p.Results.n_lat),numel(p.Results.n_lon),numel(p.Results.t))*p.Results.scale+p.Results.bias,...
          'lat',p.Results.n_lat,'lon',p.Results.n_lon,...
          varargin{:}...
        );
      end
    end
    % Creates a random model with mean 0 and std 1
    function obj=unit_randn(n_lon,n_lat,varargin)
      p=inputParser;
      p.addRequired( 'n_lon',  @num.isscalar);
      p.addRequired( 'n_lat',  @num.isscalar);
      p.addParameter('scale',1,@num.isscalar);
      p.addParameter('bias', 0,@num.isscalar);
      p.addParameter('t',datetime('now'),@isdatetime);
      p.parse(n_lon,n_lat,varargin{:});
      %initialize with random grid
      obj=simplegrid(...
        p.Results.t,...
        randn(p.Results.n_lat,p.Results.n_lon,numel(p.Results.t))*p.Results.scale+p.Results.bias...
      );
    end
    function obj=slanted(n_lon,n_lat,varargin)
      p=inputParser;
      p.addRequired( 'n_lon',  @num.isscalar);
      p.addRequired( 'n_lat',  @num.isscalar);
      p.addParameter('scale',1,@num.isscalar);
      p.addParameter('bias', 0,@num.isscalar);
      p.addParameter('t',datetime('now'),@isdatetime);
      p.parse(n_lon,n_lat,varargin{:});
      %initialize grid with sum of lon and lat in degrees
      [lon_map,lat_map]=meshgrid(simplegrid.lon_default(p.Results.n_lon),simplegrid.lat_default(p.Results.n_lat));
      %initialize with grid
      obj=simplegrid(...
        p.Results.t,...
        repmat(lon_map+lat_map*p.Results.scale+p.Results.bias,1,1,numel(p.Results.t))...
      );
    end
    %NOTICE: inpupts lon and lat can be scalar (domain size) or vectors (domain entries)
    function obj=landmask(lon,lat,varargin)
      p=inputParser;p.KeepUnmatched=true;p.PartialMatching=false; p.CaseSensitive=true;
      p.addRequired( 'lon',  @isnumeric);
      p.addRequired( 'lat',  @isnumeric);
      p.addParameter('t',     datetime('now'),@isdatetime);
      %values higher than this cuttoff are land, lower are ocean; negative values leave the interpolation unchanged
      p.addParameter('cutoff',-1, @num.isscalar); 
      p.parse(lon,lat,varargin{:});
      %load the data
      fmat=fullfile(file.orbdir('auxiliary'),'landmask.mat');
      fdat=fullfile(file.orbdir('auxiliary'),'landmask.dat');
      if exist(fmat,'file')
        load(fmat,'landmask')
      elseif exist(fdat,'file')
        landmask=simplegrid.load(fdat,t);
        save(fmat,'landmask')
      else
        error(['Cannot find landmask data file, expecting ',fmat,' or ',fdat,'.'])
      end
      %resample to requested resolution
      if isscalar(lon) && isscalar(lat)
        obj=landmask.spatial_resample(lon,lat);
      elseif isvector(lon) && isvector(lat)
        obj=landmask.spatial_interp(lon,lat);
      else
        error('Inputs ''lat'' and ''lon'' must both be scalar or vector')
      end
      %enforce cutoff
      if p.Results.cutoff>0
        land_idx=obj.y>=p.Results.cutoff;
        y_now=obj.y;
        y_now(~land_idx)=0;
        y_now( land_idx)=1;
        obj=obj.assign(y_now);
      end
    end
    function obj=oceanmask(lon,lat,varargin)
      p=inputParser;p.KeepUnmatched=true;
      p.addRequired( 'lon',  @isnumeric);
      p.addRequired( 'lat',  @isnumeric);
      p.addParameter('t',     datetime('now'),@isdatetime);
      p.addParameter('cutoff',1.4, @num.isscalar);
      p.parse(lon,lat,varargin{:});
      %load the data
      fdat=fullfile(file.orbdir('auxiliary'),'wahr.global_ocn_kernel.txt');
      oceanmask=gravity.load(fdat,'mod').grid;
      %resample to requested resolution
      if isscalar(lon) && isscalar(lat)
        obj=oceanmask.spatial_resample(lon,lat);
      elseif isvector(lon) && isvector(lat)
        obj=oceanmask.spatial_interp(lon,lat);
      else
        error('Inputs ''lat'' and ''lon'' must both be scalar or vector')
      end
      %enforce cutoff
      if p.Results.cutoff>0
        ocean_idx=obj.y>=p.Results.cutoff;
        y_now=obj.y;
        y_now(~ocean_idx)=0;
        y_now( ocean_idx)=1;
        obj=obj.assign(y_now);
      end
    end
    function out=load(filename,t)
      fid=fopen(filename);
      d=textscan(fid,'%f %f %s');
      fclose(fid);
      out=simplegrid(t,...            %time
        double(str.logical(d{3})),... %map
        'lon',d{1},...                %lon
        'lat',d{2}...                 %lat
      );
    end
    %% map add-ons
    function h=coast(varargin)
      p=inputParser;
      p.addParameter('datafile',fullfile(file.orbdir('auxiliary'),'coast.mat'));
      p.addParameter('line_color','k',@ischar);
      p.addParameter('line_width',1.5,  @num.isscalar);
      p.addParameter('lon',[-180,180],  @(i) ~num.isscalar(i));
      p.addParameter('lat',[ -90, 90],  @(i) ~num.isscalar(i));
      p.parse(varargin{:});
      coast = load(p.Results.datafile);
      keep_idx=find(...
        ( coast.long<=max(p.Results.lon) & coast.long>=min(p.Results.lon) & ...    
          coast.lat <=max(p.Results.lat) & coast.lat >=min(p.Results.lat) ) | ...
        isnan(coast.lat) | isnan(coast.long)  ...
      );
      h=plot(coast.long(keep_idx),coast.lat(keep_idx),'LineWidth',p.Results.line_width,'Color',p.Results.line_color);
    end
    %% catchment handling
    function out=catchment_idx(name)
      out=cells.strfind(lower(simplegrid.catchment_list(:,1)),lower(name));
      assert(~isempty(out),['Cannot handle catchment ''',name,'''.'])
    end
    function out=catchment_details(name,field)
      if iscellstr(name)
        out=cellfun(@(i) simplegrid.catchment_details(i,field),name,'UniformOutput',false)';
        %need to tweak this
        switch field
        case {'lonlat','latlon'}
          tmp=out;
          out=cell(numel(tmp),2);
          for i=1:numel(tmp)
            out{i,1}=tmp{i}{1};
            out{i,2}=tmp{i}{2};
          end
        end    
        return
      end
      idx=simplegrid.catchment_idx(name);
      switch field
      case 'lat';                out=simplegrid.catchment_list{idx,2};
      case 'lon';                out=simplegrid.catchment_list{idx,3};
      case {'lonlat','latlon'};  out=simplegrid.catchment_list(idx,2:3);
      otherwise; error(['Cannot handle field ''',field,'''.'])
      end    
    end
    function out=catchment_subset(names)
      %trivial call
      if iscell(names) && size(names,2)==3 && iscellstr(names(1,:))
        out=names;
        return
      end
      assert(iscellstr(names),['Input ''names'' must be a cell string, not a ',class(names),'.'])
      out=cell(numel(names),3);
      for i=1:numel(names)
        out(i,:)=simplegrid.catchment_list(simplegrid.catchment_idx(names{i}),:);
      end
    end
    %% smooting/conv aux functions
    function map=map_mirroredges(map,n)
      assert(isscalar(n),'Input ''n'' must be scalar')
      %NOTICE: n is the size of the smoothing/conv kernel, so edges with width (n-1)/2 are appended
      w=(num.odd(n)-1)/2;
      %upper edge
      map=[flipud(map(2:w+1,:));map];
      %lower edge
      map=[map;flipud(map(end-w:end-1,:))];
      %left edge
      map=[map(:,end-w+1:end),map];
      %right edge
      map=[map,map(:,w+1:2*w)];
    end
    %% catchment plotting
    function catchment=catchment_plot(catchment,varargin)
      v=varargs.wrap('sources',{....
        {...
          'plot_parametric_components',{},   @iscellstr;...
          'parametric_components_line_fmt',{'--'},@ischar;...
          'plot_lines_over_gaps_narrower_than', days(120), @isduration;...
          'time',[],@isdatetime;...
        },...
      },varargin{:});
      %add gaps (trivial call handled inside)
      if isempty(v.time)
        ct=catchment.ws.addgaps(v.plot_lines_over_gaps_narrower_than);
      else
        %align with requested time domain
        ct=catchment.ws.interp(v.time,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than);
      end
      %plot it; NOTICE: for plots with more than 7 lines, matlab's default color scheme is
      %too short, so you need to externally define a colormap and, for each line, pass the
      %argument pair: 'color',[x,y,z] (which is handled in simpledata.plot)
      catchment.plot=ct.plot(varargin{:});
      %add parametric decompositions (if requested)
      if isfield(catchment,'pws')
        %accumulate components to plot
        for i=1:numel(v.plot_parametric_components)
          if isempty(v.time)
            pwst_now=catchment.pws.(['ts_',v.plot_parametric_components{i}]);
          else
            pwst_now=pardecomp.join(catchment.pws,'time',v.time,'coeffnames',v.plot_parametric_components(i));
          end
          if i==1
            catchment.pwst=pwst_now;
          else
            catchment.pwst=catchment.pwst+pwst_now;
          end
          catchment.pwst_val.(v.plot_parametric_components{i})=catchment.pws.(v.plot_parametric_components{i}).y_masked;
        end
        %plot it
        catchment.pwst.plot(varargin{:},'line',v.parametric_components_line_fmt)
        %copy the color from the last line (and get its handle)
        llh=plotting.line_color_set_last;
        %don't show this line in the legend
        set(get(get(llh,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
      end
      title(catchment.name)
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
    function a=test(method,l)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      if ~exist('l','var') || isempty(l)
        l=[5,7,2];
      end
      switch lower(method)
      case 'all'
       for i={'reps','print','resample',...
         'plus','minus','times','rdivide',...
         'area','height','volume',...
         'crop','interp2',...
         'min','max','mean','std','rms','meanabs','stdabs','rmsabs','length','gaps',...
         'coast','plot','conv'}
         simplegrid.test(i{1},l);
       end
      case 'reps'
        dt={'vecmat','matrix','list','flatlist'};
        a=struct(...
          'lon',simplegrid.lon_default(l(2)),...
          'lat',simplegrid.lat_default(l(1)),...
          'map',reshape(1:prod(l),l),...
          't',datetime(sort(datenum(round(rand(l(3),6)*60))),'ConvertFrom','datenum')...
        );
        dd={...
          a,...
          simplegrid.vecmat2matrix(a),...
          simplegrid.vecmat2list(a),...
          simplegrid.vecmat2flatlist(a)...
        };
        for i=1:numel(dt)
          for j=1:numel(dt)
            out=simplegrid.dtc(dt{i},dt{j},dd{i});
            if any(any(out.map~=dd{j}.map))
              error([mfilename,': failed data type conversion between ''',dt{i},''' and ''',dt{j},'''.'])
            end
          end
        end
      case 'print'
        simplegrid.slanted(l(1),l(2)).print
      case 'resample'
        a=simplegrid.slanted(l(1),l(2));
        figure
        subplot(1,2,1)
        a.imagesc('t',a.t);
        title('original')
        subplot(1,2,2)
        a.spatial_resample(l(1)*3,l(2)*2).imagesc('t',a.t);
        title('resampled')
      case {'plus','minus','times','rdivide'}
        a=simplegrid.unit_randn(l(1),l(2),'t',time.rand(l(3)));
        b=simplegrid.unit(      l(1),l(2),'t',a.t,'scale',2);
        switch lower(method)
          case 'plus'
            c=a+b;
          case 'minus'
            c=a-b;
          case 'times'
            c=a.*b;
          case 'rdivide'
            c=a./b;
        end
        disp('a:')
        disp(a.map)
        disp('b:')
        disp(b.map)
        disp(['a ',lower(method),' b:'])
        disp(c.map)
      case {'area','height','volume'}
        a=simplegrid.unit(l(1),l(2),'t',time.rand(l(3)));
        disp('original')
        a.prettymap(a.t(1))
        disp(method)
        a.(method).prettymap(a.t(1))
      case 'crop'
        %TODO: update this test once the interpolant method works with multiple epochs
        %a=simplegrid.unit_randn(l(1),l(2),'t',time.rand(l(3)));
        a=simplegrid.unit_randn(l(1),l(2),'t',time.rand(1));
        disp('original')
        a.prettymap(a.t(1))
        disp('cropped')
        a.lat=transpose(-45:15:60);
        a.lon=45:45:270;
        a.prettymap(a.t(1))
      case 'interp2'
        %TODO: update this test once the interpolant method works with multiple epochs
        %a=simplegrid.unit(l(1),  l(2),  't',datetime('now')+days(  1:1:2),'scale',0.5);
        %b=simplegrid.unit(l(1)-1,l(2)-1,'t',datetime('now')+days(0.5:1:2),'scale',2);
        a=simplegrid.unit(l(1),  l(2),  't',datetime('now'),'scale',0.5);
        b=simplegrid.unit(l(1)-1,l(2)-1,'t',datetime('now'),'scale',2);
        disp('originals:a')
        a.prettymap
        disp('originals:b')
        b.prettymap
        [a,b]=interp2(a,b);
        disp('interpolated:a')
        a.prettymap
        disp('interpolated:b')
        b.prettymap
      case {'min','max','mean','std','rms','meanabs','stdabs','rmsabs','length','gaps'}
        a=simplegrid.unit_randn(l(1),l(2),'t',time.rand(l(3)));
        m=a.map;
        m(1,1,:)=0;
        m(round(l(1)/2),round(l(2)/3),:)=1;
        a.map=m;
        disp('original')
        a.prettymap
        disp(method)
        a.stats('mode',method,'period',days(inf)).prettymap
      case 'coast'
        a=simplegrid.coast;
      case 'plot'
        a=simplegrid.slanted(l(1)*10,l(2)*10).imagesc;
      end
    end
  end
  methods
    %% constructor
    function obj=simplegrid(t,map,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('t',  @(i) ischar(i) || isnumeric(i) || isdatetime(i)); %this can be char, double or datetime
      p.addRequired('map',@(i) isnumeric(i) || iscell(i));
      p.addParameter('lat',simplegrid.lat_default(size(map,1)), @isnumeric);
      p.addParameter('lon',simplegrid.lon_default(size(map,2)), @isnumeric);
      %create argument object, declare and parse parameters, save them to obj
      v=varargs.wrap('parser',p,'sources',{simplegrid.parameters('obj')},'mandatory',{t,map},varargin{:});
      % retrieve structure with list
      sl=simplegrid.dti(t,p.Results.map,p.Results.lon,p.Results.lat,'list');
      % call superclass
      vararginnow=cells.vararginclean(varargin,{'lat','lon'});
      obj=obj@simpletimeseries(sl.t,sl.map,'lon',sl.lon,'lat',sl.lat,vararginnow{:});
      % save lon and lat
      obj.loni=sl.lon;
      obj.lati=sl.lat;
      % save the arguments v into this object
      obj=v.save(obj,{'t','map','lat','lon'});
    end
    function obj=assign(obj,map,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'map',         @(i) isnumeric(map) || iscell(map));
      p.addParameter('t',   obj.t,  @(i) simpletimeseries.valid_t(i));
      p.addParameter('lon', obj.lon,@(i) isnumeric(i) && ~isempty(i));
      p.addParameter('lat', obj.lat,@(i) isnumeric(i) && ~isempty(i));
      % parse it
      p.parse(map,varargin{:});
      % simpler names
      presence=simpletimeseries.ispresent(p);
      % handle time dependency: branch on presence of arguments 't' and 'x'
      if presence.t
        t=p.Results.t;
      elseif presence.x
        t=p.Unmatched.x;
      elseif ~isempty(obj.t)
        t=obj.t;
      else
        error([mfilename,': cannot assign a map without either ''x'' or ''t''.'])
      end
      % If inputs p.Results.lon and p.Results.lat are empty, it's because
      % obj.lati and obj.loni are too. This test makes sense because
      % simpledata will call the assign method before the lon/lat fields
      % are assigned.
      if ~isempty(p.Results.lon) && ~isempty(p.Results.lat)
        %checking if lon/lat are in the argument list
        if any(strcmp(p.UsingDefaults,'lon')) && ...
           any(strcmp(p.UsingDefaults,'lat'))
          % the format may not be properly determined, so assume map is given in list-type
        else
          % lon/lat given in inputs, retrieve structure with list (format is determined in simplegrid.dti)
          sl=simplegrid.dti(t,map,p.Results.lon,p.Results.lat,'list');
          % save lon and lat
          obj.loni=sl.lon;
          obj.lati=sl.lat;
          % propagate map to be assigned below
          map=sl.map;
        end
      else
        %sanity
        if ~ismatrix(map)
          error([mfilename,': if lon/lat are empty, can only deal with matrix-type maps.'])
        end
        if numel(t) ~= size(map,1)
          error([mfilename,': input ''t'' or ''x'' must have the same nr of row as ''map''.'])
        end
      end
      % add width-reseting flag, if needed
      if obj.width~=size(map,2)
        if size(map,2)==numel(obj.lati) %or obj.loni
          varargin{end+1}='reset_width';
          varargin{end+1}=true;
        else
          error([mfilename,': trying to assign a map that is not consistent with existing lon/lat domain. ',...
            'Update that first or pass them as arguments to this member.'])
        end
      end
      % pass it upstream
      obj=assign@simpletimeseries(obj,map,varargin{:});
    end
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in,[simplegrid.parameters('list');more_parameters(:)]);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      out=metadata@simpletimeseries(obj,[simplegrid.parameters('list');more_parameters(:)]);
    end
    %% info methods
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      relevant_parameters={'lat_domain','lon_domain','sp_units','sp_tol'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print superclass
      print@simpletimeseries(obj,tab)
    end
    function prettymap(obj,t)
      if ~exist('t','var') || isempty(t)
        t=obj.t;
      end
      out=nan(obj.lat_length+1,obj.lon_length+1);
      out(1,2:end)=obj.lon;
      out(2:end,1)=obj.lat;
      for i=1:numel(t)
        out(1,1)=i;
        out(2:end,2:end)=obj.interp(t(i)).map;
        disp(t(i))
        disp(out)
      end
    end
    %NOTICE: these are temporal statistics
    function out=stats(obj,varargin)
      %call upstream method
      s=stats@simpletimeseries(obj,varargin{:});
      %retrieve domains
      lat_now=obj.lat;
      lon_now=obj.lon;
      %branch on the type of output
      switch class(s)
      case 'double'
        %create new object
        out=simplegrid(...
          mean(obj.t),...
          reshape(s,numel(lat_now),numel(lon_now)),...
          'lat',lat_now,...
          'lon',lon_now...
        );
      case 'simpletimeseries'
        %translate new object
        out=simplegrid(...
          s.t,...
          reshape(s.y,numel(lat_now),numel(lon_now),s.t),...
          'lat',lat_now,...
          'lon',lon_now...
        );
      case 'struct'
        error([mfilename,': cannot retrive multiple stats in the form a structure.'])
      otherwise
        error([mfilename,': cannot handle stats when upstream method returns class ',class(s),'.'])
      end
      %copy metadata
      out=out.copy_metadata(obj);
      %update descriptor
      out.descriptor=[varargs(varargin).mode,' of ',obj.descriptor];
    end
    %% (spatial) units convertion
    function obj=convert(obj,units,scale)
      switch lower(units)
      case {'rad','radian','radians'}
        switch obj.sp_units
        case 'deg'; scale=deg2rad(1);
        case 'rad'; scale=1;
        end
      case {'deg','degree','degrees'}
        switch obj.sp_units
        case 'deg'; scale=1;
        case 'rad'; scale=rad2deg(1);
        end
      end
      %trivial all
      if scale==1,return,end
      %scale maps
      obj=obj.scale(scale);
      %scale the lat/lon domains
      obj.lati=obj.lati*scale;
      obj.loni=obj.loni*scale;
      %update units
      obj.sp_units=units;
    end
    %% vecmat handling
    function sv=get.vecmat(obj)
      sv=simplegrid.dtc('list','vecmat',obj.list);
    end
    function obj=set.vecmat(obj,vecmat)
      simplegrid.vecmat_check(vecmat);
      obj=obj.assign(vecmat.map,...
                 't',vecmat.t,...
               'lat',vecmat.lat,...
               'lon',vecmat.lon...
      );
    end
    %% list handling 
    function sl=get.list(obj)
      sl=simplegrid.list_init(obj.t,obj.y,obj.loni,obj.lati);
    end
    function obj=set.list(obj,list)
      simplegrid.list_check(list);
      obj=obj.assign(list.map,...
                 't',list.t,...
               'lat',list.lat,...
               'lon',list.lon...
      );
    end
    %% flatlist handling 
    function sl=get.flatlist(obj)
      sl=simplegrid.dtc('list','flatlist',obj.list);
    end
    function obj=set.flatlist(obj,flatlist)
      simplegrid.flatlist_check(flatlist);
      obj=obj.assign(flatlist.map,...
                 't',flatlist.t,...
               'lat',flatlist.lat,...
               'lon',flatlist.lon...
      );
    end
    %% matrix handling
    function sm=get.matrix(obj)
      sm=simplegrid.dtc('vecmat','matrix',obj.vecmat);
    end
    function obj=set.matrix(obj,matrix)
      simplegrid.matrix_check(matrix);
      obj=obj.assign(matrix.map,...
                 't',matrix.t,...
               'lat',matrix.lat,...
               'lon',matrix.lon...
      );
    end
    %% map handling
    function out=get.map(obj)
      %convert map (internally in list) to vecmat
      sv=simplegrid.list2vecmat(simplegrid.list_init(obj.t,obj.y,obj.loni,obj.lati));
      %output
      out=sv.map;
    end
    function obj=set.map(obj,map)
      obj.vecmat=simplegrid.vecmat_init(obj.t,map,obj.lon,obj.lat);
    end
    function out=map_weighted(obj)
      w=obj.lat_weights;
      out=obj.map;
      for i=1:obj.length
        out(:,:,i)=w.*out(:,:,i);
      end
    end
    function out=map_sum(obj);                  out=sum(sum(obj.map));  end
    function out=map_sum_weighted(obj);         out=sum(sum(obj.map_weighted));  end
    function out=map_squared_sum(obj);          out=sum(sum(obj.map.^2));  end
    function out=map_squared_sum_weighted(obj)
      w=obj.lat_weights;
      out=obj.map;
      for i=1:obj.length
        out(:,:,i)=w.*out(:,:,i).^2;
      end
      out=sum(sum(out));
    end
    %% lat handling
    function obj=set.lat(obj,lat_now)
      %trivial call
      if numel(lat_now)==numel(obj.lat) && all(lat_now==obj.lat)
        return
      end
      obj=obj.spatial_interp(obj.lon,lat_now);
    end
    function out=get.lat(obj)
      if isempty(obj.lati)
        out=[];
      else
        out=obj.lat_convert('vecmat');
      end
    end
    function out=lat_convert(obj,type)
      try
        out=obj.(type).lat;
      catch
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    function out=get.latSpacing(obj)
      out=simplegrid.spacing(obj.lat);
    end
    function obj=set.latSpacing(obj,spacing)
      obj.lat=simplegrid.lat_default(181/spacing);
    end
    function out=lat_domain(obj)
      out=['[ ',...
        num2str(obj.lat(1)),' : ',...
        num2str(obj.latSpacing),' : ',...
        num2str(obj.lat(end)),...
      ' ]'];
    end
    function out=lat_length(obj)
      out=numel(obj.lat);
    end
    function out=lat_rad(obj)
      switch obj.sp_units
      case 'deg'; out=deg2rad(obj.lat);
      case 'rad'; out=obj.lat;
      end
    end
    function out=lat_weights(obj)
      out=(cos(obj.lat_rad))*ones(size(obj.lon));
    end
    %% lon handling
    function obj=set.lon(obj,lon_now)
      %trivial call
      if numel(lon_now)==numel(obj.lon) && all(lon_now==obj.lon)
        return
      end
      obj=obj.spatial_interp(lon_now,obj.lat);
    end
    function out=get.lon(obj)
      if isempty(obj.loni)
        out=[];
      else
        out=obj.lon_convert('vecmat');
      end
    end
    function obj=set.lon360(obj,lon_now)
      obj.lon=simplegrid.setlon360(lon_now);
    end
    function out=get.lon360(obj)
      out=simplegrid.setlon360(obj.lon);
    end
    function obj=set.lon180(obj,lon_now)
      obj.lon=wrapTo180(lon_now);
    end
    function out=get.lon180(obj)
      out=wrapTo180(obj.lon);
    end
    function out=lon_convert(obj,type)
      try
        out=obj.(type).lon;
      catch
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    function out=get.lonSpacing(obj)
      out=simplegrid.spacing(obj.lon);
    end
    function obj=set.lonSpacing(obj,spacing)
      obj=obj.lon_set(simplegrid.lon_default(300/spacing));
    end
    function out=lon_domain(obj)
      out=['[ ',...
        num2str(obj.lon(1)),' : ',...
        num2str(obj.lonSpacing),' : ',...
        num2str(obj.lon(end)),...
      '] '];
    end
    function out=lon_length(obj)
      out=numel(obj.lon);
    end
    function out=lon_rad(obj)
      switch obj.sp_units
      case 'deg'; out=deg2rad(obj.lon);
      case 'rad'; out=obj.lon;
      end
    end
    function [obj,changed]=add360lon(obj)
      %see if things need changing
      changed=~simplegrid.islon360complete(obj.lon);
      if changed
        %get map
        m=obj.map;
        %append longitude 360
        m(:,end+1,:)=m(:,1,:);
        %TODO: make this work for time entries (multiple maps)
        obj=obj.assign(m,'lon',[obj.lon360,360]);
      end
    end
    function [obj,changed]=del360lon(obj)
      %assume nothing changes
      changed=false;
      %easier names
      l360=obj.lon360;
      %check if lon 360 is not there
      if (l360(end)-360)^2>simplegrid.parameters('sp_tol')^2
        return
      end
      %check if current spatial domain includes longitude 0
      if (l360(1))^2>simplegrid.parameters('sp_tol')^2
        return
      end
      %if not, then things changed
      changed=true;
      %crop longitude 360
      m=obj.map;
      %TODO: make this work for multiple time entries
      obj=obj.assign(m(:,1:end-1),'lon',obj.lon360(1:end-1));
    end
    %% spatial domain
    function out=sp_limits(obj)
      out=simplegrid.getLimits(obj.lon,obj.lat);
    end
    function out=sp_spacing(obj)
      out=simplegrid.getSpacing(obj.lon,obj.lat);
    end
    function out=isglobal(obj)
      out=all( (obj.sp_limits-[0 360 -90 90]).^2 < simplegrid.parameters('sp_tol')^2 );
    end
    %% interpolant handling
    function [I,obj]=interpolant(obj)
      [lat_meshed,lon_meshed,t_meshed]=ndgrid(obj.lat,obj.lon,datenum(obj.t));
      %TODO: handle multiple time entries
      if obj.length==1
        I=griddedInterpolant(lat_meshed,lon_meshed,obj.map,'linear','none');
      else
        error('look in the code: need to debug the wrapping of obj.map below (so it looks like the above)')
        I=griddedInterpolant(lat_meshed,lon_meshed,t_meshed,obj.map,'linear','none');
      end
    end
    %% spatial interpolation
    function obj=spatial_interp(obj,lon_new,lat_new)
      if ~iscolumn(lat_new)
        error([mfilename,': illegal input ''lat''.'])
      end
      if ~isrow(lon_new)
        error([mfilename,': illegal input ''lon''.'])
      end
      %build new meshed domain
      [lat_meshed,lon_meshed,t_meshed]=ndgrid(lat_new,lon_new,datenum(obj.t));
      %fix missing 360 (if needed)
      [obj,changed]=obj.add360lon;
      %get the interpolant
      [I,obj]=obj.interpolant;
      %interpolate
      if obj.length==1
        map_new=I(lat_meshed,lon_meshed);
      else
        map_new=I(lat_meshed,lon_meshed,t_meshed);
      end
      %patch nans
      if any(isnan(map_new(:)))
        %fixing known cases of stupid nans
        if all(isnan(map_new(end,:)))
          map_new(end,:)=interp1(obj.lon,obj.map(end,:),lon_new);
        end
        if all(isnan(map_new(:,end)))
          map_new(:,end)=interp1(obj.lat,obj.map(:,end),lat_new);
        end
        if isnan(map_new(end))
          submap=map_new(end-1:end,end-1:end);
          map_new(end)=mean(submap(~isnan(submap)));
        end
        assert(~any(isnan(map_new(:))),'BUG TRAP:Still there''s NaNs in the interpolated map')
      end
      %save results
      obj=obj.assign(map_new,'t',obj.t,'lat',lat_new,'lon',lon_new);
      %unfix missing 360 (if needed)
      if changed && ~simplegrid.islon360complete(lon_new)
        obj=obj.del360lon;
      end
    end
    function obj=spatial_resample(obj,lon_n,lat_n)
      obj=obj.spatial_interp(simplegrid.lon_default(lon_n),simplegrid.lat_default(lat_n));
    end
    function obj=center_resample(obj)
      %fix missing 360 (if needed)
      obj=obj.add360lon;
      %get mean spatial domain
      lon_new=mean([obj.lon(1:end-1);obj.lon(2:end)]);
      lat_new=mean([obj.lat(1:end-1),obj.lat(2:end)],2);
      %reduced to mean domain
      obj=obj.spatial_interp(lon_new,lat_new);
    end
    %% spatial cropping
    function obj=spatial_crop(obj,lon_limits,lat_limits)
      %get the data
      m=obj.matrix;
      %get the spatial mask
      mask=struct(...
        'lat',obj.lat>=lat_limits(1) & obj.lat<=lat_limits(2), ...
        'lon',obj.lon360>=simplegrid.setlon360(lon_limits(1)) & obj.lon360<=simplegrid.setlon360(lon_limits(2)) ...
      );
%         'lon',wrapTo180(obj.lon)>=wrapTo180(lon_limits(1)) & wrapTo180(obj.lon)<=wrapTo180(lon_limits(2)) ...
      %apply the mask
      m.lat=m.lat(mask.lat,mask.lon);
      m.lon=m.lon(mask.lat,mask.lon);
      m.map=m.map(mask.lat,mask.lon,:);
      %propagate the data
      obj.matrix=m;
      %fix labels and units (NOTICE; this assumes there is no order in the labels and units (which makes sense in a grid)
      obj.labels=obj.labels(1:obj.width);
      obj.units=obj.units(1:obj.width);
    end
    function obj=spatial_mask(obj,mode,varargin)
      %NOTICE: a positive buffer reduces the land/ocean areas, a negative buffer enlarges it
      %trivial call
      if str.none(mode); return; end
      %handle special cases
      switch lower(mode)
      case 'deep ocean'
        obj=obj.spatial_mask('ocean','buffer',1000e3);
        return
      case {'land','ocean'}
        %handle optionals
        v=varargs.wrap('sources',{{...
          'buffer',   0, @num.isscalar;...
          'cutoff', 0.1,@(i) (num.isscalar(i)) || isempty(i);... %NOTICE: this value seems to work well
        }},varargin{:});
        %retrieve the land mask
        spmask=simplegrid.landmask(obj.lon,obj.lat,varargin{:});
        %apply buffer
        if v.buffer~=0
          spmask=spmask.sh(simplegrid.sh_lmax(obj.lon,obj.lat)...
          ).scale(abs(v.buffer),'gauss'...
          ).grid(varargin{:}...
          ).spatial_interp(obj.lon,obj.lat...
          );
          obj.descriptor=[num2str(v.buffer*1e-3),'km buffered ',mode,' areas of ',obj.descriptor];
        else
          obj.descriptor=[mode,' areas of ',obj.descriptor];
        end
        %enforce cutoff
        if ~isempty(v.cutoff)
           y_now=spmask.y;
           switch lower(mode)
           case 'land';  if v.buffer>0; idx=y_now>1-v.cutoff; else idx=y_now>v.cutoff; end
           case 'ocean'; if v.buffer<0; idx=y_now>1-v.cutoff; else idx=y_now>v.cutoff; end
           end
           y_now(idx)=1;
           y_now(~idx)=0;
           spmask=spmask.assign(y_now);
        end
        %enforce mode
        switch lower(mode)
        case 'land';  %do nothing
        case 'ocean'; spmask=spmask.assign(1-spmask.y);
        otherwise; error(['Cannot handle the spatial mask ''',mode,'''.'])
        end
      case {'tropical','non-tropical','polar','non-polar'}
        spmask=simplegrid.unit(obj.lon,obj.lat,varargin{:});
        latmask=true(size(obj.lat));
        switch lower(mode)
        case 'tropical'    ;latmask(abs(obj.lat)> 23.43671)=false;
        case 'non-tropical';latmask(abs(obj.lat)< 23.43671)=false;
        case 'polar'       ;latmask(abs(obj.lat)< 60      )=false;
        case 'non-polar'   ;latmask(abs(obj.lat)> 60      )=false;
        end
        spmap=spmask.map;
        spmap(~latmask,:)=0;
        spmask.map=spmap;
      otherwise
        error(['Unknown spatial mask with name ''',mode,'''.'])
      end
      %apply mask
      obj=obj.*spmask;
    end
    %% spatial ops
    function [out,m]=global_spatial_op(obj,mode)
      %NOTICE: this returns a time series, the spatial domain is reduced to a point
      %get the data
      m=obj.vecmat;
      %collapse data to a single spatial point
      m.lat=mean(m.lat);
      m.lon=mean(m.lon);
      m.t=obj.t;
      m.n=sum(sum(obj.lat_weights));
      %branch on mode
      switch lower(mode)
        case 'sum';          m.map=obj.map_sum_weighted;
        case 'mean';         m.map=obj.map_sum_weighted/m.n;
        case 'squared-sum';  m.map=obj.map_squared_sum_weighted;
        case 'rms';          m.map=sqrt(obj.map_squared_sum_weighted/m.n);
        case {'rss','norm'}; m.map=sqrt(obj.map_squared_sum_weighted);
      end
      %reduce to timeseries object
      out=simpletimeseries(obj.t,m.map(:)).copy_metadata(obj);
      out.descriptor=[mode,' of ',obj.descriptor,' (centered at ',num2str(m.lat),'deg lat by ',num2str(m.lon),'deg long)'];
      %adjust labels and units (NOTICE; this assumes there is no order in the labels and units (which makes sense in a grid)
      out.labels=obj.labels(1);
      out.units=obj.units(1);
    end
    %% utilities
    function obj=area(obj,radius)
      if ~exist('radius','var') || isempty(radius)
        radius=gravity.parameters('Rm');
      end
      original_units=obj.sp_units;
      obj=obj.convert('rad');
      lat1=obj.lat(1:end-1);
      lat2=obj.lat(2:end);
      lon1=obj.lon(1:end-1);
      lon2=obj.lon(2:end);
      map_now=abs( (sin(lat1)-sin(lat2)) * (lon1-lon2) )*radius^2;
      %new lat/lon domain
      lon_new=mean([lon1;lon2]);
      lat_new=mean([lat1,lat2],2);
      %save results
      map_c=cell(1,obj.length);
      map_c(:)={map_now};
      obj=obj.assign(map_c,'t',obj.t,'lat',lat_new,'lon',lon_new);
      %convert back to original units
      obj=obj.convert(original_units);
    end
    function obj=height(obj)
      lon_new=mean([obj.lon(1:end-1);obj.lon(2:end)]);
      lat_new=mean([obj.lat(1:end-1),obj.lat(2:end)],2);
      m=obj.map;
      %save results
      obj=obj.assign(...
        (...
          m(1:end-1,1:end-1,:)+...
          m(2:end  ,1:end-1,:)+...
          m(1:end-1,2:end  ,:)+...
          m(2:end  ,2:end  ,:)...
        )*0.25,...
        't',obj.t,...
        'lat',lat_new,...
        'lon',lon_new...
        );
    end
    function obj=volume(obj)
      obj=obj.area.*obj.height;
    end
    function [y,x,ts]=rms(obj,mode)
      h=obj.map;
      w=obj.lat_weights;
      switch mode
      case {'total','tot',0}
        y=zeros(size(h,3),1);
        for i=1:obj.length
          y(i)=num.rms(h(:,:,i),w,0);
        end
        x=0;
      case {'lat','row',1}
        y=zeros(size(h,3),size(h,2));
        for i=1:obj.length
          y(i)=num.rms(h(:,:,i),w,1);
        end
        x=h.lon;
      case {'lon','col',2}
        y=zeros(size(h,3),size(h,1));
        for i=1:obj.length
          y(i)=num.rms(h(:,:,i),w,2);
        end
        x=h.lat;
      otherwise
        error(['Cannot understand mode ''',mode,'''.'])
      end
      if nargout>1
        switch mode
        case {'total','tot',0}
          labels={'Global RMS'};
        case {'lat','row',1}
          labels=arrayfun(@(i) {str.show({'RMS at',i,'lon'})},x,'UniformOutput',false);
        case {'lon','col',2}
          labels=arrayfun(@(i) {str.show({'RMS at',i,'lat'})},x,'UniformOutput',false);
        end
        ts=simpletimeseries(...
          obj.t,...
          y,...
          'format','datetime',...
          'labels',labels,...
          'timesystem',obj.timesystem,...
          'units',obj.units(1:numel(x)),...
          'descriptor',['RMS of ',obj.descriptor]...
        );      
      end
    end
    function [y,x,ts]=mean(obj,mode)
      h=obj.map;
      w=obj.lat_weights;
      switch mode
      case {'total','tot',0}
        y=zeros(size(h,3),1);
        for i=1:obj.length
          y(i)=num.mean(h(:,:,i),w,0);
        end
        x=0;
      case {'lat','row',1}
        y=zeros(size(h,3),size(h,2));
        for i=1:obj.length
          y(i)=num.mean(h(:,:,i),w,1);
        end
        x=h.lon;
      case {'lon','col',2}
        y=zeros(size(h,3),size(h,1));
        for i=1:obj.length
          y(i)=num.mean(h(:,:,i),w,2);
        end
        x=h.lat;
      otherwise
        error(['Cannot understand mode ''',mode,'''.'])
      end
      if nargout>1
        switch mode
        case {'total','tot',0}
          labels={'Global mean'};
        case {'lat','row',1}
          labels=arrayfun(@(i) {str.show({'mean at',i,'lon'})},x,'UniformOutput',false);
        case {'lon','col',2}
          labels=arrayfun(@(i) {str.show({'mean at',i,'lat'})},x,'UniformOutput',false);
        end
        ts=simpletimeseries(...
          obj.t,...
          y,...
          'format','datetime',...
          'labels',labels,...
          'timesystem',obj.timesystem,...
          'units',obj.units(1:numel(x)),...
          'descriptor',['mean of ',obj.descriptor]...
        );      
      end
    end
    function [y,x,ts]=std(obj,mode)
      h=obj.map;
      w=obj.lat_weights;
      switch mode
      case {'total','tot',0}
        y=zeros(size(h,3),1);
        for i=1:obj.length
          y(i)=num.std(h(:,:,i),w,0);
        end
        x=0;
      case {'lat','row',1}
        y=zeros(size(h,3),size(h,2));
        for i=1:obj.length
          y(i)=num.std(h(:,:,i),w,1);
        end
        x=h.lon;
      case {'lon','col',2}
        y=zeros(size(h,3),size(h,1));
        for i=1:obj.length
          y(i)=num.std(h(:,:,i),w,2);
        end
        x=h.lat;
      otherwise
        error(['Cannot understand mode ''',mode,'''.'])
      end
      if nargout>1
        switch mode
        case {'total','tot',0}
          labels={'Global STD'};
        case {'lat','row',1}
          labels=arrayfun(@(i) {str.show({'STD at',i,'lon'})},x,'UniformOutput',false);
        case {'lon','col',2}
          labels=arrayfun(@(i) {str.show({'STD at',i,'lat'})},x,'UniformOutput',false);
        end
        ts=simpletimeseries(...
          obj.t,...
          y,...
          'format','datetime',...
          'labels',labels,...
          'timesystem',obj.timesystem,...
          'units',obj.units(1:numel(x)),...
          'descriptor',['STD of ',obj.descriptor]...
        );      
      end
    end
    %% convert to spherical harmonics
    %it's a good idea to identify the functional represented in this grid explicitly in varargin ('functional',<something>)
    function out=sh(obj,lmax,varargin)
%TODO: I'm not sure why this is needed; when plotting the std of a GRACE time series,
%there was a strong vertical artefact at zero longitude
%       % make sure longitude domain goes up to 360 degree
%       obj=obj.add360lon;
%       % need global domain
%       assert(obj.isglobal,'Need global domain to perform SH analysis')
      % make room for CS-structures
      cs(obj.length)=struct('C',[],'S',[]);
      % make spherical harmonic analysis of each grid
      c=0;
      for i=1:obj.length
        if ~obj.mask(i)
          str.say('WARNING: discarding gap of',obj.descriptor,'at',obj.t(i))
          continue
        end
        c=c+1;
        [cs(c).C,cs(c).S,msg]=mod_sh_ana(obj.lon_rad,obj.lat_rad,obj.map(:,:,i),lmax);
        if ~isempty(msg); str.say(msg,'at',obj.t(i)); end
      end
      % initialize output gravity object (no need to set the units, need 'functional' in vargain for that to work)
      out=gravity.unit(lmax,...
        't',obj.t_masked,...
        'descriptor',obj.descriptor,...
      varargin{:});
      %assign coefficients (truncate to the masked domain)
      out.cs=cs(1:c);
    end
    %% multiple operands
    function compatible(obj1,obj2,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('skip_par_check',{''},@iscellstr)
      p.parse(varargin{:});
      %call mother routine
      compatible@simpledata(obj1,obj2,varargin{:});
      %shorter names
      par=simplegrid.compatible_parameter_list;
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
    %the merge method can be called irectly
    function [lon,lat]=sp_domain_lcm(obj1,obj2)
      %get parameters
      lim1=obj1.sp_limits;
      spa1=obj1.sp_spacing;
      lim2=obj2.sp_limits;
      spa2=obj2.sp_spacing;
      %build longitude domain
      lon=(...
        min([lim1(1),lim2(1)]):...
        gcd( spa1(1),spa2(1) ):...
        max([lim1(2),lim2(2)])...
      );
      %build latitude domain
      lat=transpose(...
        min([lim1(3),lim2(3)]):...
        gcd( spa1(2),spa2(2) ):...
        max([lim1(4),lim2(4)])...
      );
    end
    function [obj1,obj2]=interp2(obj1,obj2,varargin)
      %get common spatial domain
      [lon_new,lat_new]=sp_domain_lcm(obj1,obj2);
      %interpolate in the spatial domain
      obj1=obj1.spatial_interp(lon_new,lat_new);
      obj2=obj2.spatial_interp(lon_new,lat_new);
      %interpolate in the time domain
      [obj1,obj2]=interp2@simpletimeseries(obj1,obj2);
    end
    %the append method can be called directly (only acts in the time domain)
    %% plot functions
    function out=imagesc(obj,varargin)
      v=varargs.wrap('sources',{{...
      't',obj.t_masked([],'start'), @(i) simpletimeseries.valid_t(i);...
      'center_resample',      true, @islogical;...
      'show_coast',           true, @islogical;...
      'show_colorbar',        true, @islogical;...
      'plot_spatial_step',       0, @num.isscalar;...
      'cb_loc',     'SouthOutside', @ischar;...
      'cb_title',               '', @ischar;...
      'bias',                    0, @num.isscalar;...
      'boxes',                  {}, @iscell;...s
      'boxes_fmt',         {'r--'}, @iscellstr;...
      }},varargin{:});
      %interpolate at the requested time and 
      obj_interp=obj.interp(v.t);
      %upsample if needed
      if v.plot_spatial_step>0
        obj_interp=spatial_interp(obj_interp,simplegrid.lon_stepped(v.plot_spatial_step),simplegrid.lon_stepped(v.plot_spatial_step));
      end
      %resample to center of grid if requested 
      if v.center_resample
        obj_interp=obj_interp.center_resample;
      end
      %need to have lon domain in the -180 to 180 domain, so that coast is show properly
      [lon_now,lon_idx]=sort(obj_interp.lon180);
      %build image
      out.axis_handle=imagesc(...
        lon_now,...
        obj_interp.lat,...
        obj_interp.map(:,lon_idx,:)...
      );hold on
      axis xy
      axis equal
      axis tight
      %labels
      xlabel(['lon [',obj.sp_units,']'])
      ylabel(['lat [',obj.sp_units,']'])
      %ploting coastline
      if v.show_coast
        out.coast_handle=simplegrid.coast('lon',lon_now);
      end
      %add colorbar
      if v.show_colorbar
        out.colorbar_handle=colorbar('location',v.cb_loc);
        if ~isempty(v.cb_title)
          switch lower(v.cb_loc)
          case {'north','south','northoutside','southoutside'}
            set(get(out.colorbar_handle,'xlabel'),'string',v.cb_title);
          case {'east','west','eastoutside','westoutside'}
            set(get(out.colorbar_handle,'ylabel'),'string',v.cb_title);
          end
        end
      end
      %add boxes
      if ~isempty(v.boxes)
        b=simplegrid.catchment_subset(v.boxes);
        b_fmt=cells.deal(v.boxes_fmt,numel(b));
        for i=1:size(b,1)
          plot(...
  wrapTo180([b{i,3}(1) b{i,3}(2) b{i,3}(2) b{i,3}(1) b{i,3}(1)]),...
            [b{i,2}(1) b{i,2}(1) b{i,2}(2) b{i,2}(2) b{i,2}(1)],...
            b_fmt{i},'LineWidth',3 ...
          );
        end
      end
    end
    %% spatial truncation and parameter decomposition
    function catchment=catchment_get(obj,name,varargin)
      v=varargs.wrap('sources',{....
        {...
          'parametric_decomposition', false, @islogical;...
        },...
      },varargin{:});
      %save name of this catchment
      catchment.name=name;
      %branch between different types of "catchments"
      if any(ismember(lower(simplegrid.catchment_list(:,1)),lower(name)))
        %build array with catchment limits
        catchment.lat=simplegrid.catchment_details(name,'lat');
        catchment.lon=simplegrid.catchment_details(name,'lon');
        %compute storage timeseries on this catchment
        catchment.ws=obj.spatial_crop(catchment.lon,catchment.lat).global_spatial_op('mean');
      else
        %TODO: it would be nice to set lat/lon equal to the mask boundary
        catchment.lat=[];
        catchment.lon=[];
        catchment.ws=obj.spatial_mask(name,varargin{:}).global_spatial_op('mean');
      end
      %add parametric decompositions (if requested)
      if v.parametric_decomposition
        %decompose
        catchment.pws=pardecomp.split(catchment.ws,varargin{:});
      end
    end
    function out=catchment(obj,name,varargin)
      out=simplegrid.catchment_plot(obj.catchment_get(name,varargin{:}),varargin{:});
    end
    %% export 
    function xyz(obj,filename,varargin)
      v=varargs.wrap('sources',{{...
        'header',  'default',   @ischar;...
        'columns', 1:obj.width, @isnumeric;...
        'force',   false,       @islogical;...
        'fmt',     '%10.4f %10.4f %14.8e\n', @ischar;
      }},varargin{:});
      disp([datestr(now),': start exporting ',filename])
      %make sure this directory exists
      assert(file.ensuredir(filename),['Error creating directory of file ',filename,'.'])
      %get extension and filename
      [fd,fn,fe]=fileparts(filename);
      %get the map
      m=obj.vecmat;
      for f=1:size(m.map,3)
        %build file for the current time
        f_now=fullfile(fd,[fn,'.',datestr(obj.t(f),'yyyymmddTHHMMSS'),fe]);
        if ~exist(f_now,'file') || v.force
          %open the file (sanity done inside)
          fid=file.open(f_now,'w');
          %save the data
          s.msg=['exporting ',obj.descriptor,' for ',datestr(obj.t(f),'yyyy-mm-dd HH:MM:SS'),' to file ',f_now];s.n=size(m.map,2);
          for j=1:size(m.map,2)
            for i=1:size(m.map,1)
              fprintf(fid,v.fmt,m.lon(j),m.lat(i),m.map(i,j,f));
            end
            s=time.progress(s,j);
          end
          fclose(fid);
        end
      end
    end
  end
end

%% These routines are the foundation for the grid object and were developed by or in cooperation with Pedro Inácio.
%  They remain here to acknowledge that fact.
function out_grid = grid_constructor(varargin)
% GRID=GRID_CONSTRUCTOR is the function that defines a grid 'object'. More
%   specifically GRID is not a MATLAB object, but simply a structure with
%   standard fields which allow the set of function grid_* to operate on it.
%   Typically objects are much more versatile as they allow the definition
%   public and private fields. This implementation of the grid object lacks
%   these features but is nontheless considered appropriate for most uses.
%
%   There are several possible uses of this function:
%   GRID=GRID_CONSTRUCTOR with no arguments returns an empty grid
%   structure with default parameters.
%
%   GRID=GRID_CONSTRUCTOR(XLIST,YLIST,ZLIST) returns the GRID object for the
%   list of points in xlist,ylist,zlist. The grid limits and spacing are
%   computed from the x- and ylists.
%   GRID=GRID_CONSTRUCTOR(XLIST,YLIST) to get an zero-grid is also possible.
%
%   GRID=GRID_CONSTRUCTOR(...,'Units',UNITS) additionally lets you specify the
%   units of the list of points. UNITS is a 1x3 cell with strings
%   specifying the units of x,y and z values. Eg: {'deg','deg','m'}.
%
%   GRID=GRID_CONSTRUCTOR(...,'Date','yyyy-mm-dd HH:MM:SS') lets you
%   specify the date to which the list of points refers to.
%
%   GRID=GRID_CONSTRUCTOR(...,'Description',DESCRIPTION) lets you specify a
%   description of the data.
%
%   Finally notice that the inputs XLIST,YLIST,ZLIST can be replaced with
%   an existing grid, meaning that one should use grid_constructor to edit
%   the field of an existing structure, eg.
%       a_grid = grid_constructor(a_grid,'Units',{'deg','deg','kg'}

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>

% TODO: Simplify limits and spacing handling.
% TODO: How to handle slightly unregular grids?

  % Defaults
  THRS = 1E-8; % Allow input grid to deviate from regular one by this relative amount.

  % Define default inputs
  p=inputParser;
  % xlist  can be many this, so no check is made.
  p.addOptional('xlist',[],@(x)isgrid(x) | isnumeric(x) | ischar(x) | iscell(x));
  p.addOptional('ylist',[],@(x)isnumeric(x));
  p.addOptional('zlist',[],@(x)isnumeric(x));
  p.addParameter('Date','',@(x)ischar(x));
  p.addParameter('Description','',@(x)ischar(x));
  p.addParameter('Units',{'','',''},@(x)iscell(x));
  p.addParameter('GridLimits',grid_getLimits,@(x)isvector(x) && length(x) == 4)
  p.addParameter('GridSpacing',grid_getSpacing,@(x)isvector(x) && length(x) == 2)

  % Parse inputs
  p.parse(varargin{:});
  % Convert to double precision
  if isnumeric(p.Results.xlist)
      xlist = double(p.Results.xlist);
  else
      xlist = p.Results.xlist;
  end
  ylist = double(p.Results.ylist);
  zlist = double(p.Results.zlist);
  Units = p.Results.Units;
  Date  = p.Results.Date;
  Description = p.Results.Description;
  GridLimits = p.Results.GridLimits;
  GridSpacing = p.Results.GridSpacing;

  % Flag inputs coming from the parser
  % NOTE: Maybe this can be done with the parser.UsingDefaults property
  isDefaultXlist = any(strcmp(p.UsingDefaults,'xlist'));
  isDefaultUnits = any(strcmp(p.UsingDefaults,'Units'));
  isDefaultDate  = any(strcmp(p.UsingDefaults,'Date'));
  isDefaultDescription = any(strcmp(p.UsingDefaults,'Description'));
  isDefaultGridLimits = any(strcmp(p.UsingDefaults,'GridLimits')) || all(isnan(GridLimits));
  isDefaultGridSpacing = any(strcmp(p.UsingDefaults,'GridSpacing')) || all(isnan(GridSpacing));

  if nargin == 0
      % Return the default grid
      out_grid.x   = xlist;
      out_grid.nx  = numel(xlist);
      out_grid.y   = xlist;
      out_grid.ny  = numel(ylist);
      out_grid.Map = zlist;
      out_grid.GridLimits  = grid_getLimits(xlist,ylist);
      out_grid.GridSpacing = grid_getSpacing(xlist,ylist);
      out_grid.Date = Date;
      out_grid.Description = Description;
      out_grid.Units = Units;
      return;

  elseif isnumeric(xlist)
      % Inputs are xlist, ylist and zlist.

      % Get default grid
      out_grid = grid_constructor;

      % if zlist is undefined, relax and set the grid to zeros
      if isempty(zlist)
          zlist = zeros(numel(ylist),numel(xlist));
      end

      if isvector(xlist) && isvector(ylist) && numel(xlist) == numel(ylist) && numel(xlist) == numel(zlist)
          % inputs are lists of points

          % Calculate the spacing and limits
          if isDefaultGridLimits || isempty(GridLimits)
              out_grid.GridLimits = grid_getLimits(xlist(:),ylist(:));
          else
              out_grid.GridLimits = GridLimits;
          end
          if isDefaultGridSpacing || isempty(GridSpacing)
              out_grid.GridSpacing = grid_getSpacing(xlist(:),ylist(:));
          else
              out_grid.GridSpacing = GridSpacing;
          end

          % Create full x and y coordinate vectors.
          if isempty(xlist) && isempty(ylist)
              out_grid.x = [];
              out_grid.y = [];
          elseif isempty(xlist) || isempty(ylist)
              error('%s, ERROR: Grid must have values for x- and y',mfilename)
          else
              % TODO: notice that this doesnt not necessarily span all the
              % points in xlist and ylist. Only if they
              % already defined in an (in)complete regular grid.
              out_grid.x  = out_grid.GridLimits(1):out_grid.GridSpacing(1):out_grid.GridLimits(2);
              out_grid.y  = (out_grid.GridLimits(3):out_grid.GridSpacing(2):out_grid.GridLimits(4))';
          end
          out_grid.nx = numel(out_grid.x);
          out_grid.ny = numel(out_grid.y);

          % Missing data on the provided x,y,z list become NaN's.
          nindx = NaN(size(xlist));
          nindy = NaN(size(ylist));
          for i = 1:out_grid.nx
              %THRS This takes into account numerical errors in a regular
              %input grid
              aux = num_equal(xlist,out_grid.x(i),THRS,THRS);
              nindx(aux) = i;
          end
          if ( sum(isnan(nindx)) ~= 0 )
              error('%s: ERROR: x-coordinate is not regular',mfilename)
          end
          for i = 1:out_grid.ny
              %THRS This takes into account numerical errors in a regular
              %input grid
              aux = num_equal(ylist,out_grid.y(i),THRS,THRS);
              nindy(aux) = i;
          end
          if ( sum(isnan(nindy)) ~= 0 )
              error('%s: ERROR: y-coordinate is not regular',mfilename)
          end
          out_grid.Map = NaN(out_grid.ny,out_grid.nx);
          nind = sub2ind(size(out_grid.Map),nindy,nindx);
          out_grid.Map(nind) = zlist;

      elseif ismatrix_(xlist) && ismatrix_(ylist) && all(size(xlist) == size(ylist)) &&...
             all(size(xlist) == size(zlist))
          % inputs are matrices of coordinates and values

          % Check if grid is regular
          if all(diff(xlist(1,:),1,2) == 0)
              % Transpose everything
              xlist = xlist';
              ylist = ylist';
              zlist = zlist';
          end
          if any(any(diff(xlist,1,2) ~= mean(diff(xlist(1,:),1,2))))
              error('%s, ERROR: x,in mat form, must be a regular grid.',mfilename);
          end
          if any(any(diff(ylist) ~= mean(diff(ylist(:,1)))))
              error('%s, ERROR: y- in mat form must be a regular grid.',mfilename);
          end

          % Calculate the spacing and limits
          if isDefaultGridLimits || isempty(GridLimits)
              out_grid.GridLimits = grid_getLimits(xlist(:),ylist(:));
          else
              out_grid.GridLimits = GridLimits;
          end
          if isDefaultGridSpacing || isempty(GridSpacing)
              out_grid.GridSpacing = grid_getSpacing(xlist(:),ylist(:));
          else
              out_grid.GridSpacing = GridSpacing;
          end

          % Create full x and y coordinate vectors.
          if isempty(xlist) && isempty(ylist)
              out_grid.x = [];
              out_grid.y = [];
          elseif isempty(xlist) || isempty(ylist)
              error('%s, ERROR: Grid must have values for x- and y',mfilename)
          else
              % TODO: notice that this doesnt not necessarily span all the
              % points in xlist and ylist. Only if they
              % already defined in an (in)complete regular grid.
              out_grid.x  = out_grid.GridLimits(1):out_grid.GridSpacing(1):out_grid.GridLimits(2);
              out_grid.y  = (out_grid.GridLimits(3):out_grid.GridSpacing(2):out_grid.GridLimits(4))';
          end
          out_grid.nx = numel(out_grid.x);
          out_grid.ny = numel(out_grid.y);

          % Copy zlist to map variable.
          out_grid.Map = zlist;

      elseif all([numel(ylist) numel(xlist)] == size(zlist))
          % inputs are set of x- y- coordinates and mat of values

          % Copy values to the empty grid and check if grid is consistent
          if isvector(xlist)
              if size(xlist,1) > 1
                  out_grid.x = xlist';
              else
                  out_grid.x = xlist;
              end
          else
              error('%s: ERROR: Input ylist should be a vector.',mfilename)
          end
          if isvector(ylist)
              if size(ylist,2) > 1
                  out_grid.y = ylist';
              else
                  out_grid.y = ylist;
              end
          else
              error('%s: ERROR: Input ylist should be a vector.',mfilename)
          end
          out_grid.Map = zlist;
          out_grid.nx = numel(out_grid.x);
          out_grid.ny = numel(out_grid.y);

          % Calculate the spacing and limits
          if isDefaultGridLimits || isempty(GridLimits)
              out_grid.GridLimits = grid_getLimits(xlist(:),ylist(:));
          else
              out_grid.GridLimits = GridLimits;
          end
          if isDefaultGridSpacing || isempty(GridSpacing)
              out_grid.GridSpacing = grid_getSpacing(xlist(:),ylist(:));
          else
              out_grid.GridSpacing = GridSpacing;
          end

      else
          error('%s: ERROR: Grid elements with inconsistent sizes.',mfilename);
      end

  elseif isgrid(xlist)
      fprintf(1,['$s: WARNING: This use of this function will be removed.\n'...
                 '             Use grid_set* functions instead.'],mfilename);
      % Simply copy grid. Remaining fields will be edited after this if.
      out_grid = xlist;
  elseif ischar(xlist) || iscell(xlist)
      % Read files specified in xlist
      out_grid = grid_read(xlist);
      % exit function
      return
  else
      error('%s: ERROR: invalid input xlist.',mfilename)
  end
  % -> Here we are sure a valid grid structure is defined with. From now on
  % it should be possible to use isgrid.
  %   out.grid . x
  %            . nx
  %            . y
  %            . ny
  %            . Map
  %            . GridLimits
  %            . GridSpacing
  %            . Units
  %            . Date
  %            . Description

  %% Now we add the remaining fields optional fields.
  % % Calculate and check the spacing and limits
  % if ~isDefaultGridLimits && any(GridLimits ~= out_grid.GridLimits)
  %     % different limits specified, verify
  %     out_grid.GridLimits = GridLimits;
  %     if ~isgrid(out_grid)
  %         limStr = sprintf('[%G,%G,%G,%G]',grid_get_limits(out_grid.x,out_grid.y));
  %         error('%s: ERROR: Incorrect limits. Should be a 4x1 vector at least spaning the range: %s',...
  %             mfilename,limStr);
  %     end
  %     hasNewLimits = true;
  % else
  %     hasNewLimits = false;
  % end
  %
  % if ~isDefaultGridSpacing && any(GridSpacing ~= out_grid.GridSpacing)
  %     % spacing specified, verify
  %     out_grid.GridSpacing = GridSpacing;
  %     if ~isgrid(out_grid)
  %         spStr = sprintf('[%G,%G]',grid_get_spacing(out_grid.x,out_grid.y));
  %         error('%s: ERROR: Incorrect spacing. Should be a 2x1 vector with maximum values: %s',...
  %             mfilename,spStr);
  %     end
  %     hasNewSpacing = true;
  % else
  %     hasNewSpacing = false;
  % end
  %
  % if hasNewLimits || hasNewSpacing
  %     % Expand grid domain.
  %     newx =  out_grid.GridLimits(1):out_grid.GridSpacing(1):out_grid.GridLimits(2);
  %     newy = (out_grid.GridLimits(3):out_grid.GridSpacing(2):out_grid.GridLimits(4))';
  %
  %     % map old coordinates into new coordinates
  %     [aux,aux,idxx] = intersect(out_grid.x,newx);
  %     [aux,aux,idxy] = intersect(out_grid.y,newy);
  %
  %     % New grid
  %     out_grid.x  = newx;
  %     out_grid.y  = newy;
  %     out_grid.nx = numel(out_grid.x);
  %     out_grid.ny = numel(out_grid.y);
  %
  %     % Expand mat with empty values
  %     newMap = nan(out_grid.ny,out_grid.nx);
  %     % Map old mat into new one
  %     newMap(idxy,idxx) = out_grid.Map;
  %     out_grid.Map = newMap;
  % end

  % Add Units to structure
  if ~isDefaultUnits
      % If Units were specified check and add them to the structure
      out_grid = grid_setUnits(out_grid,Units);
  end

  % Check the date
  if ~isDefaultDate && ~isempty(Date)
      % If Date was specified check and add them to the structure
      out_grid = grid_setDate(out_grid,Date);
  end

  % Add Description
  if ~isDefaultDescription
      % If Description was specified check and add them to the structure
      out_grid = grid_setDescription(out_grid,Description);
  end
end
function GridLimits = grid_getLimits(xlist,ylist)
% LIMITS = GRID_GETLIMITS(GRID) computes the limits of the grid. The
% limits are a 4x1 vectory with the maximum and minimum values of the x-
% and y- coordinates of the grid.
%
%   GRID_GETLIMITS(XLIST,YLIST) also accepts a list of points to compute
%   the limits from.

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>

  if nargin == 1 && isgrid(xlist)
      grid  = xlist;
      xlist = grid.x;
      ylist = grid.y;
  elseif nargin == 2 && isvector(xlist) && isvector(ylist)
      % continue
  elseif nargin == 0 || ( nargin == 2 && isempty(xlist) && isempty(ylist) )
      % possible empty inputs
      xlist = NaN;
      ylist = NaN;
  else
      error('%s: ERROR: Inputs must be either a valid grid or a list of x- and y- points.',mfilename)
  end

  GridLimits = [ min(xlist) max(xlist) min(ylist) max(ylist) ];
end
function GridSpacing = grid_getSpacing(xlist,ylist)
% SPACING = GRID_GETSPACING(GRID) computes the spacing of the input grid.
% The spacing is the 2x1 vectory with the minimum step between the values
% of the x- and y- coordinates of the grid.
%
%   GRID_GETSPACING(xlist,ylist) also accepts a list of points to compute
%   the limits from.

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>

  if nargin == 1 && isgrid(xlist)
      grid  = xlist;
      xlist = grid.x;
      ylist = grid.y;
  elseif nargin == 2 && isvector(xlist) && isvector(ylist)
      % continue
  elseif nargin == 0 || ( nargin == 2 && isempty(xlist) && isempty(ylist) )
      % possible empty inputs
      xlist = NaN(2,1);
      ylist = NaN(2,1);
  else
      error('%s: ERROR: Inputs must be either a valid grid or a list of x- and y- points.',mfilename)
  end

  GridSpacing = [ min(diff(sort(unique(xlist)))) min(diff(sort(unique(ylist)))) ];
end
function grid = grid_setDate(grid,Date)
% GRID = GRID_SETDATE(GRID,DATE) set the date field of the grid structure.
%
%   DATE should be a valid date in the format 'dd-mm-yyyy HH:MM:SS'
%   If DATE is unspecified, the Date is set to the current date.

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>
  DATEFMT = 'yyyy/mm/dd HH:MM:SS';

  if nargin == 0
      grid = DATEFMT;
  elseif isgrid(grid)
      if ~exist('Date','var') || isempty(Date)
          grid.Date = datestr(now,DATEFMT);
      else
          % Check the input grid format
          try
              datenum(Date,DATEFMT);
          catch exception
              if strcmp(exception.identifier,'MATLAB:datenum:ConvertDateString')
                  error('%s: ERROR: Date must be in format: %s',mfilename,DATEFMT);
              else
                  throw(exception)
              end
          end
          grid.Date = Date;
      end
  else
      error('%s: ERROR: Input must be a valid grid structure.',mfilename)
  end
end
function grid = grid_setDescription(grid,Description)
% GRID = GRID_SETDESCRIPTION(GRID,DESC) set the description of the grid
% structure.
%
%   DESC should be a string. If DESC is unspecified it is deleted.

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>
  if isgrid(grid)
      if ~exist('Description','var') || isempty(Description)
          %do nothing
      elseif ischar(Description)
          grid.Description = Description;
      else
          error('%s: ERROR: DESC must be a valid string.',mfilename)
      end
  else
      error('%s: ERROR: Input must be a valid grid structure.',mfilename)
  end
end
function grid = grid_setUnits(grid,Units)
% GRID = GRID_SETUNITS(GRID,UNITS) set the units field of the grid structures.
%
%   UNITS should be a three string cell specifying the units,
%       e.g.: {'px','px','m'}
%   If UNITS is unspecified, calling this function will erase existing units.

% Created by P.Inacio <P.M.G.Inacio@tudelft.nl>
  if isgrid(grid)
      % Check input units
      if ~exist('Units','var') || isempty(Units)
          grid.Units = '';
      elseif iscell(Units) && all(size(Units) == [1 3])
          grid.Units = Units;
      else
          error('%s: ERROR: Must provide valid Units to constructor. Should be a 1x3 cell array.',mfilename);
      end
  else
      error('%s: ERROR: Input must be a valid grid structure.',mfilename)
  end
end
function [h,cbh]=plot_grid(glon,glat,gvals,projection,origin_lon,origin_lat,MapLonLimit,...
            MapLatLimit,coastFlag,colorAxis,colorbarMode,h,visible,titleFlag)
% PLOT_GRID plots a grid defined in different formats.
%
%   PLOT_GRID(GRID) plots GRID, a structure created with grid_constructor.
%   PLOT_GRID(X,Y,MAP) plots a grid defined by latitude,longitude and
%   mat of values.
%
%   PLOT_GRID(...,PROJ) to chose the kind of projection to be used in the
%   plot. Available modes are:
%       - 'image','surf','contour' and all map projections supported by
%         MATLAB. Use command 'maps' for a list of them.
%
%   PLOT_GRID(...,ORIGINLON,ORIGINLAT,MAPLONLIMIT,MAPLATLIMIT) allows the
%   user to specify the origin and limits of the plotted map.
%
%   PLOT_GRID(...,COASTFLAG) to specify whether the coastline is plotted
%   along the data. Defualt is true;
%
%   PLOT_GRID(...,COLORAXIS,COLORBARMODE) also allow you to specify the
%   color axis. COLORBAR_MODE allows you to configure the colorbar
%   behaviour: no colorbar (0); symmetric colorbar (1); natural colorbar
%   (2). If the same values are used but negative a horizontal colorbar is
%   placed instead of the vertical one.
%
%   PLOT_GRID(...,H,VISIBLE,TITLEFLAG) to plot in the figure/subplot with
%   handle H, a flag to control visibility and a flag to control if some
%   statistics are displayed in the title of the plot.

% Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>
% List of changes:
%
% P.Inacio <P.M.G.Inacio@tudelft.nl>, adapted to only use the standard grid
%   structure defined by grid_constructor.
% P.I, added color axis and colorbar mode inputs.
% P.I, added figurehandle, visible and titlestring inputs.

% TODO: use input parser
% TODO: use coast flag input to set line thickness
% TODO: handle multiple grids
% TODO: handle invisible plots properly.

  %% Defaults
  line_width=1.5;
  line_color='k';
  isMap = false;

  %% Check inputs
  if ischar(glon)
      glon=grid_read(glon);
  end
  if isgrid(glon)
      % save input grid
      grid = glon;

      % Change variables contents to the right place. This can be avoided by
      % using input parser.
      if (nargin > 12); error('%s: ERROR: Too many input arguments.',mfilename); end
      if exist('h','var'); titleFlag = h;  end
      if exist('colorbarMode','var'); visible = colorbarMode;  end
      if exist('colorAxis','var'); h = colorAxis;  end
      if exist('coastFlag','var'); colorbarMode = coastFlag;  end
      if exist('MapLatLimit','var'); colorAxis  = MapLatLimit;end
      if exist('MapLonLimit','var'); coastFlag  = MapLonLimit;end
      if exist('origin_lat','var'); MapLatLimit = origin_lat; end
      if exist('origin_lon','var'); MapLonLimit = origin_lon; end
      if exist('projection','var'); origin_lat  = projection; end
      if exist('gvals','var'); origin_lon = gvals; end
      if exist('glat','var'); projection = glat; end

      % extract grid values
      % This does not make much sense, latitude is always [-90;90]
  %     glat =ang_fix_domain(grid.y,-90,90);
      glat =grid.y;
      glon =ang_fix_domain(grid.x,0,360);
      gvals=grid.Map;
  elseif ~exist('glat','var') || ~exist('gvals','var') || ~isnumeric(glon) || ~isnumeric(glat) || ~isnumeric(gvals)
      error([mfilename,': ERROR: Input grid must be either valid grid structure or a set of numeric x,y,z vals.']);
  end

  % Define default behaviour
  if ~exist('MapLonLimit','var') || isempty(MapLonLimit)
      MapLonLimit = [0 360];
  end
  if ~exist('MapLatLimit','var') || isempty(MapLatLimit)
      MapLatLimit = [-90 90];
  end
  if ~exist('origin_lat','var') || isempty(origin_lat)
      origin_lat = 0;
  end
  if ~exist('origin_lon','var') || isempty(origin_lon)
      origin_lon = 0;
  end
  if ~exist('projection','var') || isempty(projection)
      projection = 'robinson';
  end
  if ~exist('coastFlag','var') || isempty(coastFlag)
      coastFlag = true;
  end
  if ~exist('colorAxis','var')
      colorAxis = [];
  end
  if ~exist('colorbarMode','var')
      colorbarMode = [];
  end
  if ~exist('h','var') || isempty(h) || ~ishandle(h)
  %     h = plot_format('a4landscape',figure('Visible',visible));
      h = gcf;
  end
  if ~exist('visible','var') || isempty(visible)
      visible = get(h,'Visible');
  end
  if ~exist('titleFlag','var') || isempty(titleFlag)
      titleFlag = true;
  end

  %% Start Plotting

  if exist('grid','var')
      % if it is a grid structure, check the units and decide on the type of
      % plot accordingly.
      % 'deg' and 'rad' trigger map defaults. Others trigger image or surface.
      if all(strcmpi(grid.Units(1:2),'deg'))

          % Add longitude 360 if it is not defined in the dataset, and copy
          % values from longitude 0.
          if glon(end)-360 ~= glon(1)
              glon(end+1) = glon(1)+360;
              gvals(:,end+1) = gvals(:,1);
          end

          % P.Inacio, 09-2011, with grid structure it is possible to assume the
          %   shape of glat and glon instead of using the old code.
          % Create grids for plotting
          glat = repmat(glat,[1 size(gvals,2)]);
          glon = repmat(glon,[size(gvals,1) 1]);

      elseif all(strcmpi(grid.Units(1:2),'rad'))
          % Convert to degrees and use recursion for plotting
          plot_grid(grid_convert(grid,{'deg','deg',''}),projection,origin_lon,origin_lat,MapLonLimit,MapLatLimit,coastFlag,colorAxis,colorbarMode,h)
          return;
      elseif any(strcmp(grid.Units(1:2),''))
          fprintf(1,'%s: WARNING: Unspecified x- and y- grid units.\n',mfilename);
          % Handle generic units
          % TODO: apply these only if they are not specified by caller.
          projection = 'image';
          coastFlag  = false;
      else
          % Handle generic units
          % TODO: apply these only if they are not specified by caller.
          projection = 'image';
          coastFlag  = false;
      end
  else
      % It is not a grid, use provided settings
      [glon,glat,gvals]=grid_check_sizes(glon,glat,gvals,'mat');

      %patching longitude 180
      if ~any(glon(1,:)==180)
          glon_new=sort([glon(1,:),180]);
          [glon_new,glat_new]=meshgrid(glon_new,glat(:,1));
          [glon,glat,gvals]=griddata(glon,glat,gvals,glon_new,glat_new);
      end
      %patching longitude 360
      if glon(1,end) ~= 360
          glon(:,end+1)=360;
          glat(:,end+1)=glat(:,end);
          gvals(:,end+1)=gvals(:,1)';
      end

      glon=ang_fix_domain(glon,-180,180);
      %patching longitude -180
      if ~any(glon(1,:)==-180)
          patch_idx=find(glon(1,:)==180);
          %bug trap
          if isempty(patch_idx)
              error([mfilename,': ERROR: isempty(patch_idx), debug needed.'])
          end
          glon=[glon(:,1:patch_idx),-180*ones(size(glon,1),1),glon(:,patch_idx+1:end)];
          glat=[glat,glat(:,1)];
          gvals=[gvals(:,1:patch_idx),gvals(:,patch_idx),gvals(:,patch_idx+1:end)];
      end
  end


  % TODO: handle invisible plots properly.
  try
      % If it is a subplot, this will make it active
      subplot(h);
  catch me
      % If figure, this will make it active
      if strcmp(me.identifier,'MATLAB:subplot:SubplotIndexOutOfRange')
          figure(h); set(h,'visible',visible)
      end
  end


  %% Plotting
  switch projection
      case {'img','image'}
          %need to sort
          imagesc(grid.x,grid.y,grid.Map), hold on, axis xy
          %invert y-axis
          set(gca,'YDir','reverse');
          %ploting coastline
          if coastFlag
              load coast
              plot3(long,lat,ones(size(lat))*max(gvals(:))+1,line_color,'LineWidth',line_width)
          end
      case {'surf','surface'}
          surf(glon,glat,gvals), hold on
          shading interp, view([0 90]), axis tight
          %ploting coastline
          if coastFlag
              load coast
              plot3(long,lat,ones(size(lat))*max(gvals(:)),line_color,'LineWidth',line_width)
          end
      case 'contour'
          contourf(glon,glat,gvals), hold on
          %ploting coastline
          if coastFlag
              load coast
              plot(long,lat,line_color,'LineWidth',line_width)
          end
      case '3d'
          % Scale the radial axis
          sc = 0.25;
          R = ones(size(glon))+sc*gvals/max(abs(gvals(:)));
          X = R.*cosd(glon).*cosd(glat);
          Y = R.*sind(glon).*cosd(glat);
          Z = R.*sind(glat);
          hold on
          surf(X,Y,Z,gvals,'EdgeColor','none','FaceColor','interp','FaceLighting','gouraud')
          if coastFlag
              offset = 0.02;
              line_width=2;
              load coast
              % Find the nearest points in the grid
              idxLon = num_map(glon(1,:),long);
              idxLat = num_map(glat(:,1),lat);
              idx=sub2ind(size(R),idxLat,idxLon);
              X = (R(idx)+offset).*cosd(long).*cosd(lat);
              Y = (R(idx)+offset).*sind(long).*cosd(lat);
              Z = (R(idx)+offset).*sind(lat);
              plot3(X,Y,Z,'color',line_color,'LineWidth',line_width)
          end
          axis equal
          return
      otherwise
          % set map flag
          isMap = true;

          ax=worldmap(MapLatLimit,MapLonLimit);
          setm(ax,'MapProjection',projection);
          % P.Inacio, sometimes changing the projection overrides the
          % previously selected lat lon limits.
          setm(ax,'MapLatLimit',MapLatLimit);
          setm(ax,'MapLonLimit',MapLonLimit);
          setm(ax,'Frame','on');
          setm(ax,'Origin',[origin_lat origin_lon 0]);

          setm(ax,...
                  'MeridianLabel','off',...
                  'ParallelLabel','off')

          % NOTE: replacing NaN's with 0's. This might not always be desired.
          % The best is to set transparency on the NaN values, find out how
          % to do this.
          gvals(isnan(gvals)) = 0;

          % plot grid on map
          geoshow(ax,glat,glon,gvals,'DisplayType','texturemap')
          %ploting coastline
          if coastFlag
              load coast
              geoshow(ax,lat, long,'LineWidth',line_width,'Color',line_color)
          end
  end

  %colorbar domain: caxis must be called before <plot_colorbar_center> or
  %<plot_colorbar_opt>, otherwise the operation that those routines do is
  %no longer the desired one (the white color is no longer at zero, for
  %example).
  if ~isempty(colorAxis)
      caxis(colorAxis);
  end
  %nice colormap
  if isempty(colorbarMode)
      if any(gvals(:) > 0.001*mean(gvals(:))) && any(gvals(:) < 0.001*mean(gvals(:)))
          colorbarMode=-1;
      else
          colorbarMode=-2;
      end
  end
  switch abs(colorbarMode)
      case 0
          %no colorbar
      case 1
          % make colorbar symmetric if no color limits have been specified.
          if isempty(colorAxis)
              caxis([ -max(abs(caxis)) max(abs(caxis)) ]);
          end
          plot_colorbar_center(0,jet(256))
      case 2
          plot_colorbar_opt(jet(256))
      otherwise
          error([mfilename,': unknown colorbar_mode ',num2str(colorbarMode),'.'])
  end
  % P.Inacio - This is not useful here because one does not need to freeze
  % the colors when plotting a given grid. If any other function used this
  % function to plot a grid, then that function should freeze the colors.
  % % This freezes the colormap
  % plot_freezeColors
  switch sign(colorbarMode)
      case 0
          %no colorbar
      case 1
  %         %default vertical colobar
            cbh = colorbar;
  %         % P.Inacio <P.M.G.Inacio@tudelft.nl> - plot_freezeColors no
  %         %   longer works on colorbars.
  %         cbh = plot_cbfreeze(colorbar);
  %         % P.Inacio - workaround because labels are normally wrong after
  %         %   freezing the colorbar.
  %         set(cbh,'YTickLabelMode','auto')
      case -1
            cbh = colorbar('Location','SouthOutside');
  %         % P.Inacio <P.M.G.Inacio@tudelft.nl> - plot_freezeColors no
  %         %   longer works on colorbars.
  %         cbh = plot_cbfreeze(colorbar('Location','SouthOutside'));
  %         % P.Inacio - workaround because labels are normally wrong after
  %         %   freezing the colorbar.
  %         set(cbh,'XTickLabelMode','auto')
  end

  % Title
  if titleFlag
      non_zero_idx=(gvals~=0);
      title(['RMS=',num2str(num_rms(gvals(non_zero_idx))),'; max=',num2str(max(abs(gvals(:))))])
      if ~isempty(grid.Description)
          plot_appendtotitle(grid.Description);
      end
  end
  %strecht map
  if isMap; tightmap; end

  %handle output
  if nargout == 0
      clear h
  end
end
function out = ang_fix_domain(in,lower,upper,string)
  % ANG_FIX_DOMAIN(IN,LOWER,UPPER) is the same as <check_domain>, except that
  % the bounds are taken to be connected, much like angular quantities that
  % should change smoothly betweem 0 and 360 degrees.
  %
  %   ANG_FIX_DOMAIN(IN,LOWER,UPPER,STRING) to pass the name of the angle
  %   being checked, in which case a warning message is displayed.

  % Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>

  if ~exist('string','var')
      string='';
  end
  if ~isscalar(lower) || ~isscalar(upper)
      error([mfilename,': inputs <lower> and <upper> must be scalar.'])
  end

  while any(in(:)>upper) || any(in(:)<lower)
      idx = (in(:) > upper);
      if any(idx)
          in(idx) = in(idx) - (upper-lower);
          if ~isempty(string) && isscalar(in)
              disp([mfilename,': ',string,'=',num2str(in),' (was ',num2str(in),...
                  ') needed fixing because upper circular bound is ',num2str(upper),'.'])
          end
      end
      idx = (in(:) < lower);
      if any(idx)
          in(idx) = in(idx) + (upper-lower);
          if ~isempty(string) && isscalar(in)
              disp([mfilename,': ',string,'=',num2str(in),' (was ',num2str(in),...
                  ') needed fixing because lower circular bound is ',num2str(lower),'.'])
          end
      end
  end
  out = in;
end

%% Spherical Harminc analysis
function [c_out,s_out,msg]=mod_sh_ana(long,lat,grid,N)
% [C,S]=MOD_SH_ANA(LONG,LAT,GRID) is the low-level spherical harmonic
% analysis routine. It should be a self-contained script.
%
%   Input LONG is a horizontal vector [0:2*pi].
%   Input LAT is a vertical vector [-pi/2:pi/2].
%   Input GRID is a matrix with size [length(lat),length(long)]
%   MOD_SH_ANA(LONG,LAT,GRID,N) with the optional input N sets the maximum
%   degree of the spherical harmonic expansion.
%
%   Outputs C and S are the cosine and sine coefficients, organized in
%   lower triangle matrices, with constant degree in each row and constant
%   order in each column, sectorial coefficients are in the diagonals, zonal
%   coefficients are in the first column of matrix C (first column of S
%   matrix is zeros).
%
%   Due to matlab indexing, degree/order 0 is in position 1, i.e.
%   degree/order N is in index N+1.
%
%   All angular quantities are in radians.

% Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>

% List of changes:
%   P.Inacio (p.m.g.inacio@tudelft.nl), 09/2011, Added check that grid
%       values at 0 and 2*pi latitudes are the same and use only 0 in the
%       analysis.
%   P.Inacio, 10/2011, Added check that grid values at latitudes ??pi/2 must
%       all be the same. Also added a more robust check that grids are
%       regular taking into account a numerical error threshold.

% TODO: Function doesnt seem to care which long,lat values were provided.
%       It assumes the observations are equally distributed over the
%       sphere. This could be a problem if, for example, one provides a grid
%       where the latitude values are from -60 to 60 deg.

  % Defaults
  NUM_THRS = 1E-14; % Numerical errors threshold to check for regular grid.

  %check dimensions
  if size(long,1) ~= 1
      error([mfilename,': input <long> must be a horizontal vector.'])
  end
  if size(lat,2) ~= 1
      error([mfilename,': input <lat> must be a vertical vector.'])
  end
  %check regular grid
  if any(diff(lat) <= 0) || any(abs(diff(lat,2)) > max(abs(lat))*NUM_THRS)
      error([mfilename,': input <lat> must be monotonically increasing set of equally-spaced latitudes.'])
  end
  if any(diff(long) <= 0) || any(abs(diff(long,2)) > max(abs(long))*NUM_THRS)
      % if long starts in pi->2pi|0->pi there is a "false" discontinuity
      aux = ang_fix_pi(long);
      if any(diff(aux) <= 0) || any(abs(diff(aux,2)) > max(abs(aux))*NUM_THRS)
          error([mfilename,': input <long> must be monotonically increasing set of equally-spaced longitudes.'])
      end
  end
  % check domain
  %  if max(lat) > pi/2 || min(lat) < -pi/2
  if max(lat)-pi/2 > 1E-14 || min(lat)+pi/2 < -1E-14
      error([mfilename,': input <lat> does not seem to be in radians or outside legal domain [-pi/2,pi/2].'])
  end
  if max(long) > 2*pi || min(long) < 0
      error([mfilename,': input <long> does not seem to be in radians or outside legal domain [0,2*pi].'])
  end
  % check repeated data
  if long(end) == 2*pi && long(1) == 0
      if any( grid(:,1) ~= grid(:,end) )
          error([mfilename,': Invalid grid, different values for same meridian (0 and 2*pi).'])
      else
          % do not use the 360 longitude
          grid(:,end) = [];
          long(end) = [];
      end
  end
  % check pole singularities
  if lat(1) == pi/2 && any(diff(grid(1,:)) ~= 0)
      error([mfilename,': input <lat> has different values for the same point with lat=pi/2.'])
  end
  if lat(end) == -pi/2 && any(diff(grid(end,:)) ~= 0)
      error([mfilename,': input <lat> has different values for the same point with lat=-pi/2.'])
  end
  % check grid size
  if any(size(grid) ~= [length(lat),length(long)])
      error([mfilename,': size of input <grid> is not compatible with size of inputs <lat> or <long>.'])
  end

  % parameters
  I = length(lat);
  J = length(long);

  %optional inputs
  N_max = min([I-1,floor((J-1)/2)]);
  if ~exist('N','var') || isempty(N)
      N = N_max;
  elseif N > N_max
      error([mfilename,': <N> is too large. The maximum for this grid is ',num2str(N_max),'.'])
  end

  %need co-latitude
  lat=pi/2-lat;

  % determine fourier coefficients
  a=zeros(I,N+1);
  b=zeros(I,N+1);
  for i=1:I
      for m=0:N
          j=1:J;
          a(i,m+1) = sum(grid(i,j).*cos(m*long(j)))/J*(2-delta_kronecher(0,m));
          b(i,m+1) = sum(grid(i,j).*sin(m*long(j)))/J*2;
      end
  end

  % calculating Legendre polynomials, P{degree}
  P=legendre_degree(N,lat);

  % building design matrices
  A=cell(1,N+1);
  for m=0:N
      for n=m:N
          A{m+1}(:,n-m+1)=P{n+1}(m+1,:)';
  %         disp(['m=',num2str(m),' n=',num2str(n),' n-m+1=',num2str(n-m+1),...
  %               ' size(P{n+1})=',num2str(size(P{n+1})),...
  %               ' size(A{m+1})=',num2str(size(A{m+1}))])
      end
  end

  % estimate spherical harmonics
  c=cell(1,N+1);
  s=cell(1,N+1);
  for m=0:N
      Nm=A{m+1}'*A{m+1};
      P=A{m+1}'*a(:,m+1);
      c{m+1}=Nm\P;
      P=A{m+1}'*b(:,m+1);
      s{m+1}=Nm\P;
  end

  %putting into matricial form, lower triangle, same degree in lines, same order in
  %columns, coefficients (k,k) in diagonal
  c_out = zeros(N+1);
  s_out = zeros(N+1);
  for i=0:N
      c_out(i+1:end,i+1) = c{i+1};
      s_out(i+1:end,i+1) = s{i+1};
  end

  %init message
  msg='';
  
  %NaNs and infs are set to zero
  if any(~isfinite(c_out(:)))
      msg=['WARNING: found ',num2str(sum(~isfinite(c_out(:)))),' NaN or inf cosine coefficients'];
      c_out(~isfinite(c_out))=0;
  end
  if any(~isfinite(s_out(:)))
      msg=['WARNING: found ',num2str(sum(~isfinite(s_out(:)))),' NaN or inf sine coefficients'];
      s_out(~isfinite(s_out))=0;
  end
end
function out = legendre_degree(N,lat)
  if min(size(lat)) ~= 1
      error([mfilename,': input <lat> must be a vector.'])
  end
  %getting legendre coefficients, per degree
  out=cell(1,N+1);
  for n=0:N
     out{n+1}=legendre(n,cos(lat'),'norm')*2;
     out{n+1}(1,:)=out{n+1}(1,:)/sqrt(2);
  end
end
function out = delta_kronecher(i,j)
  out = double(i==j);
end

