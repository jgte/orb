classdef grid < simpletimeseries
  %static
  properties(Constant)
    %default value of some internal parameters
    default_list=struct(...
      'tolerance',1e-8,...
      'lon_units','deg',...
      'lat_units','deg'...
    );
    parameter_list=struct(...
      'tolerance', struct('default',grid.default_list.tolerance, 'validation',@(i) isnumeric(i) && isscalar(i)),...
      'lon_units', struct('default',grid.default_list.lon_units, 'validation',@(i) ischar(i)),...
      'lat_units', struct('default',grid.default_list.lat_units, 'validation',@(i) ischar(i))...
    );
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'lon_units','lat_units'};
  end
  %read only
  properties(SetAccess=private)
    lon
    lat
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    
  end
  properties(GetAccess=public,SetAccess=public)
    tolerance
    lon_units
    lat_units
  end
  %calculated only when asked for
  properties(Dependent)
    map
    vecmat
    list
    flatlist
    matrix
    latSpacing
    lonSpacing
    latLimit
    lonLimit
    interpolant
  end
  methods(Static)
    function out=parameters
      out=fieldnames(grid.parameter_list);
    end
    function out=default
      out=grid.default_list;
    end
    %% spacing
    function out=getSpacing(lon,lat)
      out = [ min(diff(unique(lon(:),'sorted'))) min(diff(unique(lat(:),'sorted'))) ];
    end
    function out=lat_spacing_valid(spacing)
      out= isnumeric(spacing) && isscalar(spacing) && mod(180,spacing)==0;
    end
    function out=lon_spacing_valid(spacing)
      out= isnumeric(spacing) && isscalar(spacing) && mod(360,spacing)==0;
    end
    %% lon and lat handlers
    function out=getLimits(lon,lat)
      out = [ min(lon(:)) max(lon(:)) min(lat(:)) max(lat(:)) ];
    end
    function out=lon_default(n)
      out=linspace(0,360,n);
    end
    function out=lat_default(n)
      out=transpose(linspace(-90,90,n));
    end
    %% vector/matrix representation
    function out=vecmat_valid(vecmat)
      try
        grid.vecmat_check(vecmat);
        out=true;
      catch
        out=false;
      end
    end
    function out=vecmat_check(vecmat)
      %easier names
      type='vectmat'; 
      %need structure
      if ~isstruct( vecmat);       error([mfilename,': invalid ',type,': not a structure.']); end
      %need to check lon/lat fields existence now, before determining if they are uniform
      if ~isfield(    vecmat,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end 
      if ~isfield(    vecmat,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(vecmat.lon) && ~isempty(vecmat.lat)
        % Calculate the spacing and limits
        GridLimits  = grid.getLimits( vecmat.lon(:),vecmat.lat(:));
        GridSpacing = grid.getSpacing(vecmat.lon(:),vecmat.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon',          GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',transpose(GridLimits(3):GridSpacing(2):GridLimits(4))...
        );
        % sanity
        if numel(out.lon) ~= numel(vecmat.lon) || any( (out.lon(:)-vecmat.lon(:)).^2>grid.default_list.tolerance^2 )
          error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.'])
        end
        if numel(out.lat) ~= numel(vecmat.lat) || any( (out.lat(:)-vecmat.lat(:)).^2>grid.default_list.tolerance^2 )
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
    function out=vecmat_size(vecmat)
      out=size(vecmat.map);
    end
    function vectmat=vecmat_init(t,map,lon,lat)
      vectmat=struct('t',t,'lon',lon,'lat',lat,'map',map);
    end
    %% list representation
    function out=list_valid(list)
      try
        grid.list_check(list);
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
        GridLimits  = grid.getLimits( list.lon(:),list.lat(:));
        GridSpacing = grid.getSpacing(list.lon(:),list.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon', GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',(GridLimits(3):GridSpacing(2):GridLimits(4))'...
        );
        % get indexes
        idx.lon=nan(1,numel(list.lon));
        for i=1:numel(out.lon)
          idx.lon( (out.lon(i)-list.lon).^2<grid.default_list.tolerance^2 )=i;
        end
        idx.lat=nan(1,numel(list.lat));
        for i=1:numel(out.lat)
          idx.lat( (out.lat(i)-list.lat).^2<grid.default_list.tolerance^2 )=i;
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
            size(     list.map,2); error([mfilename,': invalid ',type,': field ''lon'' has different elements than rows of ''map''.']);end %#ok<ALIGN>
      end
      %lat
      if ~isempty(    list.lat)
        if ~isrow(    list.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not a row vector.']);end
        if ~isnumeric(list.lat);   error([mfilename,': invalid ',type,': field ''lat'' is not numeric.']);end
        if ~isempty(  list.map) && ...
            numel(    list.lat) ~= ...
            size(     list.map,2); error([mfilename,': invalid ',type,': field ''lat'' has different sizes than rows of ''map''.']);end %#ok<ALIGN>
      end
      %t
      if ~isfield(   list,'t');   error([mfilename,': invalid ',type,': field ''t'' is missing.']);end
      if ~isvector(  list.t);     error([mfilename,': invalid ',type,': field ''t'' is not a vector.']);end
      if ~isdatetime(list.t  ) && ...
         ~isnumeric( list.t  );   error([mfilename,': invalid ',type,': field ''t'' is not datetime or numeric.']);end
      if ~isempty(   list.map) && ...
          numel(     list.t)   ~= ...
          size(      list.map,1); error([mfilename,': invalid ',type,': field ''t'' has different nr of elements than pages of ''map''.']);end %#ok<ALIGN>
    end
    function out=list_size(list)
      out = [ numel(unique(list.lon(:))) numel(unique(list.lat(:))) numel(list.t) ];
    end
    function list=list_init(t,map,lon,lat)
      list=struct('t',t,'lon',lon,'lat',lat,'map',map);
    end
    %% flatlist representation
    function out=flatlist_valid(flatlist)
      try
        grid.flatlist_check(flatlist);
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
        GridLimits  = grid.getLimits( flatlist.lon(:),flatlist.lat(:));
        GridSpacing = grid.getSpacing(flatlist.lon(:),flatlist.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon', GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',(GridLimits(3):GridSpacing(2):GridLimits(4))',...
            't', unique(flatlist.t,'sorted')...
        );
        % get indexes
        idx.lon=nan(1,numel(flatlist.lon));
        for i=1:numel(out.lon)
          idx.lon( (out.lon(i)-flatlist.lon).^2<grid.default_list.tolerance^2 )=i;
        end
        if any(isnan(idx.lon)); error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.']);end
        idx.lat=nan(1,numel(flatlist.lat));
        for i=1:numel(out.lat)
          idx.lat( (out.lat(i)-flatlist.lat).^2<grid.default_list.tolerance^2 )=i;
        end
        if any(isnan(idx.lat)); error([mfilename,': invalid ',type,': field ''lat'' is not an uniform domain.']);end
        idx.t=nan(1,numel(flatlist.t));
        fltdn=datenum(flatlist.t);
        for i=1:numel(out.t)
          idx.t( (datenum(out.t(i))-fltdn).^2<grid.default_list.tolerance^2 )=i;
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
    function out=flatlist_size(list)
      out = [ numel(unique(list.lon(:))) numel(unique(list.lat(:))) numel(unique(list.t)) ];
    end
    function list=flatlist_init(t,map,lon,lat)
      list=struct('t',t,'lon',lon,'lat',lat,'map',map);
    end
    %% matrix representation
    function out=matrix_valid(matrix)
      try
        grid.matrix_check(matrix);
        out=true;
      catch
        out=false;
      end
    end
    function out=matrix_check(matrix)
      %easier names
      type='matrix';
      %need structure
      if ~isstruct( matrix);       error([mfilename,': invalid ',type,': not a structure.']); end
      %need to check lon/lat fields existence now, before determining if they are uniform
      if ~isfield(    matrix,'lon'); error([mfilename,': invalid ',type,': field ''lon'' is missing.']);end 
      if ~isfield(    matrix,'lat'); error([mfilename,': invalid ',type,': field ''lat'' is missing.']);end
      if ~isempty(matrix.lon) && ~isempty(matrix.lat)
        % Calculate the spacing and limits
        GridLimits  = grid.getLimits( matrix.lon(:),matrix.lat(:));
        GridSpacing = grid.getSpacing(matrix.lon(:),matrix.lat(:));
        % build equal-spaced lon and lat
        out=struct(...
          'lon',          GridLimits(1):GridSpacing(1):GridLimits(2),...
          'lat',transpose(GridLimits(3):GridSpacing(2):GridLimits(4))...
        );
        %build meshed lon and lat
        [lon_m,lat_m]=meshgrid(out.lon,out.lat);
        % sanity
        if numel(lon_m) ~= numel(matrix.lon) || any( (lon_m(:)-matrix.lon(:)).^2>grid.default_list.tolerance^2 )
          error([mfilename,': invalid ',type,': field ''lon'' is not an uniform domain.'])
        end
        if numel(lat_m) ~= numel(matrix.lat) || any( (lat_m(:)-matrix.lat(:)).^2>grid.default_list.tolerance^2 )
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
    function out=matrix_size(matrix)
      out=size(matrix.map);
    end
    function matrix=matrix_init(t,map,lon,lat)
      matrix=struct('t',t,'lon',lon,'lat',lat,'map',map);
    end
    %% vecmat and list convertions
    function [vecmat,idx]=list2vecmat(list)
      % inputs are lists of points in rows, changing with times along the columns
      [vecmat,idx]=grid.list_check(list);
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
      e=grid.vecmat_check(vecmat);
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
      list=grid.vecmat2list(vecmat);
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
      [vecmat,idx]=grid.flatlist_check(flatlist);
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
      vecmat=grid.matrix_check(matrix);
      % copy the rest
      vecmat.map = matrix.map;
      vecmat.t   = matrix.t;
    end
    function matrix=vecmat2matrix(vecmat)
      e=grid.vecmat_check(vecmat);
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
            case 'list';     out=grid.vecmat2list(in);
            case 'flatlist'; out=grid.vecmat2flatlist(in);
            case 'matrix';   out=grid.vecmat2matrix(in);
          end
        case 'list'
          switch lower(to)
            case 'vecmat';   out=                     grid.list2vecmat(in);
            case 'flatlist'; out=grid.vecmat2flatlist(grid.list2vecmat(in));
            case 'matrix';   out=grid.vecmat2matrix(  grid.list2vecmat(in));
          end
        case 'flatlist'
          switch lower(to)
            case 'vecmat'; out=                   grid.flatlist2vecmat(in);
            case 'list';   out=grid.vecmat2list(  grid.flatlist2vecmat(in));
            case 'matrix'; out=grid.vecmat2matrix(grid.flatlist2vecmat(in));
          end
        case 'matrix'
          switch lower(to)
            case 'vecmat';   out=                     grid.matrix2vecmat(in);
            case 'flatlist'; out=grid.vecmat2flatlist(grid.matrix2vecmat(in));
            case 'list';     out=grid.vecmat2list(    grid.matrix2vecmat(in));
          end
      end
    end
    %data type validity
    function c=dtv(type,in)
      switch lower(type)
      case 'vecmat';   c=grid.vecmat_valid(in);
      case 'list';     c=grid.list_valid(in);
      case 'flatlist'; c=grid.flatlist_valid(in);
      case 'matrix';   c=grid.matrix_valid(in);
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    %data type check
    function dtcheck(type,in)
      switch lower(type)
      case 'vecmat';grid.vecmat_check(in);
      case 'list';  grid.list_check(in);
      case 'flatlist';  grid.flatlist_check(in);
      case 'matrix';grid.matrix_check(in);
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    %data type length
    function l=dtl(type,in)
      switch lower(type)
      case 'vecmat';l=grid.vecmat_length(in);
      case 'list';  l=grid.list_length(in);
      case 'flatlist';  l=grid.flatlist_length(in);
      case 'matrix';l=grid.mat_length(in);
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
      elseif isnumeric(map) && size(map,2) == numel(t)
        %do nothing, this will be handled later
      else
        error([mfilename,': error in sizes of inputs ''t'' and ''map''.'])
      end
      % build structure with inputs
      si=struct(...
        't',t,...
        'lon',lon,...
        'lat',lat,...
        'map',map...
      );
      % determine format of inputs
      if grid.vecmat_valid(si)
        fmt='vecmat';
      elseif grid.list_valid(si)
        fmt='list';
      elseif grid.flatlist_valid(si)
        fmt='flatlist';
      elseif grid.matrix_valid(si)
        fmt='matrix';
      else
        error([mfilename,': cannot handle the format of inputs'])
      end
      % convert to requested type
      s=grid.dtc(fmt,type,si);
    end
    %% constructors
    function obj=unit(n_lon,n_lat,varargin)
      p=inputParser;
      p.addRequired( 'n_lon',  @(i) isscalar(i) && isnumeric(i));
      p.addRequired( 'n_lat',  @(i) isscalar(i) && isnumeric(i));
      p.addParameter('scale',1,@(i) isscalar(i) && isnumeric(i));
      p.addParameter('bias', 0,@(i) isscalar(i) && isnumeric(i));
      p.parse(n_lon,n_lat,varargin{:});
      %initialize with unitary grid
      obj=grid(...
        datetime('now'),...
        ones(p.Results.n_lat,p.Results.n_lon)*p.Results.scale+p.Results.bias...
      );
    end
    % Creates a random model with mean 0 and std 1
    function obj=unit_randn(n_lon,n_lat,varargin)
      p=inputParser;
      p.addRequired( 'n_lon',  @(i) isscalar(i) && isnumeric(i));
      p.addRequired( 'n_lat',  @(i) isscalar(i) && isnumeric(i));
      p.addParameter('scale',1,@(i) isscalar(i) && isnumeric(i));
      p.addParameter('bias', 0,@(i) isscalar(i) && isnumeric(i));
      p.parse(n_lon,n_lat,varargin{:});
      %initialize with random grid
      obj=grid(...
        datetime('now'),...
        randn(p.Results.n_lat,p.Results.n_lon)*p.Results.scale+p.Results.bias...
      );
    end
    function obj=slanted(n_lon,n_lat,varargin)
      p=inputParser;
      p.addRequired( 'n_lon',  @(i) isscalar(i) && isnumeric(i));
      p.addRequired( 'n_lat',  @(i) isscalar(i) && isnumeric(i));
      p.addParameter('scale',1,@(i) isscalar(i) && isnumeric(i));
      p.addParameter('bias', 0,@(i) isscalar(i) && isnumeric(i));
      p.parse(n_lon,n_lat,varargin{:});
      %initialize grid with sum of lon and lat in degrees
      [lon_map,lat_map]=meshgrid(grid.lon_default(p.Results.n_lon),grid.lat_default(p.Results.n_lat));
      %initialize with grid
      obj=grid(...
        datetime('now'),...
        lon_map+lat_map*p.Results.scale+p.Results.bias...
      );
    end
%     function out=load(filename)
%       %TODO
%     end
    %% tests
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
        l=[3,5,2];
      end
      
      dt={'vecmat','matrix','list','flatlist'};
      a=struct(...
        'lon',grid.lon_default(l(2)),...
        'lat',grid.lat_default(l(1)),...
        'map',reshape(1:prod(l),l),...
        't',datetime(sort(datenum(round(rand(l(3),6)*60))),'ConvertFrom','datenum')...
      );
      dd={...
        a,...
        grid.vecmat2matrix(a),...
        grid.vecmat2list(a),...
        grid.vecmat2flatlist(a)...
      };
      for i=1:numel(dt)
        for j=1:numel(dt)
          out=grid.dtc(dt{i},dt{j},dd{i});
          if any(any(out.map~=dd{j}.map))
            error([mfilename,': failed data type conversion between ''',dt{i},''' and ''',dt{j},'''.'])
          end
        end
      end
      
     a=grid.slanted(l(1),l(2));
     a.spatial_resample(l(1)*3,l(2)*2).imagesc(a.t)
      
    end
  end
  methods
    %% constructor
    function obj=grid(t,map,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('t',  @(i) ischar(i) || isnumeric(i) || isdatetime(i)); %this can be char, double or datetime
      p.addRequired('map',@(i) isnumeric(map) || iscell(map));
      p.addParameter('lat',grid.lat_default(size(map,1)), @(i) isnumeric(i));
      p.addParameter('lon',grid.lon_default(size(map,2)), @(i) isnumeric(i));
      %declare parameters
      for j=1:numel(grid.parameters)
        %shorter names
        pn=grid.parameters{j};
        %declare parameters
        p.addParameter(pn,grid.parameter_list.(pn).default,grid.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(t,map,varargin{:});
      % retrieve strucutre with list
      sl=grid.dti(t,p.Results.map,p.Results.lon,p.Results.lat,'list');
      % call superclass
      obj=obj@simpletimeseries(p.Results.t,sl.map,varargin{:});
      % convert to vecmat (to get lat/lon domains)
      sv=grid.dtc('list','vecmat',sl);
      % save lon and lat
      obj.lon=sv.lon;
      obj.lat=sv.lat;
      % save parameters
      for i=1:numel(grid.parameters)
        %shorter names
        pn=grid.parameters{i};
        if ~isscalar(p.Results.(pn))
          %vectors are always lines (easier to handle strings)
          obj.(pn)=transpose(p.Results.(pn)(:));
        else
          obj.(pn)=p.Results.(pn);
        end
      end
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
      else
        error([mfilename,': cannot assign a map without either ''x'' or ''t''.'])
      end
      % If inputs p.Results.lon and p.Results.lat are empty, it's because
      % obj.lat and obj.lon are too. This test makes sense because
      % simpledata will call the assign method before the lon/lat fields
      % are assigned.
      if ~isempty(p.Results.lon) && ~isempty(p.Results.lat)
        %checking if lon/lat are in the argument list
        if any(strcmp(p.UsingDefaults,'lon')) && ...
           any(strcmp(p.UsingDefaults,'lat'))
          % the format may not be properly determined, so assume map is given in list-type
        else
          % lon/lat give in inputs, retrieve structure with list (format is determined in grid.dti)
          sl=grid.dti(t,map,p.Results.lon,p.Results.lat,'list');
          % convert to vecmat (to get lat/lon domains)
          sv=grid.dtc('list','vecmat',sl);
          % save lon and lat
          obj.lon=sv.lon;
          obj.lat=sv.lat;
          % propagate map to be assigned below
          map=sl.map;
        end
      else
        %sanity
        if ~ismatrix(map)
          error([mfilename,': if lon/lat are empty, can only deal with list-type maps.'])
        end
        if numel(t) ~= size(map,1)
          error([mfilename,': input ''t'' or ''x'' must have the same nr of row as ''map''.'])
        end
      end
      % add width-reseting flag, if needed
      if obj.width~=size(map,2)
        if size(map,2)==numel(obj.lat)*numel(obj.lon)
          varargin{end+1}='reset_width';
          varargin{end+1}=true;
        else
          error([mfilename,': tryng to assign a map that is not consistent with existing lon/lat domain. ',...
            'Update those first or pass them as arguments to this member.'])
        end
      end
      % pass it upstream
      obj=assign@simpletimeseries(obj,map,varargin{:});
    end
    function obj=copy_metadata(obj,obj_in)
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in);
      %propagate parameters of this object
      parameters=grid.parameters;
      for i=1:numel(parameters)
        if isprop(obj,parameters{i}) && isprop(obj_in,parameters{i})
          obj.(parameters{i})=obj_in.(parameters{i});
        end
      end
    end
    %% print
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      relevant_parameters={'lat','lat_units','lon','lon_units','tolerance'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print superclass
      print@simpletimeseries(obj,tab)
    end
    %% vecmat handling
    function sv=get.vecmat(obj)
      sv=grid.vecmat_init(obj.t,obj.map,obj.lon,obj.lat);
    end
    function obj=set.vecmat(obj,vecmat)
      grid.vecmat_check(vecmat);
      obj=obj.assign(vecmat.map,...
                 't',vecmat.t,...
               'lat',vecmat.lat,...
               'lon',vecmat.lon...
      );
    end
    %% list handling 
    function sl=get.list(obj)
      sl=grid.dtc(obj.vecmat,'vecmat','list');
    end
    function obj=set.list(obj,list)
      grid.list_check(list);
      obj=obj.assign(list.map,...
                 't',list.t,...
               'lat',list.lat,...
               'lon',list.lon...
      );
    end
    %% flatlist handling 
    function sl=get.flatlist(obj)
      sl=grid.dtc(obj.vecmat,'vecmat','flatlist');
    end
    function obj=set.flatlist(obj,flatlist)
      grid.flatlist_check(flatlist);
      obj=obj.assign(flatlist.map,...
                 't',flatlist.t,...
               'lat',flatlist.lat,...
               'lon',flatlist.lon...
      );
    end
    %% matrix handling
    function sm=get.matrix(obj)
      sm=grid.dtc('vecmat','matrix',obj.vecmat);
    end
    function obj=set.matrix(obj,matrix)
      grid.matrixt_check(matrix);
      obj=obj.assign(matrix.map,...
                 't',matrix.t,...
               'lat',matrix.lat,...
               'lon',matrix.lon...
      );
    end
    %% map handling
    function out=get.map(obj)
      %convert lon/lat domain (internally in vecmat) to list
      sl=grid.vecmat2list(grid.vecmat_init(obj.t,[],obj.lon,obj.lat));
      %convert map (internally in list) to vecmat
      sv=grid.list2vecmat(grid.vecmat_init(obj.t,obj.y,sl.lon,sl.lat));
      %output
      out=sv.map;
    end
    function obj=set.map(obj,map)
      obj.vecmat=grid.vecmat_init(obj.t,map,obj.lon,obj.lat);
    end
    %% lat handling
    function obj=lat_set(obj,lat_now)
      %trivial call
      if numel(lat_now)==numel(obj.lat) && all(lat_now==obj.lat)
        return
      end
      obj=obj.spatial_interp(obj.lon,lat_now);
    end
    function out=lat_get(obj,type)
      switch lower(type)
      case 'list';  out=obj.list.lat;
      case 'vecmat';out=obj.lat;
      case 'matrix';out=obj.matrix.lat;
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    function out=get.latSpacing(obj)
      spacing=grid.getSpacing(obj.lon,obj.lat);
      out=spacing(2);
    end
    function obj=set.latSpacing(obj,spacing)
      obj=obj.lat_set(grid.lat_default(180/spacing));
    end
    %% lon handling
    function obj=lon_set(obj,lon_now)
      %trivial call
      if numel(lon_now)==numel(obj.lon) && all(lon_now==obj.lon)
        return
      end
      obj=obj.spatial_interp(lon_now,obj.lat);
    end
    function out=lon_get(obj,type)
      switch lower(type)
      case 'list';  out=obj.list.lon;
      case 'vecmat';out=obj.lon;
      case 'matrix';out=obj.matrix.lon;
      otherwise
        error([mfilename,': unknown data type ''',type,'''.'])
      end
    end
    function out=get.lonSpacing(obj)
      spacing=grid.getSpacing(obj.lon,obj.lat);
      out=spacing(1);
    end
    function obj=set.lonSpacing(obj,spacing)
      obj=obj.lon_set(grid.lon_default(300/spacing));
    end
    %% interpolant handling
    function out=get.interpolant(obj)
      [lat_meshed,lon_meshed,t_meshed]=ndgrid(obj.lat,obj.lon,datenum(obj.t));
      if obj.length==1
        out=griddedInterpolant(lat_meshed,lon_meshed,obj.map,'linear','none');
      else
        out=griddedInterpolant(lat_meshed,lon_meshed,t_meshed,obj.map,'linear','none');
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
      %get the interpolant
      I=obj.interpolant;
      %build new meshed domain and interpolate
      if obj.length==1
        [lat_meshed,lon_meshed]=ndgrid(lat_new,lon_new);
        map_new=I(lat_meshed,lon_meshed);
      else
        [lat_meshed,lon_meshed,t_meshed]=ndgrid(lat_new,lon_new,datenum(obj.t));
        map_new=I(lon_meshed,lat_meshed,t_meshed);
      end
      %save results 
      obj.lat=lat_new;
      obj.lon=lon_new;
      obj.map=map_new;
    end
    function obj=spatial_resample(obj,lon_n,lat_n)
      obj=obj.spatial_interp(grid.lon_default(lon_n),grid.lat_default(lat_n));
    end
    function obj=center_resample(obj)
      lon_new=mean([obj.lon(1:end-1);obj.lon(2:end)]);
      lat_new=mean([obj.lat(1:end-1),obj.lat(2:end)],2);
      obj=obj.spatial_interp(lon_new,lat_new);
    end
    %% multiple operands
    function compatible(obj1,obj2)
      %This method checks if the objectives are referring to the same
      %type of data, i.e. the data length is not important.
      parameters=grid.compatible_parameter_list;
      for i=1:numel(parameters)
        % if a parameter is empty, no need to check it
        if ( iscell(obj1.(parameters{i})) && isempty([obj1.(parameters{i}){:}]) ) || ...
           ( ischar(obj1.(parameters{i})) && isempty( obj1.(parameters{i})    ) ) || ...
           ( iscell(obj2.(parameters{i})) && isempty([obj2.(parameters{i}){:}]) ) || ...
           ( ischar(obj2.(parameters{i})) && isempty( obj2.(parameters{i})    ) )
          continue
        end
        if ~isequal(obj1.(parameters{i}),obj2.(parameters{i}))
          error([mfilename,': discrepancy in parameter ',parameters{i},'.'])
        end 
      end
      %call mother routine
      compatible@simpledata(obj1,obj2);
    end
    function [obj1,obj2]=consolidate(obj1,obj2)
      %match the compatible parameters 
      parameters=grid.compatible_parameter_list;
      for i=1:numel(parameters)
        p=(parameters{i});
        if ~isequal(obj1.(p),obj2.(p))
          obj2=obj2.scale(p,obj1.(p));
        end 
      end
      %extend the time-domain of both objects to be in agreement with the each other.
      [obj1,obj2]=consolidate@simpletimeseries(obj1,obj2);
    end
    %% plot functions
    function out=imagesc(obj,t,varargin)
      p=inputParser;
      p.addParameter('show_coast',   false,@(i) islogical(i));
      p.addParameter('show_colorbar',false,@(i) islogical(i));
      p.addParameter('bias', 0,@(i) isscalar(i) && isnumeric(i));
      p.parse(varargin{:});
      %interpolate at the requested time and resample to center of grid
      obj_interp=obj.interp(t).center_resample;
      %build image
      out.axis_handle=imagesc(obj_interp.lon,obj_interp.lat,obj_interp.map);hold on
      axis xy
      %invert y-axis
      set(gca,'YDir','reverse');
      %labels
      xlabel(['lon [',obj.lon_units,']'])
      ylabel(['lat [',obj.lat_units,']'])
      %ploting coastline
      if p.Results.show_coast
        out.coast_handle=obj.coast;
      end
      %add colorbar
      if p.Results.show_colorbar
        out.colorbar_handle=colorbar;
      end
      
    end
    function h=coast(obj,varargin)
      p=inputParser;
      p.addParameter('line_color','k',@(i) ischaracter(i));
      p.addParameter('line_width',1,  @(i) isscalar(i) && isnumeric(i));
      p.parse(varargin{:});
      coast = load('coast');
      h=plot3(...
        coast.long,...
        coast.lat,...
        ones(size(coast.lat))*max(obj.y(:))+1,...
        p.Results.line_color,...
        'LineWidth',p.Results.line_width...
      );
    end
%     function out=plot(obj,method,varargin)
%       % Parse inputs
%       p=inputParser;
%       p.KeepUnmatched=true;
%       % optional arguments
%       p.addParameter('something',false,@(i)islogical(i));
%       % parse it
%       p.parse(varargin{:});
% 
%     end
  end
end

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