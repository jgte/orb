 classdef simpledata
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list={...
      'x_tol',        1e-6,     @num.isscalar;...
      'peeklength',   10,       @num.isscalar;...
      'peekwidth',    10,       @num.isscalar;...
      'labels',       {''},     @iscell;...
      'units',        {''},     @iscellstr;...
      'x_units',      '',       @ischar;...
      'descriptor',   '',       @ischar;...
      'plot_zeros',   true,     @(i) islogical(i) && isscalar(i);...
      'invalid',      999999999,@num.isscalar;...
      'outlier_sigma',4,        @num.isscalar;...
      'cdate',datetime('now'),  @isdatetime;...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'units','x_units'};
  end
  %read only
  properties(SetAccess=private)
    length
    width
    x
    y
    mask
  end
  %These parameters should not modify the data in any way; they should
  %only describe the data or the input/output format of it.
  %NOTE: edit this if you add a new parameter (if read/write)
  properties(GetAccess=public,SetAccess=public)
    labels
    units
    x_units
    descriptor
    peeklength
    peekwidth
    x_tol
    plot_zeros
    invalid
    outlier_sigma
    cdate
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(simpledata.parameter_list); end
      out=v.picker(varargin{:});
    end
    function out=valid_x(x)
      out=(isnumeric(x) && ~isempty(x) && isvector(x)) || simpletimeseries.valid_t(x);
    end
    function [x,y,mask]=monotonic(x,y,mask,mode,debug)
      if ~exist('mode','var') || isempty(mode)
        mode='check';
      end
      if ~exist('debug','var') || isempty(debug)
        debug=true;
      end
      if ~issorted(x)
        [x,idx]=sort(x);
        y=y(idx,:);
        if ~isempty(mask)
          mask=mask(idx);
        end
      end
      %data should be soft monotonic now
      assert(~any(diff(x)<0),'sorting failed')
      %get index of duplicates (last one is not duplicate by definition)
      bad_idx=[diff(x)==0;false];
      if any(bad_idx)
        switch lower(mode)
        case 'check'
          %get index+1 of duplicates
          bad_idx1=circshift(bad_idx,1,1);
          %make sure the values are the same: check along columns
          for i=1:size(y,2)
            assert(~any(y(bad_idx,i)~=y(bad_idx1,i)),['cannot monotonize the data because column ',...
              num2str(i),' has different data in common epochs.'])
          end
        case {'remove','clear','clean'}
          disp([' WARNING: blindly removing ',num2str(sum(bad_idx)),' epoch(s) to make abcissae monotonic.'])
          %this is done below and is valid for all modes
        case {'average','mean'}
          disp([' WARNING: need to average ',num2str(sum(bad_idx)),' epoch(s) to make abcissae monotonic.'])
          %convert logical array to indexes
          bi=find(bad_idx);
          %loop over all indexes
          for i=1:numel(bi)
            %find idx in the data that have the same x-values
            idx=find(x==x(bi(i)));
            %get the index of the data entry where the average is going to be save into (c of 'collapsed')
            cidx=setdiff(idx,bi);
            %compute the collapsed y values
            cy=mean(y(idx,:));
            %noisy debug
            if debug
              tab=20;
              disp('DEBUG: from the following bad_idx:')
              disp(str.tablify(tab,bi(i),x(bi(i)),y(bi(i),:)))
              disp(['DEBUG: from bad_idx(',num2str(i),')=',num2str(bi(i)),', averaged the following data:'])
              disp(str.tablify(tab,'idx','x(idx)','y(idx)'))
              for j=1:numel(idx)
                disp(str.tablify(tab,idx(j),x(idx(j)),y(idx(j),:)))
              end
              disp('DEBUG: into the following single epoch:')
              disp(str.tablify(tab,cidx,x(cidx),cy))
            end
            %save the averaged y value
            y(cidx,:)=cy;
          end
        otherwise
          error(['unknown mode ''',mode,'''.'])
        end
        %remove duplicate epochs
        x=x(~bad_idx);
        y=y(~bad_idx,:);
        if ~isempty(mask)
          mask=mask(~bad_idx);
        end
      end
    end
    function out=valid_y(y)
      out=isnumeric(y) && ~isempty(y); % && size(y,1) > size(y,2);
    end
    function out=valid_mask(mask)
      out=islogical(mask) && ~isempty(mask) && isvector(mask);
    end
    function out=isx(mode,x1,x2,tol)
      if ~exist('tol','var') || isempty(tol)
        tol=simpledata.parameters('x_tol');
      end
      switch mode
      case {'=','==','equal'}
        if numel(x1)==numel(x2)
         out=(x1(:)-x2(:)).^2<tol.^2;
        elseif isscalar(x1)
          out=(x1-x2(:)).^2<tol.^2;
        elseif isscalar(x2)
          out=(x1(:)-x2).^2<tol.^2;
        else
          out=false;
        end
        return
      case {'<','less','smaller'}
        out=x1<x2;
        out(simpledata.isx('==',x1,x2,tol))=false;
      case {'<=','lessorequal'}
        out=x1<x2;
        out(simpledata.isx('==',x1,x2,tol))=true;
      case {'>','more','larger'}
        out=x1>x2;
        out(simpledata.isx('==',x1,x2,tol))=false;
      case {'>=','moreorequal','largerorequal'}
        out=x1>x2;
        out(simpledata.isx('==',x1,x2,tol))=true;
      otherwise
        error(['unknown mode ''',mode,'''.'])
      end
    end
    function out=transmute(in)
      if isa(in,'simpledata')
        %trivial call
        out=in;
      else
        %transmute into this object
        if isprop(in,'t')
          out=simpledata(simpletimeseries.time2num(in.t),in.y,in.varargin{:});
        elseif isprop(in,'x')
          out=simpledata(in.x,in.y,in.varargin{:});
        else
          error('Cannot find ''t'' or ''x''. Cannot continue.')
        end
      end
    end
    function [x,varargout] = union(x1,x2,tol)
      %If x1 and x2 are (for example) time domains, with each entry being a time
      %value, x is the vector of entries in either <x1> or <x2>.
      %
      %<indx1> and <indx2> have the indexes that relate <x> to <x1> and <x2>, so
      %that:
      %
      %x(ind1) = x1;
      %x(ind2) = x2;
      %
      %All outputs are column vectors.
      %
      %This is an OR operation, x has elements of any x1 or x2.
      %
      %A practical application for this would be to reconstruct a signal that is
      %measured with two sensors, each in its own time domain.
      %
      % t1 = [1 2 3 4 5 6];
      % s1 = sin(t1);
      % t2 = [0 1.5 3 4.5 6 7.5];
      % s2 = sin(t2);
      % [t,i1,i2] = union(t1,t2);
      % s(i1) = s1;
      % s(i2) = s2;
      % plot(t1,s1,'bo',t2,s2,'ro',t,s)

      %trivial call
      if isempty(x1) || isempty(x2)
        x = [];
        varargout={[],[]};
        return
      end
      if ~exist('tol','var') || isempty(tol)
        tol=1e-6;
      end
      %sanity
      if ~isvector(x1) || ~isvector(x2)
        error('inputs <x1> and <x2> must be 1D vectors.');
      end
      if (isnumeric(x1) && any(isnan(x1))) || ...
         (isnumeric(x2) && any(isnan(x2)))
        error('cannot handle NaNs in the input arguments. Clean them first.')
      end
      %shortcuts
      if numel(x1)==numel(x2) && all(x1(:)==x2(:))
        x=x1;
        varargout={true(size(x1)),true(size(x2))};
        return
      end
      %get first output argument (need to use tolerance to captures time system conversion errors)
      xt=[x1(:);x2(:)];
      tol=tol/max(abs(xt));
      x=uniquetol(xt,tol);
      %maybe we're done here
      if nargout==1
        return
      end
      %get second argument
      varargout(1)={ismembertol(x,x1,tol)};
      %maybe we're done here
      if nargout==2
        return
      end
      %get third argument
      varargout(2)={ismembertol(x,x2,tol)};
    end
    function [out,idx]=rm_outliers(in,varargin)
      % NOTE: greatly simplified the version of this routine described below:
      %
      % [OUT,IDX,STD_ERR]=NUM_RM_OUTLIERS(IN,NSIGMA,SMOOTH_METHOD,SPAN_FACTOR,PLOTFLAG)
      % is an outlier detection routine. It looks at the points in IN and flags
      % the ones that deviate from the reference by NSGIMA*std as outliers.
      % OUT is equal to IN, with NaN's replacing the flagged outliers. IDX is
      % a logical mask which contains indexes of IN flagged as outliers.
      %
      %   NUM_RM_OUTLIERS(IN) removes the outliers in input IN using default
      %   settings.
      %
      %   NUM_RM_OUTLIERS(IN,NSIGMA) additionally specifies how many
      %   times bigger than the STD are the outliers. Default value is 4.
      %
      %   NUM_RM_OUTLIERS(IN,NSIGMA,METHOD,SPAN_FACTOR) can be used to specify
      %   different methods of outlier detection. In general an outlier can be
      %   seen as a data-point which deviates from a certain reference by a
      %   statistical limit. The reference taken depends greatly on the type of
      %   input signal IN. If IN is a random zero-mean process than outliers
      %   can be computed by the distance to the x-axis. However if IN is a noisy
      %   sinusoidal signal the above reference (x-axis) result in cropping
      %   the wave at min and max regions. For that reason different methods
      %   are available to defined the reference points,
      %   The following methods are available,
      %   'zero'    - Outliers referenced to the x-axis.
      %   'mean'    - Outliers referenced to the mean of the signal [default]
      %   'detrend' - Outliers referenced to the best fitting line to IN
      %   'polyfit' - Outliers referenced to the best fitting polynom of degree
      %               POLYDEGREE. This method requires the additional input:
      %       NUM_RM_OUTLIERS(IN,NSIGMA,METHOD,POLYDEGREE)
      %
      %   Additional methods are available throught the internal use of the
      %   'smooth' function. There are 'moving','lowess','loess','sgolay',
      %   'rlowess','rloess'. For all these methods an additional input SPAN is
      %   required,
      %       [OUT,IDX]=NUM_RM_OUTLIERS(IN,NSIGMA,METHOD,SPAN) which specifies
      %   the length SPAN of the smoothing window relative to the lenght of
      %   the input. The default span is 0.1. Negative values of SPAN are also
      %   possible and they are implemented to allow the user to set a SPAN
      %   which is independent of the signal length.
      %
      %   NUM_RM_OUTLIERS(...,PLOTFLAG) with PLOTFLAG set to 'true' or to a
      %   figure handle will plot the flagged outliers and the reference signal.

      % Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>
      % List of Changes:
      %   Pedro Inacio <p.m.g.inacio@tudelft.nl>, reorganized outlier detection
      %       routine.

      % Handle Inputs
      p=machinery.inputParser;
      p.addRequired( 'in',                       @isnumeric);
      p.addParameter('outlier_sigma',...
         simpledata.parameters('outlier_sigma'), @num.isscalar);
      p.addParameter('outlier_value', nan,       @num.isscalar);
      % parse it
      p.parse(in,varargin{:});
      %assume there are no outliers
      out=in;
      idx=false(size(in));
      %ignore non-positive outlier_sigmas
      if p.Results.outlier_sigma>0 && numel(in)>1
        % ignore NaNs
        nanidx = isnan(in(:));
        % make aux containers with the same size as input
        err=zeros(size(in));
        % compute error
        err(~nanidx) = in(~nanidx)-mean(in(~nanidx));
        std_err = std(err(~nanidx));
        %trivial call
        if std_err>0
          %outlier indexes
          idx = abs(err)>p.Results.outlier_sigma*std_err;
          %flag outliers
          out(idx(:))=p.Results.outlier_value;
        end
      end
    end
    function obj_list=op_multiple(op,obj_list,msg)
      %vector mode
      if iscellstr(op)
        if ~exist('msg','var');msg='';end
        for i=1:numel(op)
          obj_list=simpledata.op_multiple(op{i},obj_list,msg);
        end
        return
      end
      %sanity
      if ~iscell(obj_list)
        error(['need input ''obj_list'' to be a cell array, not a ',class(obj_list),'.'])
      end
      % save timing info
      if ~exist('msg','var') || isempty(msg)
        s.msg=[op,' op'];
      else
        s.msg=[op,' op for ',msg];
      end
      s.n=2*(numel(obj_list)-1); c=0;
      %forward op
      for i=1:numel(obj_list)-1
        [obj_list{i},obj_list{i+1}]=obj_list{i}.(op)(obj_list{i+1});
        c=c+1;s=time.progress(s,c);
      end
      %backward op
      for i=numel(obj_list):-1:2
        [obj_list{i},obj_list{i-1}]=obj_list{i}.(op)(obj_list{i-1});
        c=c+1;s=time.progress(s,c);
      end
    end
    function out=isequal_multiple(obj_list,columns,msg)
      %sanity
      if ~iscell(obj_list)
        error(['need input ''obj_list'' to be a cell array, not a ',class(obj_list),'.'])
      end
      % save timing info
      if ~exist('msg','var') || isempty(msg)
        s.msg='isequal op';
      else
        s.msg=['isequal op for ',msg];
      end
      s.n=numel(obj_list)-1; c=0;
      %declare type of output argument
      out=cell(numel(obj_list)-1,1);
      %loop over obj list
      for i=1:numel(obj_list)-1
        out{i}=obj_list{i}.isequal(obj_list{i+1},columns);
        c=c+1;s=time.progress(s,c);
      end
    end
    %% vector operators
    function out=vmean(obj_list,varargin)
      p=machinery.inputParser;
      p.addRequired( 'obj_list', @iscell);
      p.addParameter('weights',ones(size(obj_list))/numel(obj_list),@(i) isnumeric(i) && all(size(obj_list)==size(i)))
      p.parse(obj_list,varargin{:})
      out=obj_list{1}.*p.Results.weights(1);
      for i=2:numel(obj_list)
        out=out+obj_list{i}.*p.Results.weights(i);
      end
    end
    function out=vtimes(obj_list1,obj_list2)
      p=machinery.inputParser;
      p.addRequired( 'obj_list1', @iscell);
      p.addRequired( 'obj_list2', @(i) iscell(i) && all(size(obj_list1)==size(obj_list2)) )
      p.parse(obj_list1,obj_list2)
      out=cell(size(obj_list1));
      for i=1:numel(obj_list1)
        out{i}=obj_list1{i}.*obj_list2{i};
      end
    end
    function out=vscale(obj_list,scales)
      p=machinery.inputParser;
      p.addRequired( 'obj_list', @iscell);
      p.addRequired( 'scales',   @(i) isnumeric(i) && all(numel(obj_list)==numel(scales)) )
      p.parse(obj_list,scales)
      out=cell(size(obj_list));
      for i=1:numel(obj_list)
        out{i}=obj_list{i}.*scales(i);
      end
    end
    function out=vsqrt(obj_list)
      p=machinery.inputParser;
      p.addRequired( 'obj_list', @iscell);
      p.parse(obj_list)
      out=cell(size(obj_list));
      for i=1:numel(obj_list)
        out{i}=obj_list{i}.sqrt;
      end
    end
    %% constructors
    function out=one(x,width,varargin)
      out=simpledata(x(:),ones(numel(x),width),varargin{:});
    end
    function out=zero(x,width,varargin)
      out=simpledata(x(:),zeros(numel(x),width),varargin{:});
    end
    function out=randn(x,width,varargin)
      out=simpledata(x(:),randn(numel(x),width),varargin{:});
    end
    function out=sinusoidal(x,w,varargin)
      y_now=cell2mat(arrayfun(@(i) sin(x(:)*i),w,'UniformOutput',false));
      out=simpledata(x(:),y_now,varargin{:});
    end
    %% aux routines
    function out=parametric_component_fieldnames(varargin)
      %converts {'polynomial',[1 1 1],'sinusoidal',[n1, n2, n3]} into {p0,p1,p2,s1,s2,s3,c1,c2,c3}
      p=machinery.inputParser;
      p.addParameter('polynomial',[], @(i) isnumeric(i) || isempty(i));
      p.addParameter('sinusoidal',[], @(i) isnumeric(i) || isduration(i) || isempty(i));
      p.parse(varargin{:});
      out=cell(1,numel(p.Results.polynomial)+2*numel(p.Results.sinusoidal));
      for i=1:numel(p.Results.polynomial)
        out{i}=['p',num2str(i-1)];
      end
      for i=1:numel(p.Results.sinusoidal)
        out{numel(p.Results.polynomial)+i}=['s',num2str(i)];
      end
      for i=1:numel(p.Results.sinusoidal)
        out{numel(p.Results.polynomial)+numel(p.Results.sinusoidal)+i}=['c',num2str(i)];
      end
    end
    %% general test for the current object
    function out=test_parameters(field,varargin)
      %basic parameters
      switch field
      case 'l'; out=100; return
      case 'w'; out=4;   return
      end
      %optional parameters
      switch numel(varargin)
      case 0
        l=simpledata.test_parameters('l');
        w=simpledata.test_parameters('w');
      case 1
        l=varargin{1};
        w=simpledata.test_parameters('w');
      case 2
        l=varargin{1};
        w=varargin{2};
      end
      %more parameters
      switch field
      case 'args'
        out={...
          'units',   strcat(cellstr(repmat('y_unit-',w,1)),cellstr(num2str((1:w)'))),...
          'labels',    strcat(cellstr(repmat('label-', w,1)),cellstr(num2str((1:w)'))),...
          'x_unit',    'x_unit',...
          'descriptor','original'...
        };
      case 'x';             out=transpose(0:l-1);
      case 'T';             out=l./[3 5];
      case 'y_randn_scale'; out=0.2;
      case 'y_poly_scale';  out=[1 1 1]./[1 l l^2];
      case 'y_sin_scale';   out=[0.8 0.5];
      case 'y_cos_scale';   out=[0.2 1.2];
      case 'y_randn'
        out=simpledata.test_parameters('y_randn_scale')*randn(l,w);
      case 'y_poly'
        out=pardecomp(...
          simpledata.test_parameters('x',l),...
          [],...
          'np',numel(simpledata.test_parameters('y_poly_scale')),...
          'T',[],...
          'x',simpledata.test_parameters('y_poly_scale',l)'...
        ).y_sum*ones(1,w);
        % % The code below is the same as the code above
        % c=simpledata.test_parameters('y_poly_scale',l);
        % x1=simpledata.test_parameters('x',l);
        % out=zeros(l,w);
        % for i=1:numel(c)
        %   out=out+c(i)*(x1).^(i-1)*ones(1,w);
        % end
      case 'y_sin_T'
        out=pardecomp(...
          simpledata.test_parameters('x',l),...
          [],...
          'np',0,...
          'T',simpledata.test_parameters('T',l),...
          'x',[simpledata.test_parameters('y_sin_scale'),zeros(size(simpledata.test_parameters('y_cos_scale')))]...
        ).y_sum*ones(1,w);
        % % The code below is the same as the code above
        % c=simpledata.test_parameters('y_sin_scale');
        % T=simpledata.test_parameters('T',l);
        % x1=simpledata.test_parameters('x',l);
        % out=zeros(l,w);
        % for i=1:numel(T)
        %   out=out+c(i)*sin( ...
        %     2*pi/T(i)*x1*ones(1,w)...
        %   );
        % end
      case 'y_cos_T'
        out=pardecomp(...
          simpledata.test_parameters('x',l),...
          [],...
          'np',0,...
          'T',simpledata.test_parameters('T',l),...
          'x',[zeros(size(simpledata.test_parameters('y_sin_scale'))),simpledata.test_parameters('y_cos_scale')]...
        ).y_sum*ones(1,w);
        % % The code below is the same as the code above
        % c=simpledata.test_parameters('y_cos_scale');
        % T=simpledata.test_parameters('T',l);
        % x1=simpledata.test_parameters('x',l);
        % out=zeros(l,w);
        % for i=1:numel(T)
        %   out=out+c(i)*cos( ...
        %     2*pi/T(i)*x1*ones(1,w)...
        %   );
        % end
      case 'y_sin'
        c=simpledata.test_parameters('y_sin_scale');
        T=simpledata.test_parameters('T',l);
        x1=simpledata.test_parameters('x',l);
        out=zeros(l,w);
        for i=1:numel(T)
          out=out+c(i)*sin(2*pi/T(i)*x1*randn(1,w) + ones(l,1)*randn(1,w));
        end
      case 'y_all'
        out=...
          simpledata.test_parameters('y_randn',l,w)+...
          simpledata.test_parameters('y_poly', l,w)+...
          simpledata.test_parameters('y_sin',  l,w);
      case 'y_all_T'
        out=...
          simpledata.test_parameters('y_randn',l,w)+...
          simpledata.test_parameters('y_poly', l,w)+...
          simpledata.test_parameters('y_sin_T',l,w)+...
          simpledata.test_parameters('y_cos_T',l,w);
      case {'randn','trend','sin','all','all_T'}
        args=simpledata.test_parameters('args',l,w);
        out=simpledata(...
          simpledata.test_parameters('x',l,w),...
          simpledata.test_parameters(['y_',field],l,w),...
          args{:}...
        );
      case 'mask'
        out=rand(l,1)<0.9;
      case 'no-mask'
        out=true(l,1);
      case 'columns'
        %NOTICE: this return from 1 to w columns
        out=[1,find(randn(1,w-1)>0)+1];
      case 'obj'
        args=simpledata.test_parameters('args',l,w);
        out=simpledata(...
          simpledata.test_parameters('x',    l,w),...
          simpledata.test_parameters('y_all',l,w),...
          'mask',simpledata.test_parameters('mask',l,w),...
          args{:}...
        );
      otherwise
        error(['cannot understand field ''',field,'''.'])
      end
    end
    function test(method,l,w)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      if ~exist('l','var') || isempty(l)
        l=simpledata.test_parameters('l');
      end
      if ~exist('w','var') || isempty(w)
        w=simpledata.test_parameters('w');
      end
      %get common parameters
         args=simpledata.test_parameters('args',l,w);
      columns=simpledata.test_parameters('columns',l,w);
      %init object
      a=simpledata.test_parameters('obj',l,w);

      switch method
        case 'all'
          for i={'plot','append','median','trim','slice','interp','detrend','outlier','plus','pardecomp','augment','vmean'}
           simpledata.test(i{1},l);
          end
        case 'plot'
          figure
          a.plot('columns', columns)
        case 'append'
          b=a.append(...
            simpledata(-l-1:-1,simpledata.test_parameters('y_all',l+1,w),...
              'mask',simpledata.test_parameters('mask',l+1,w),...
              args{:}...
            )...
          );
          figure
          subplot(2,1,1)
          a.plot('columns', columns)
          title('original')
          subplot(2,1,2)
          b.plot('columns', columns)
          title(method)
        case 'median'
          lines1=cell(size(columns));lines1(:)={'-o'};
          lines2=cell(size(columns));lines2(:)={'-x'};
          figure
          a.median(10 ).plot('columns', columns,'line',lines1);
          a.medfilt(10).plot('columns', columns,'line',lines2);
          legend('median','medfilt')
          title(method)
        case 'trim'
          b=a.trim(...
            round(-l/2),...
            round( l/3)...
          );
          figure
          subplot(2,1,1)
          a.plot('columns', columns)
          title('original')
          subplot(2,1,2)
          b.plot('columns', columns)
          title(method)
        case 'slice'
          b=a.slice(...
            round(l*0.3),...
            round(l*0.4)...
          );
          figure
          a.plot('columns', columns,'line',repmat({'+-'},numel(columns),1))
          b.plot('columns', columns,'line',repmat({'o-'},numel(columns),1))
          legend_str={};
          for i=1:numel(columns);legend_str{end+1}=['original-',num2str(i)]; end %#ok<AGROW>
          for i=1:numel(columns);legend_str{end+1}=['sliced-',  num2str(i)]; end %#ok<AGROW>
          %%
          %%
          legend(legend_str)
          title(method)
        case 'interp'
          figure
          legend_str={};
          %define interpolated time domain
          t_interp=round(-l*0.1):0.4:round(l*0.7);
          %add some gaps to the data
          mask=a.mask;
          mask(a.idx(round(l*0.45)):a.idx(round(l*0.5)))=false;
          a=a.mask_and(mask);
          %plot it
          legend_str{end+1}='original';
          a.plot('columns',1,'line',{'o:'});
          legend_str{end+1}='linear interp over gaps narrower than=0';
          a.interp(...
            t_interp,'interp_over_gaps_narrower_than',0 ...
          ).plot('columns',1,'line',{'x-'});
          legend_str{end+1}='spline interp over gaps narrower than 3';
          a.interp(...
            t_interp,'interp_over_gaps_narrower_than',3 ...
          ).plot('columns',1,'line',{'d-'});
          legend_str{end+1}='spline interp over gaps narrower than \infty';
          a.interp(...
            t_interp,'interp_over_gaps_narrower_than',inf ...
          ).plot('columns',1,'line',{'s-'});
          legend(legend_str)
        case 'detrend'
          b=a.detrend;
          figure
          subplot(2,1,1)
          a.plot('columns', columns)
          title('original')
          subplot(2,1,2)
          b.plot('columns', columns)
          title(method)
        case 'outlier'
          %TODO: need to revise the outlier method, this test is giving strange results
          outlier_iter=2;
          figure
          a=a.assign(simpledata.test_parameters('y_randn',l,w));
          a.plot('columns', 1,'line',{'x-'});
          legend_str{1}=a.descriptor;
          [a,o]=a.outlier('outlier_iter',outlier_iter,'outlier_sigma',3);
          a.plot(...
            'columns', 1,...
            'masked',true,...
            'line',{'o'}...
          );
          legend_str{2}=[a.descriptor,' w/out outliers'];
          for i=1:outlier_iter
            if any(o{i}.mask)
              o{i}.plot(...
                'columns', 1,...
                'masked',true,...
                'line',{'d'}...
              );
              legend_str{i+2}=o{i}.descriptor;
            end
          end
          legend(legend_str)
          title(method)
        case 'plus'
          a=simpledata(0:l,simpledata.test_parameters('y_all',l+1,w),...
            'mask',[true(0.2*l+1,1);false(0.2*l,1);true(0.6*l,1)],...
            args{:}...
          );
          b=simpledata(-l/2:2:3/2*l,simpledata.test_parameters('y_all',l+1,w),...
            'mask',[true(0.3*l+1,1);false(0.2*l,1);true(0.2*l,1);false(0.1*l,1);true(0.2*l,1)],...
            args{:}...
          );
          figure
          a.plot('columns',1,'line',{'+-'})
          b.plot('columns',1,'line',{'x-'})
          c=a+b;
          c.plot('columns',1,'line',{'o-'})
          legend('a','b','a+b')
          title(method)
        %the tests below do not require the object 'a', defined above
        case 'pardecomp'
          %NOTICE: this uses a deprecated method
          t=simpledata.test_parameters('x',l);
          y=simpledata.test_parameters('y_all_T',l,w);
          a=num.pardecomp(t,y(:,1),...
            'polynomial',numel(simpledata.test_parameters('y_poly_scale'))-1,...
            'sinusoidal',simpledata.test_parameters('T',l)...
          );
          disp(['sin_periods : ',num2str(simpledata.test_parameters('T',l))])
          disp(['poly_coeffs : ',num2str(simpledata.test_parameters('y_poly_scale'))])
          disp(['sin_coeffs  : ',num2str(simpledata.test_parameters('y_sin_scale'))])
          disp(['cos_coeffs  : ',num2str(simpledata.test_parameters('y_cos_scale'))])
          disp(['randn_scale : ',num2str(simpledata.test_parameters('y_randn_scale'))])
          num.plot_pardecomp(a)
        case 'augment'
          a=simpledata([1:4,5 7 9 11],[11 12 13 14 15 17 19 31]',...
            'mask',logical([1 1 0 0 1 1 0 0])...
          );
          b=simpledata([1:4,6,8,10 12],[21 22 23 24 26 28 40 42]',...
            'mask',logical([1 0 1 0 1 0 1 0])...
          );
          disp('Augment test: a')
          a.peek
          disp('Augment test: b')
          b.peek
          disp('Augment test: a.augment(b,common=[true,false],new=true,old=true)')
          a.augment(b,'quiet',true,'common',true,'new',true,'old',true,'skip_gaps',false).peek
          disp('---')
          a.augment(b,'quiet',true,'common',false,'new',true,'old',true,'skip_gaps',false).peek

          disp('Augment test: a.augment(b,common=[true,false],new=false,old=true)')
          a.augment(b,'quiet',true,'common',true,'new',false,'old',true,'skip_gaps',false).peek
          disp('---')
          a.augment(b,'quiet',true,'common',false,'new',false,'old',true,'skip_gaps',false).peek

          disp('Augment test: a.augment(b,common=true,new=true,old=true,skip_gaps=[false,true])')
          a.augment(b,'quiet',true,'common',true,'new',true,'old',true,'skip_gaps',false).peek
          disp('---')
          a.augment(b,'quiet',true,'common',true,'new',true,'old',true,'skip_gaps',true).peek
        case 'vmean'
          figure('visible','on');
          a=cell(1,w*2);
          legend_str=cell(1,numel(a)+1);
          for i=1:numel(a)
            a{i}=simpledata.randn(1:l,w);
            a{i}.plot('columns',1,'line',{'.'})
            legend_str{i}=['a',num2str(i)];
          end
          b=simpledata.vmean(a);
          b.plot('columns',1)
          legend_str{end}='mean(a_i)';
          legend(legend_str)
        otherwise
          error(['Cannot handle test method ''',method,'''.'])
      end

    end
  end
  methods
    %% constructor
    function obj=simpledata(x,y,varargin)
      % input parsing
      p=machinery.inputParser;
      p.addRequired( 'x'      ,                  @(i) simpledata.valid_x(i));
      p.addRequired( 'y'      ,                  @(i) simpledata.valid_y(i));
      p.addParameter('mask'   ,true(size(x(:))), @(i) simpledata.valid_mask(i));
      %create argument object, declare and parse parameters, save them to obj
      [~,~,obj]=varargs.wrap('sinks',{obj},'parser',p,'sources',{simpledata.parameters('obj')},'mandatory',{x,y},varargin{:});
      %assign (this needs to come before the parameter check, so that sizes are known)
      varargin=cells.vararginclean(varargin,'t');
      obj=obj.assign(y,'x',x,varargin{:});
      %rowize labels and units
      obj.labels =transpose(obj.labels( :));
      obj.units=transpose(obj.units(:));
      % check parameters
      for i=1:simpledata.parameters('length')
        obj=check_annotation(obj,simpledata.parameters('name',i));
      end
    end
    %NOTICE: this method blindly assigns data
    function obj=assign_x(obj,x)
      %propagate x
      obj.x=x(:);
      %update length
      obj.length=numel(obj.x);
    end
    %NOTICE: this method blindly assigns data
    function obj=assign_y(obj,y)
      %propagate y
      obj.y=y;
      %update width
      obj.width=size(y,2);
    end
    %NOTICE: this method blindly assigns data
    function obj=assign_mask(obj,mask)
      obj.mask=mask;
    end
    %NOTICE: unlike the assign_* methods above, this method does some checking
    function obj=assign(obj,y,varargin)
      p=machinery.inputParser;
      p.addRequired( 'y'          ,         @(i) simpledata.valid_y(i));
      p.addParameter('x'          ,obj.x,   @(i) simpledata.valid_x(i));
      p.addParameter('mask'       ,obj.mask,@(i) simpledata.valid_mask(i));
      p.addParameter('reset_width',false,   @isscalar);
      p.addParameter('monotonize','none',   @ischar);
      % parse it
      p.parse(y,varargin{:});
      % ---- x ----
      if numel(p.Results.x)~=size(y,1)
        error(['number of elements of input ''x'' (',num2str(numel(p.Results.x)),...
          ') must be the same as the number of rows of input ''y'' (',num2str(size(y,1)),').'])
      end
      %propagate x
      obj=obj.assign_x(p.Results.x(:));
      % ---- y ----
      if ~logical(p.Results.reset_width) && ~isempty(obj.width) && size(y,2) ~= obj.width
        error(['data width changed from ',num2str(obj.width),' to ',num2str(size(y,2)),'.'])
      end
      if obj.length~=size(y,1)
        error(['data length different than size(y,2), i.e. ',...
          num2str(obj.length),' ~= ',num2str(size(y1m)),'.'])
      end
      %propagate y
      obj=obj.assign_y(y);
      % ---- mask ----
      %check if explicit mask was given
      if ~any(strcmp(p.UsingDefaults,'mask'))
        %make sure things make sense
        assert(obj.length==numel(p.Results.mask),[' ',...
          'number of elements of input ''mask'' (',num2str(numel(p.Results.mask)),...
          ') must be the same as the data length (',num2str(obj.length),').'])
        %propagate mask
        obj=obj.assign_mask(p.Results.mask(:));
      else
        %using existing mask (the default), data length may have changed: set mask from y
        obj=obj.mask_reset;
      end
      %sanitize done inside mask_update
      obj=mask_update(obj);
      %monotonize the data if requested
      switch lower(p.Results.monotonize)
      case 'none'
        %do nothing
      otherwise
        obj=obj.monotonize(p.Results.monotonize);
      end
    end
    function obj=assign_tx_mask(obj,y,tx,mask,varargin)
      if isprop(obj,'t'); obj=obj.assign(  y,'mask',mask,'t',tx,varargin{:});
      else              ; obj=obj.assign_x(y,'mask',mask,'x',tx,varargin{:});
      end
    end
    %NOTICE: obj2=simpledata(t,y,obj1.varargin{:}) is preferable to obj2=simpledata(t,y).copy_metadata(obj1)
    %TODO: check if it's possible to ditch the copy_metadata members because of the reason above
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      pn=[simpledata.parameters('list');more_parameters(:)];
      for i=1:numel(pn)
        if isprop(obj,pn{i}) && isprop(obj_in,pn{i})
          obj.(pn{i})=obj_in.(pn{i});
        end
      end
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      warning off MATLAB:structOnObject
      out=structs.filter(struct(obj),[simpledata.parameters('list');more_parameters(:)]);
      warning on MATLAB:structOnObject
      out=structs.rm_empty(out);
    end
    function out=varargin(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      out=varargs(obj.metadata(more_parameters)).varargin;
    end
    %% info methods
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      for i=1:simpledata.parameters('length')
        obj.disp_field(simpledata.parameters('name',i),tab);
      end
      %some parameters are not listed
      more_parameters={'length','width'};
      for i=1:numel(more_parameters)
        obj.disp_field(more_parameters{i},tab);
      end
      %add some more info
      obj.disp_field('nr gaps',    tab,sum(~obj.mask))
      obj.disp_field('first datum',tab,sum(obj.x(1)))
      obj.disp_field('last datum', tab,sum(obj.x(end)))
      switch class(obj)
        case 'simplegrid'
          %do nothing, the str-nl mode does not work with simplegrid.stats
          %TODO: figure out a way to show grid statistics with the print method
        otherwise
        obj.disp_field('statistics', tab,[10,stats(obj,...
          'mode','str-nl',...
          'tab',tab,...
          'period',seconds(inf),...
          'columns',1:min([obj.peekwidth,obj.width])...
        )])
      end
      obj.peek
    end
    function disp_field(obj,field,tab,value)
      if ~exist('value','var') || isempty(value)
        value=obj.(field);
        if isnumeric(value) || iscell(value)
          value=value(1:min([numel(value),obj.peekwidth]));
        end
      end
      disp([str.tabbed(field,tab),' : ',str.show(transpose(value(:)))])
    end
    function peek(obj,idx,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      if ~exist('idx','var') || isempty(idx)
        idx=[...
          1:min([           obj.peeklength,  round(0.5*obj.length)  ]),...
            max([obj.length-obj.peeklength+1,round(0.5*obj.length)+1]):obj.length...
        ];
      elseif strcmp(idx,'all')
        idx=1:obj.length;
      end
      %check if formatted time stamp is possible
      formatted_time=ismethod(obj,'t') || isfield(obj,'t') || isprop(obj,'t');
      %adapt tab to formatted time
      if tab<19 && formatted_time
        tab_x=19;
      else
        tab_x=tab;
      end
      for i=1:numel(idx)
        out=cell(1,min([obj.peekwidth,obj.width]));
        for j=1:numel(out)
          out{j}=str.tabbed(num2str(obj.y(idx(i),j)),tab,true);
        end
        %use formatted dates is object supports them
        if formatted_time
          x_str=datestr(obj.t(idx(i)),'yyyy-mm-dd HH:MM:SS');
        else
          x_str=num2str(obj.x(idx(i)));
        end
        disp([...
          str.tabbed(x_str,tab_x,true),' ',...
          strjoin(out),' ',...
          str.show(obj.mask(idx(i)))...
        ])
      end
    end
    function out=stats(obj,varargin)
      p=machinery.inputParser;
      p.addParameter('mode',          'struct',          @ischar);
      p.addParameter('frmt',          '%-16.3g',         @ischar);
      p.addParameter('tab',           8,                 @isscalar);
      p.addParameter('minlen',        2,                 @isscalar);
      p.addParameter('columns',       1:obj.width,       @(i) isnumeric(i) || iscell(i));
      p.addParameter('struct_fields',...
        {'min','max','mean','std','rms','meanabs','stdabs','rmsabs','length','gaps','norm'},...
        @(i) iscellstr(i)...
      )
      % parse it
      p.parse(varargin{:});
      %detrend and remove outliers (if there inclinded according to outlier_iter, outlier_sigma and detrend)
      %NOTICE: need to set 'outlier_iter' to 0, otherwise the default in obj.outlier is used in the (common)
      %        case varargin{:} is omissive in the value of this option.
      obj=obj.outlier('outlier_iter',0,varargin{:});
      % type conversions
      if iscell(p.Results.columns)
        columns=cell2mat(p.Results.columns);
      else
        columns=p.Results.columns;
      end
      %trivial call
      switch p.Results.mode
      case {'min','max','mean','std','rms','meanabs','stdabs','rmsabs'}
        %bad things happend when deriving statistics of zero or one data points
        if size(obj.y_masked,1)<max([2,p.Results.minlen])
          out=zeros(1,numel(columns));
          return
        end
      end
      %branch on mode
      switch p.Results.mode
      case 'min';     out=min(     obj.y_masked([],columns));
      case 'max';     out=max(     obj.y_masked([],columns));
      case 'ampl';    out=max(     obj.y_masked([],columns))-min(obj.y_masked([],columns));
      case 'mean';    out=mean(    obj.y_masked([],columns));
      case 'std';     out=std(     obj.y_masked([],columns));
      case 'rms';     out=rms(     obj.y_masked([],columns));
      case 'meanabs'; out=mean(abs(obj.y_masked([],columns)));
      case 'stdabs';  out=std( abs(obj.y_masked([],columns)));
      case 'rmsabs';  out=rms( abs(obj.y_masked([],columns)));
      case 'length';  out=obj.nr_valid*ones(1,numel(columns));
      case 'gaps';    out=obj.nr_gaps* ones(1,numel(columns));
      %one line, two lines, three lines, many lines
      case 'norm';    out=norm(     obj.y_masked([],columns));
      case {'str','str-2l','str-3l','str-nl'}
        out=cell(size(p.Results.struct_fields));
        for i=1:numel(p.Results.struct_fields)
          out{i}=[...
            str.tabbed(p.Results.struct_fields{i},p.Results.tab),' : ',...
            num2str(...
              stats@simpledata(...
                obj,varargin{:},...
                'mode',p.Results.struct_fields{i}...
              ),...
            p.Results.frmt)...
          ];
        end
        %concatenate all strings according to the required mode
        switch p.Results.mode
        case 'str'
          out=strjoin(out,'; ');
        case 'str-2l'
          mid_idx=ceiling(numel(out)/2);
          out=[...
            strjoin(out(1:mid_idx),'; '),...
            10,...
            strjoin(out(mid_idx+1:end),'; ')...
            ];
        case 'str-3l'
          mid_idx=ceiling(numel(out)/3);
          out=[...
            strjoin(out(1:mid_idx),'; '),...
            10,...
            strjoin(out(mid_idx+1:2*mid_idx),'; '),...
            10,...
            strjoin(out(2*mid_idx+1:end),'; ')...
            ];
        case 'str-nl'
          out=strjoin(out,'\n');
        end
      case 'struct'
        for f=p.Results.struct_fields
          out.(f{1})=stats@simpledata(obj,varargin{:},'mode',f{1});
        end
      case 'obj'
        %get structure with requested stats
        s=stats@simpledata(obj,varargin{:},'mode','struct');
        % use correct abcissae
        x_now=mean(obj.tx_masked);
        % use the correct object constructor
        init=str2func(class(obj));
        %loop over all structure fields and create an object
        fields=fieldnames(s);
        if numel(fields)==1
          out=init(x_now,s.(fields{1}),obj.varargin{:});
        else
          for i=1:numel(fields)
            out.(fields{i})=init(x_now,s.(fields{i}),obj.varargin{:});
          end
        end
      otherwise
        error(['unknown mode ''',p.Results.mode,'''.'])
      end
      %bug traps
      if isnumeric(out)
        if any(isnan(out(:)))
          error('BUG TRAP: detected NaN in <out>. Debug needed.')
        end
        if size(out,2)~=numel(columns) && ~strcmp(p.Results.mode,'norm')
          error('BUG TRAP: data width changed. Debug needed.')
        end
      end
    end
    function out=stats2(obj1,obj2,varargin)
      p=machinery.inputParser;
      p.addParameter('mode',       'struct', @ischar);
      p.addParameter('minlen',            2, @isscalar);
      p.addParameter('columns', 1:min([obj1.width,obj2.width]), @(i) isnumeric(i) || iscell(i));
      p.addParameter('struct_fields',...
        {'cov','corrcoef','length','rms'},...
        @(i) iscellstr(i)...
      )
      % parse it
      p.parse(varargin{:});
      %need compatible objects
      compatible(obj1,obj2,varargin{:})
      %need to make the mask match to make sure x_masked is common
      [obj1,obj2]=obj1.mask_match(obj2,'x-domain discrepancy, interpolate objects before calling this method');
      %detrend and remove outliers (if there inclinded according to outlier_iter, outlier_sigma and detrend)
      %NOTICE: need to set 'outlier_iter' to 0, otherwise the default in obj.outlier is used in the (common)
      %        case varargin{:} is omissive in the value of this option.
      obj1=obj1.outlier('outlier_iter',0,varargin{:});
      obj2=obj2.outlier('outlier_iter',0,varargin{:});
      %trivial call
      switch p.Results.mode
      case {'cov','corrcoef'}
        %bad things happend when deriving statistics of zero or one data points
        if size(obj1.y_masked,1)<=max([2,p.Results.minlen])
          out=nan(1,obj1.width);
          return
        end
      end
      % type conversions
      if iscell(p.Results.columns)
        columns=cell2mat(p.Results.columns);
      else
        columns=p.Results.columns;
      end
      %branch on mode
      switch p.Results.mode
      case 'cov'
        out=num.cov(obj1.y_masked([],columns),obj1.y_masked([],columns));
      case {'corrcoef','corrcoeff'}
        out=num.corrcoef(obj1.y_masked([],columns),obj2.y_masked([],columns));
      case {'rms','std'}
        o=obj1.cols(columns)-obj2.cols(columns);
        out=o.stats('mode',strrep(p.Results.mode,'diff',''));
      case 'length'
        out=(obj1.nr_valid+obj2.nr_valid)*ones(1,numel(columns));
      case 'gaps'
        out=(obj1.nr_gaps+obj2.nr_gaps)*ones(1,numel(columns));
      case 'struct'
        for f=p.Results.struct_fields
          out.(f{1})=obj1.stats2(obj2,varargin{:},'mode',f{1});
        end
      case 'obj'
        %get structure with requested stats
        s=stats2@simpledata(obj1,obj2,varargin{:},'mode','struct');
        % use correct abcissae (obj2 should produce the same x_now, since the time simpledata.merge method was used
        x_now=mean(obj1.tx_masked);
        % use the correct object constructor
        init=str2func(class(obj1));
        %loop over all structure fields and create an object (unless there's only one field)
        fields=fieldnames(s);
        if numel(fields)==1
          out=init(x_now,s.(fields{1}),obj1.varargin{:});
        else
          for i=1:numel(fields)
            out.(fields{i})=init(x_now,s.(fields{i}),obj1.varargin{:});
          end
        end
      otherwise
        error(['unknown mode ''',p.Results.mode,'''.'])
      end
      %bug traps
      if isnumeric(out)
        if any(isnan(out(:)))
          error('BUG TRAP: detected NaN in <out>. Debug needed.')
        end
        if size(out,2)~=numel(columns)
          error('BUG TRAP: data width changed. Debug needed.')
        end
      end
    end
    function out=size(obj)
      out=[obj.length,obj.width];
    end
%     function out=numel(obj)
%       out=obj.length*obj.width;
%     end
%     function out=class(obj)
%       out=class(obj);
%     end
    %% x methods
    function out=x_masked(obj,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      end
      out=obj.x(mask);
    end
    function out=idx(obj,x_now,varargin)
      % sanity
      assert(isvector( x_now),'Input x_now must be a vector.')
      assert(isnumeric(x_now),['Input x_now must be numeric, not ',class(x_now),'.'])
      if numel(x_now)>1
        %making room for outputs
        out=NaN(size(x_now));
        %recursive call
        for i=1:numel(x_now)
          out(i)=obj.idx(x_now(i),varargin{:});
        end
      else
        if isfinite(x_now)
          %compute distance to input
          distance=(obj.x-x_now).^2;
          %handle default find arguments
          if isempty(varargin)
            out=find(distance==min(distance),1,'first');
          else
            out=find(distance==min(distance),varargin{:});
          end
        else
          out=find(~isfinite(obj.x));
        end
      end
    end
    function obj=at_idx(obj,i,varargin)
      if i<0
        i=obj.length+i+1;
      end
      obj=obj.assign(...
        obj.y(i,:),...
        'x',obj.x(i,:),...
        'mask',obj.mask(i,:),...
        varargin{:}...
      );
    end
    function obj=at(obj,x_now,varargin)
      obj=obj.at_idx(obj.idx(x_now,varargin{:}),varargin{:});
    end
    function [obj,idx_add,idx_old,x_old]=x_merge(obj,x_add,y_new)
      if ~exist('y_new','var') || isempty(y_new)
        y_new=NaN;
      end
      %add epochs given in x_add, the corresponding y are NaN and mask are set as gaps
      if nargout==4
        %save x_old
        x_old=obj.x;
      end
      %build extended x-domain
      [x_total,idx_old,idx_add] = simpledata.union(obj.x,x_add);
      %shortcut
      if all(idx_old) && all(idx_add)
        return
      end
      %build extended y and mask
      y_total=ones(numel(x_total),obj.width)*y_new;
      y_total(idx_old,:)=obj.y;
      mask_total=true(size(x_total)); %assume valid mask on all entries, so that y_new can propagate
      mask_total(idx_old)=obj.mask;
      %sanitize monotonicity and propagate
      [x_total,y_total,mask_total]=simpledata.monotonic(x_total,y_total,mask_total);
      %propagate to obj
      obj=obj.assign(y_total,'x',x_total,'mask',mask_total);
    end
    function out=isxavail(obj,x)
      if isscalar(x)
        out=any(obj.x==x);
      else
        out=arrayfun(@(i) obj.isxavail(i),t);
      end
    end
    function out=tx(obj)
      if isprop(obj,'t'); out=obj.t; %#ok<MCNPN>
      else              ; out=obj.x;
      end
    end
    function out=tx_masked(obj,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      end
      if ismethod(obj,'t_masked'); out=obj.t_masked(mask);
      else                     ; out=obj.x_masked(mask);
      end
    end
    %% y methods
    function obj=cols(obj,columns)
      if ~isvector(columns)
        error('input ''columns'' must be a vector')
      end
      %save relevant labels and units
      labels_now =obj.labels( columns);
      units_now=obj.units(columns);
      %retrieve requested columns
      obj=obj.assign(obj.y(:, columns),'reset_width',true);
      obj.labels =labels_now;
      obj.units=units_now;
    end
    function obj=get_cols(obj,columns)
      obj=obj.cols(columns);
    end
    function obj=set_cols(obj,columns,values,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=true(obj.length,1);
      end
      y_now=obj.y;
      if isnumeric(values)
        y_now(mask,columns)=values;
      else
        y_now(mask,columns)=values.y;
      end
      obj=obj.assign(y_now);
    end
    function out=y_masked(obj,mask,columns)
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      end
      if ~exist('columns','var') || isempty(columns)
        columns=1:obj.width;
      end
      if ~isvector(columns)
        error('input ''columns'' must be a vector')
      end
      if any(columns>obj.width)
        error('requested column indeces exceed object width.')
      end
      out=obj.y(mask,columns);
    end
    function obj=set_at(obj,x,y,f)
      if ~exist('f','var') || isempty(f)
        f='assign';
      end
      switch f
      case {'increment','add'  ,'+'}
        f=@(i,j) i+j;
      case {'decrement','minus','-'}
        f=@(i,j) i-j;
      case {'assign'   ,'equal','='}
        f=@(i,j) j;
      otherwise
        error(['Cannot understand input ''f''with value ''',f,'''.'])
      end
      %get availability indeces
      avail_idx=obj.isxavail(x);
      %operate on the entries that are already there
      if any(avail_idx)
        %get the relevant indeces in the inputs
        idx_in=find(avail_idx);
        %get the relevant indeces in the object
        idx_obj=obj.idx(x(idx_in));
        %get the data
        y_obj=obj.y;
        %operate
        y_obj(idx_obj,:)=f(y_obj(idx_obj,:),y(idx_in,:));
        %back-propagate
        obj=obj.assign(y_obj);
      end
      %append the new entries
      if any(~avail_idx)
        x=x(~avail_idx);
        y=y(~avail_idx,:);
        switch class(obj)
        case 'simpledata'
          obj_new=simpledata(x,y,obj.varargin{:});
        case 'simpletimeseries'
          obj_new=simpletimeseries(obj.x2t(x),y,obj.varargin{:});
        otherwise
          error(['Cannot handle object of class ''',class(obj),''', implementation needed (it''s quick, go and do it).'])
        end
        obj=obj.append(obj_new);
      end
    end
    function out=iszero(obj,mask)
      if ~exist('mask','var')
        o=obj.masked;
      else
        o=obj.masked(mask);
      end
      if ~isempty(o)
        out=all(o.y(:)==0);
      else
        out=true;
      end
    end
    function out=iszero_cols(obj,columns)
      if ~exist('columns','var') || isempty(columns)
        columns=1:obj.width;
      end
      out=all(obj.y_masked([],columns)==0);
    end
    %% mask methods
    %NOTICE: these functions only deal with *explicit* gaps, i.e. those in
    %the form of:
    %x=[1,2,3,4,5..], mask=[T,T,F,T,T...] (x=3 being an explict gap).
    %in comparison, implicit gaps are:
    %x=[1,2,4,5..], mask=[T,T,T,T...] (x=3 being an implict gap).
    function obj=mask_and(obj,mask_now)
      obj.mask=obj.mask & mask_now(:);
      obj=obj.mask_update;
    end
    function obj=mask_or(obj,mask_now)
      obj.mask=obj.mask | mask_now(:);
      obj=obj.mask_update;
    end
    function obj=mask_update(obj)
      %sanity
      if numel(obj.mask) ~= size(obj.y,1)
        error('sizes of mask and y are not in agreement.')
      end
      %propagate gaps from and to mask
      obj.mask=obj.mask & all(~isnan(obj.y),2);
      obj.y(~obj.mask,:)=NaN;
      %sanitize
      obj.check_sd
    end
    function obj=mask_reset(obj)
      obj.mask=all(~isnan(obj.y),2);
    end
    function [out,nr_epochs]=gap_length(obj)
      %NOTICE: this needs explicit gaps to work as expected
      %inits
            out=zeros(obj.length,1);
      nr_epochs=zeros(obj.length,1);
      gap_length=0;
      gaps_start_idx=0;
      %go through the mask
      for i=1:obj.length
        %check for gap
        if ~obj.mask(i)
          %increment counters
          gap_length=gap_length+1;
          %init start index (unless already done)
          if gaps_start_idx==0
            gaps_start_idx=i;
          end
        else
          %check if coming out of a gap
          if gap_length>0
            idx=gaps_start_idx:i-1;
            %save gap length and nr of gaps
                  out(idx)=obj.x(i)-obj.x(max([1,gaps_start_idx-1]));
            nr_epochs(idx)=gap_length;
            %reset counters
            gap_length=0;
            gaps_start_idx=0;
          end
        end
      end
      %save gap length at the end of the x-domain
      if gap_length>0
        idx=gaps_start_idx:obj.length;
        %save gap length and nr of gaps
              out(idx)=obj.x(end)-obj.x(gaps_start_idx-1);
        nr_epochs(idx)=gap_length;
      end
      %sanity
      if any(out(obj.mask)~=0)
        error('found non-zero gaps lengths outside of gaps. Debug needed!')
      end
      if any(out(~obj.mask)==0)
        error('found zero lengths inside of gaps. Debug needed!')
      end
    end
    function out=nr_gaps(obj)
      out=sum(~obj.mask);
    end
    function out=is_all_gaps(obj)
      out=all(~obj.mask);
    end
    function out=nr_valid(obj)
      out=sum(obj.mask);
    end
    function out=is_all_valid(obj)
      out=all(obj.mask);
    end
    function [obj1,obj2]=mask_match(obj1,obj2,errmsg)
      if ~exist('errmsg','var') || isempty(errmsg)
        errmsg='x-domain discrepancy, cannot match masks';
      end
      assert(obj1.length==obj2.length && all(simpledata.isx('==',obj1.x,obj2.x,min([obj1.x_tol,obj2.x_tol]))),errmsg)
      obj1.mask=obj1.mask & obj2.mask;
      obj2.mask=obj2.mask & obj1.mask;
      obj1=obj1.mask_update;
      obj2=obj2.mask_update;
    end
    function obj=masked(obj,mask)
      if ~exist('mask','var') || isempty(mask)
        mask=obj.mask;
      else
        mask=obj.mask&&mask;
      end
      if ~any(mask)
        obj=[];
      else
        obj=obj.assign(obj.y(mask,:),'x',obj.x(mask));
        %sanity
        assert(all(obj.mask) && ~any(isnan(obj.y(:))),...
          'making operation failed: found non-unitary mask entries and/or NaNs in the data.')
      end
    end
    function out=ismaskcommon(obj1,obj2)
      out=obj1.length==obj2.length && obj1.nr_gaps==obj2.nr_gaps && all(obj1.x_masked==obj2.x_masked);
    end
    function obj=noedgegaps(obj)
      %NOTICE: cannot set the loop criteria to obj.length>0 because cannot handle empty y
      while obj.length>1 && ~obj.mask(1)
        m=[false;true(obj.length-1,1)];
        obj=obj.assign(obj.y(m,:),'x',obj.x(m));
      end
      while obj.length>1 && ~obj.mask(end)
        m=[true(obj.length-1,1);false];
        obj=obj.assign(obj.y(m,:),'x',obj.x(m));
      end
    end
    function obj_nan=nan(obj)
      %duplicates an object, setting y to nan
      obj_nan=obj.and(false(size(obj.mask))).mask_update;
    end
    %% invalid methods
    function obj=demasked(obj,invalid)
      if ~exist('invalid','var') || isempty(invalid)
        invalid=obj.invalid;
      end
      %replace NaNs with invalid entries
      y_now=obj.y;
      y_now(isnan(y_now))=invalid;
      %back-propagate
      obj=obj.assign(y_now);
      %paranoid sanity
      assert(obj.nr_gaps == 0,['BUG TRAP: ',...
        'there are still explicit gaps in the data. Debug needed!'])
    end
    function obj=remasked(obj,invalid)
      if ~exist('invalid','var') || isempty(invalid)
        invalid=obj.invalid;
      end
      %replace invalids with NaNs entries
      y_now=obj.y;
      y_now(y_now==invalid)=Nan;
      %back-propagate
      obj=obj.assign(y_now);
    end
    %% management methods
    function check_sd(obj)
      %sanitize
      if obj.length~=numel(obj.x)
        error('discrepancy in length of x')
      end
      if obj.length~=numel(obj.mask)
        error('discrepancy in length of mask')
      end
      if obj.length~=size(obj.y,1)
        error('discrepancy in length of y')
      end
      if obj.width~=size(obj.y,2)
        error('discrepancy in width of y')
      end
      if any(~all(isnan(obj.y),2)~=obj.mask)
        error('discrepancy between mask and NaNs in y.')
      end
      if any(isnan(obj.y),2)~=all(isnan(obj.y),2)
        error('found epochs with some components of y as NaNs.')
      end
    end
    function obj=check_annotation(obj,annotation_name)
      check=true;
      if ~isempty(obj.width) && ...
          iscell(obj.(annotation_name)) && ...
          obj.width~=numel(obj.(annotation_name))
        if  numel(obj.(annotation_name))==1 && ...
          isempty(obj.(annotation_name){1})
          %annotation is cell, with unit length and that entry is empty,
          %i.e. the default is still there. Reset according to width
          obj.(annotation_name)=cell(obj.width,1);
          %populate with empty strings (some sanity checks need cellstr)
          obj.(annotation_name)(:)={''};
        else
          check=false;
        end
      end
      if ~check
        error(['discrepancy in between length of ',annotation_name,' and of width of y'])
      end
    end
    function obj=monotonize(obj,mode)
      [x_now,y_now,mask_now]=simpledata.monotonic(obj.x,obj.y,obj.mask,mode);
      %need blind assigns since this method can be called during initialization
      obj=obj.assign_y(y_now);
      obj=obj.assign_x(x_now);
      obj=obj.assign_mask(mask_now);
    end
    %% edit methods
    function obj=remove(obj,rm_mask)
      %remove epochs given in rm_mask, by eliminating them from the data,
      %not by setting them as gaps.
      obj=obj.assign(obj.y(~rm_mask,:),'x',obj.x(~rm_mask),'mask',obj.mask(~rm_mask));
    end
    function obj=trim(obj,start,stop)
      %remove data outside start/stop
      assert(start<=stop,'input ''start'' must refer to an abcissae before input ''stop''.')
      %trivial call
      if stop<obj.x(1) || obj.x(end)<start || all (obj.x < start | stop < obj.x)
        obj=[];
      else
        obj=obj.remove(obj.x < start | stop < obj.x);
      end
    end
    function obj=slice(obj,start,stop)
      %delete data between start/stop
      if start>stop
        error('input ''start'' must refer to a abcissae before input ''stop''.')
      end
      obj=obj.remove(start < obj.x & obj.x < stop);
      %add gap extremities
      gap_extremities=obj.idx([start,stop]);
      obj.y(   gap_extremities,:)=NaN;
      obj.mask(gap_extremities)=false;
    end
    function obj=interp(obj,x_now,varargin)
      % parse inputs
      p=machinery.inputParser;
      p.addParameter('interp_over_gaps_narrower_than',0,@(i) num.isscalar(i) && ~isnan(i))
      p.addParameter('interp1_args',cell(0),@iscell);
      % parse it
      p.parse(varargin{:});
      % need to add requested x-domain, so that implicit gaps can be handled
      %NOTICE: x_now = obj.x(idx_now)
      [obj,idx_now]=obj.x_merge(x_now);
      %interpolate over all gaps, use all available information
      switch numel(obj.x_masked)
      case 0
        y_now=nan(numel(x_now),obj.width);
      case 1
        y_now=ones(numel(x_now),1)*obj.y_masked;
      otherwise
        y_now=interp1(obj.x_masked,obj.y_masked,x_now(:),p.Results.interp1_args{:});
      end
      % create mask with requested gap interpolation scheme
      if p.Results.interp_over_gaps_narrower_than <= 0
        %this interpolates over all gaps, i.e. keeps all interpolated data,
        %bridging over gaps of any size
        mask_now=true(size(x_now(:)));
      elseif ~isfinite(p.Results.interp_over_gaps_narrower_than)
        %this keeps all gaps, i.e. only those points wich are common on
        %both x-domains AND are not gaps will appear in the ouput
        mask_now=obj.mask(idx_now);
      else
        %NOTICE: this needs explicit gaps to work as expected!
        %compute gap length and create a mask from that, i.e. it is
        %possible to control the maximum length of a gap in the output
        mask_add=obj.gap_length<p.Results.interp_over_gaps_narrower_than;
        mask_now=mask_add(idx_now);
      end
      % handle empty data
      if sum(mask_now)<1
        obj=obj.assign(nan(size(obj.y)));
        return
      end
      %propagate requested x-domain only (mask is rebuilt from NaNs in
      obj=obj.assign(y_now,'x',x_now,'mask',mask_now);
    end
    function [obj,S]=polyfit(obj,order,x_out)
      %parameters
      S_fieldnames={'R','df','normr','p'};
      %defaults
      if ~exist('x_out','var') || isempty(x_out)
        x_out=obj.x;
      end
      %branch on class of order
      switch class(order)
      case 'double'
        %make room for polyfit struct
        S=cell(0);
        %set operational flag
        polyfit_flag=true;
      case 'cell'
        %rename
        S=order;
        %sanity
        assert(cellfun(@(i) isstruct(i) && all(isfield(i,S_fieldnames)),S),...
          ['Cannot handle cell array ''S'' unless it contains fields ''',...
          strjoin(S_fieldnames,''','''),'''.'])
        %patch order
        order=numel(S{1}.p);
        %set operational flag
        polyfit_flag=false;
      end
      %copy the data
      obj.descriptor=['order ',num2str(order),' polyfit of ',obj.descriptor];
      if polyfit_flag
        %get the masked data
        y_now=obj.y_masked;
        x_now=obj.x_masked;
      end
      %make room for data
      y_polyfitted=zeros(numel(x_out),obj.width);
      %polyfit for all columns
      for i=1:obj.width
        if polyfit_flag
          [p,S{i}]=polyfit(x_now,y_now(:,i),order);
          S{i}.p=p;
        end
        y_polyfitted(:,i)=polyval(S{i}.p,x_out);
      end
      %propagate the data
      obj=obj.assign(y_polyfitted,'x',x_out);
    end
    function obj=detrend(obj,mode,disp)
      if ~exist('mode','var') || isempty(mode)
        mode='linear';
      end
      if ~exist('disp','var') || isempty(disp)
        disp=false;
      end
      %trivial call
      if str.none(mode)
        return
      end
      %inform user
      if disp
        str.say(mode,'detrending of',obj.descriptor)
      end
      if ischar(mode)
        switch mode
        case 'cubic';     n=3;
        case 'quadratic'; n=2;
        case 'linear';    n=1;
        case 'constant';  n=0;
        otherwise
          if contains(mode,'poly')
            %determine polynomial order to be fitted
            o=str2double(strrep(mode,'poly',''));
            %polyfit the data
            obj_polyfitted=obj.polyfit(o);
            %subtract polyfitted data from input data
            obj=obj-obj_polyfitted;
          else
            error(['unknown mode ''',mode,'''.'])
          end
        end
      else
        n=mode;
      end
      assert(isnumeric(n) && isscalar(n),['Input ''mode'' must be char or integer, not ',class(n),'.'])
      %detrend the whole time series (NOTICE: this does not handle ungaped segments independently)
      obj.y=detrend(obj.y,n,'omitnan','Continuous',false,'SamplePoints',obj.x);
      % sanitize
      obj.check_sd
    end
    %NOTICE: also does detrending (defaults to none) because that operation is often needed to isolated outliers
    function [obj,outliers]=outlier(obj,varargin)
      v=varargs.wrap('sources',{....
        {...
          'outlier_sigma',    obj.outlier_sigma, @isnumeric;...
          'outlier_iter',     1,                 @isfinite;...
          'detrend',     'none',                 @ischar;...
          'disp_flag',     true,                 @islogical;...
        },...
      },varargin{:});
      %handle iterative outlier removal and trivial call
      switch v.outlier_iter
      case 0
        %maybe detrending is requested anyway
        obj=obj.detrend(v.detrend,v.disp_flag);
        outliers=[];
        return;
      case 1
      % do nothing
      otherwise
        %optional outputs: do things quicker if the value of the outliers is not needed.
        if nargout>1
          outliers=cell(1,v.outlier_iter);
          for i=1:v.outlier_iter
            [obj,outliers{i}]=obj.outlier(varargin{:},'outlier_iter',1);
            outliers{i}.descriptor=[outliers{i}.descriptor,' iteration ',num2str(i)];
          end
        else
          for i=1:v.outlier_iter
            obj=obj.outlier(varargin{:},'outlier_iter',1);
          end
        end
        return
      end
      %handle inputs
      if isscalar(v.outlier_sigma)
        v.outlier_sigma=v.outlier_sigma*ones(1,obj.width);
      else
        assert(numel(v.outlier_sigma) ~= obj.width,...
          ['input <outlier_sigma> must have the same number of entries as data width (',...
          num2str(obj.witdh),'), not ',num2str(numel(v.outlier_sigma)),'.']...
        )
      end
      %detrend data, handles the 'detrend' option internally
      obj=obj.detrend(v.detrend,v.disp_flag);
      %inform user
      if v.disp_flag
        str.say('removing outliers larger than ',mean(v.outlier_sigma),'-STD from',obj.descriptor)
      end
      %create tmp container
      y_data=zeros(obj.length,obj.width);
      if nargout>1; y_outliers=zeros(size(y_data)); end
      %loop over data width, do things quicker if the value of the outliers is not needed.
      if nargout>1
        for i=1:obj.width
          [y_data(:,i),idx]=simpledata.rm_outliers(obj.y(:,i),'outlier_sigma',v.outlier_sigma(i),'outlier_value',0);
          y_outliers(idx,i)=obj.y(idx,i);
        end
      else
        for i=1:obj.width
          y_data(:,i)=simpledata.rm_outliers(obj.y(:,i),'outlier_sigma',v.outlier_sigma(i),'outlier_value',0);
        end
      end
      %sanity
      if nargout>1 && any(any(y_data(obj.mask,:)+y_outliers(obj.mask,:)~=obj.y_masked))
        error('BUG TRAP: failed the consistency check: obj.y=y_data+y_outliers. Debug needed!')
      end
      %optional outputs
      if nargout>1
        %copy object to preserve metadata
        outliers=obj;
        %propagate gaps
        y_outliers(~obj.mask,:)=NaN;
        %propatate (mask is derived from gaps in y_outliers)
        outliers=outliers.assign(y_outliers,'x',outliers.x,'mask',all(y_outliers~=0,2));
        %update descriptor
        outliers.descriptor=['outliers of ',obj.descriptor];
        %do not plot zeros, since they refer to other components, not the
        %component where the outlier was detected
        outliers.plot_zeros=false;
      end
      %propagate (mask is updated inside)
      obj=obj.assign(y_data,'mask',all(y_data~=0,2));
    end
    function obj=mean(obj,n)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      obj=obj.segstat(@mean,n);
    end
    function obj=median(obj,n)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      obj=obj.segstat(@median,n);
    end
    function obj=rms(obj,n)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      obj=obj.segstat(@rms,n);
    end
    function obj=std(obj,n)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      obj=obj.segstat(@std,n);
    end
    function obj=segstat(obj,op,n)
      if ~exist('n','var') || isempty(n)
        n=obj.length;
      end
      %sanity on op
      if ~isa(op,'function_handle')
        op=str2func(op);
      end
      valid_fun={'std','rms','mean','median','min','max','var'};
      assert(any(strcmp(func2str(op),valid_fun)),...
        ['Need input ''op'' to be (a function handle or string to) one of ',strjoin(valid_fun,','),'.'])
      %create x-domain of segmented data, cutting it into segments of
      %length n, putting each segment in one column of matrix t (x
      %increases first row-wise, then apply the statistic in 'mode' column-wise)
      x_temp=nan(n,ceil(obj.length/n));
      x_temp(1:obj.length)=obj.x;
      x_mean=transpose(mean(x_temp,'omitnan'));
      %make room for medianed data
      y_out=nan(ceil(obj.length/n),obj.width);
      %compute media of the data
      s.msg=[' computing mean and decimating every ',num2str(n),' data points.'];s.n=obj.width;
      for i=1:obj.width
        %cut the data into segments of length n, putting each segment in
        %one column of matrix y_seg (x increases first row-wise, then
        %column-wise)
        y_seg=nan(n,ceil(obj.length/n));
        y_seg(1:obj.length)=obj.y(:,i);
        switch func2str(op)
        case {'min','max'}
          y_out(:,i)=transpose(op(y_seg,[],'omitnan'));
        otherwise
          y_out(:,i)=transpose(op(y_seg,'omitnan'));
        end
        %user feedback
        s=time.progress(s,i);
      end
      %sanity
      if numel(x_mean) ~=size(y_out,1)
        error('x-domain length inconsistent with data length, debug needed!')
      end
      %propagate the x-domain and data
      obj=assign(obj,y_out,'x',x_mean);
      %update descriptor
      obj.descriptor=[func2str(op),' of ',obj.descriptor];
    end
    %% multiple object manipulation
    function out=isxequal(obj1,obj2)
      %NOTICE: this also handles the single-object operation
      if ismethod(obj1,'istequal') && ( isdatetime(obj2) || ismethod(obj2,'istequal'))
        out=obj1.istequal(obj2);
      elseif isnumeric(obj2)
        out=obj1.length==numel(obj2) && ~any(~simpledata.isx('==',obj1.x,obj2(:),obj1.x_tol));
      else
        out=obj1.length==obj2.length && ~any(~simpledata.isx('==',obj1.x,obj2.x,min([obj1.x_tol,obj2.x_tol])));
      end
    end
    function compatible(obj1,obj2,varargin)
      %This method checks if the objects are referring to the same
      %type of data, i.e. the data length is not important.
      % parse inputs
      p=machinery.inputParser;
      p.addParameter('compatible_parameters',simpledata.compatible_parameter_list,@iscellstr)
      p.addParameter('check_width',true,@islogical);
      p.addParameter('skip_par_check',{''},@iscellstr)
      % parse it
      p.parse(varargin{:});
      %basic sanity
%       if ~strcmp(class(obj1),class(obj2))
%         error(['incompatible objects: different classes'])
%       end
      if p.Results.check_width && (obj1.width ~= obj2.width)
        error('incompatible objects: different number of columns')
      end
      %shorter names
      par=p.Results.compatible_parameters;
      for i=1:numel(par)
        % if a parameter is empty, no need to check it
        if ( iscell(obj1.(par{i})) && isempty([obj1.(par{i}){:}]) ) || ...
           ( ischar(obj1.(par{i})) && isempty( obj1.(par{i})    ) ) || ...
           ( iscell(obj2.(par{i})) && isempty([obj2.(par{i}){:}]) ) || ...
           ( ischar(obj2.(par{i})) && isempty( obj2.(par{i})    ) )
          continue
        end
        if ~cells.isincluded(p.Results.skip_par_check,par{i}) && ~isequal(obj1.(par{i}),obj2.(par{i}))
          error(['discrepancy in parameter ',par{i},'.'])
        end
      end
    end
    function [obj1_out,obj2_out,idx1,idx2]=merge(obj1,obj2,y_new)
      %NOTICE:
      % - idx1 contains the index of the x in obj1 that were added from obj2
      % - idx2 contains the index of the x in obj2 that were added from obj1
      % - no data is propagated between objects, only the time domain is changed!
      %add as gaps in obj2 those x that are in obj1 but not in obj2
      % - y_new sets the value of the data at the new entries of x, both obj1
      %   and obj2 (default to NaN)
      %NOTE: no interpolation is done between the objects, only
      %      the x-domain is made in agreement between them
      if ~exist('y_new','var') || isempty(y_new)
        y_new=NaN;
      end
      [obj1_out,idx2]=obj1.x_merge(obj2.x,y_new);
      if nargout>1
        %add as gaps in obj2 those x that are in obj1 but not in obj2
        [obj2_out,idx1]=obj2.x_merge(obj1.x,y_new);
        %sanity on outputs
        if obj1_out.length~=obj2_out.length && obj1_out.isxsame(obj2_out)
          error('BUG TRAP: merge operation failed.')
        end
      end
    end
    function [obj1,obj2]=interp2(obj1,obj2,varargin)
      %extends the x-domain of both objects to be in agreement
      %with the each other. The resulting x-domains possibly have
      %numerous gaps, which are interpolated over (interpolation
      %scheme and other options can be set in varargin).
      compatible(obj1,obj2,varargin{:})
      %trivial call
      if isxequal(obj1,obj2)
        return
      end
      %build extended x-domain
      x_total=simpledata.union(obj1.x,obj2.x);
      %interpolate obj1 over x_total
      obj1=obj1.interp(x_total,varargin{:});
      %interpolate obj2 over x_total
      obj2=obj2.interp(x_total,varargin{:});
    end
    function [obj,idx1,idx2]=append(obj1,obj2,over_write_flag,varargin)
      if ~exist('over_write_flag','var') || isempty(over_write_flag)
        over_write_flag=false;
      end
      %trivial call
      if isempty(obj2)
        idx1=true( obj1.length,1);
        idx2=false(obj1.length,1);
        obj=obj1;
        return
      end
      %sanity
      compatible(obj1,obj2,varargin{:})
      %append with awareness
      if all(obj1.x(end) < obj2.x)
        x_now=[obj1.x;obj2.x];
        y_now=[obj1.y;obj2.y];
        mask_now=[obj1.mask;obj2.mask];
        idx1=[true( obj1.length,1);false(obj2.length,1)];
        idx2=[false(obj1.length,1);true( obj2.length,1)];
      elseif all(obj2.x(end) < obj1.x)
        x_now=[obj2.x;obj1.x];
        y_now=[obj2.y;obj1.y];
        mask_now=[obj2.mask;obj1.mask];
        idx1=[false(obj2.length,1);true( obj1.length,1)];
        idx2=[true( obj2.length,1);false(obj1.length,1)];
      else
        %retrieve union of time domains
        [x_now,idx1,idx2] = simpledata.union(obj1.x,obj2.x);
        %build appended data
        y_now(idx1,:)=obj1.y;
        y_now(idx2,:)=obj2.y;
        %build mask
        mask_now(idx1)=obj1.mask;
        mask_now(idx2)=obj2.mask;
        %get common data entries
        ic=idx1 & idx2;
        %check data integrity (on y and mask)
        if any(ic)
          %build mask
          m1=true(numel(ic),1);m2=m1;
          m1(idx1)=obj1.mask;
          m2(idx2)=obj2.mask;
          %retrieve common data
          y1(idx1,:)=obj1.y;
          y2(idx2,:)=obj2.y;
          %mask common data entries
          ic=ic & m1 & m2;
          %check it (unless not needed)
          if ~over_write_flag
            assert(all(all(y1(ic,:)==y2(ic,:))),...
              'cannot append objects that have different values at common epochs')
          end
        end
      end
      %update object
      obj=obj1.assign(y_now,'x',x_now,'mask',mask_now);
    end
    function obj1_out=augment(obj1,obj2,varargin)
      %NOTICE:
      % - obj1 receives the data from obj2, at those epochs defined in obj2
      % - no explicit gaps are ever copied (i.e. obj1 valid data is not turned to gaps)
      % - 'common' copies the common valid entries in obj2 to obj1 (default true)
      % - 'new' copies the new entries in obj2 to obj1 (default true)
      % - 'old' makes sure the existing data is not deleted (basically over-rides
      %   the previous two args, default false)
      % - 'skip_gaps' does not copy explicit gaps from obj2 to obj1 (only valid data)
      % - notice that common=true,new=*,old=true is the same as common=false,new=*,old=true: old data
      %   is not deleted in both cases, so the common data is never copied anyway
      % Parse inputs
      p=machinery.inputParser;
      % optional arguments
      p.addParameter('quiet',       false, @(i)islogical(i) && isscalar(i));
      p.addParameter('common',      true,  @(i)islogical(i) && isscalar(i));
      p.addParameter('new',         true,  @(i)islogical(i) && isscalar(i));
      p.addParameter('old',         false, @(i)islogical(i) && isscalar(i));
      p.addParameter('skip_gaps',   false, @(i)islogical(i) && isscalar(i));
      % parse it
      p.parse(varargin{:});
      %need compatible data
      obj1.compatible(obj2,varargin{:});
      %merge x domain, also get idx2, which is needed to propagate data
      if p.Results.skip_gaps
        [obj1_out,obj2_out,idx1,idx2]=merge(obj1,obj2.masked);
      else
        [obj1_out,obj2_out,idx1,idx2]=merge(obj1,obj2);
      end
      %get common data
      common_idx= (idx1==idx2) & obj1_out.mask & obj2_out.mask;
      %get new data
      new_idx=idx2 & obj2_out.mask;
      %get old data
      old_idx=idx1 & obj1_out.mask;
      %idx to be considered
      idx=false(size(idx1));
      %branch on method mode
      if p.Results.common;   idx=idx |     common_idx; end
      if p.Results.new;      idx=idx |        new_idx; end
      if p.Results.old;      idx=idx &       ~old_idx; end

%       [common_idx,new_idx,old_idx,idx]

      %count data in obj2 which has the same epochs as in obj1 but is different
      common_idx=find( (idx1==idx2) & idx);
      if any(common_idx) && ( any(obj1_out.mask(common_idx)) ||  any(obj1_out.mask(common_idx)) )
        common_diff_idx=common_idx(any(obj1_out.y(common_idx,:)~=obj2_out.y(common_idx,:),2));
        if ~isempty(common_diff_idx) && ~p.Results.quiet
          disp(['WARNING:BEGIN',10,'There are ',num2str(numel(common_diff_idx)),...
              ' entries in obj1 that are defined at the same epochs as in obj2 but with different values:'])
          for i=1:min([10,numel(common_diff_idx)])
            if isprop(obj1,'t')
              var_name='t';
            else
              var_name='x';
            end
            x1=obj1_out.tx( common_diff_idx(i));
            x2=obj2_out.tx( common_diff_idx(i));
            idx_str=num2str(common_diff_idx(i));
            y1=obj1_out.y(  common_diff_idx(i),:);
            y2=obj2_out.y(  common_diff_idx(i),:);
            if numel(y1)>8; y1=[y1(1:4), y1(end-5:end)]; end
            if numel(y2)>8; y2=[y2(1:4), y2(end-5:end)]; end
            disp(str.tablify([16,20,16,12],['obj1.',var_name,'(',idx_str,  ')='],x1,   ['obj1.y(',idx_str,',:)='],y1   ))
            disp(str.tablify([16,20,16,12],['obj2.',var_name,'(',idx_str,  ')='],x2,   ['obj2.y(',idx_str,',:)='],   y2))
            disp(str.tablify([16,20,16,12],['diff.',var_name,'(',idx_str,  ')='],x1-x2,['diff.y(',idx_str,',:)='],y1-y2))
          end
          disp(['The data of obj1 have been over-written by the values in obj2.',10,'WARNING:END'])
        end
      end
      if all(~idx)
        %do nothing
      else
        %retrieve original data at these entries, don't get gaps if requested
        if p.Results.skip_gaps
          obj2_new=obj2.masked.at(obj2_out.x(idx));
        else
          obj2_new=obj2.at(obj2_out.x(idx));
        end
        %propagate data
        obj1_out.y(idx,:)=obj2_new.y;
        %propagate mask
        obj1_out.mask(idx)=obj2_new.mask;
        obj1_out=obj1_out.mask_update;
      end
    end
    function out=isequal(obj1,obj2,columns)
      if ~exist('columns','var') || isempty(columns)
        columns=1:obj1.width;
      end
      if any(columns>obj1.width) || any(columns>obj2.width)
        error('requested column indices exceed width.')
      end
      %assume objects are not equal
      out=false;
      %check valid data length
      if obj1.nr_valid~=obj2.nr_valid
        return
      end
      %check values themselves
      if any(any(obj1.y_masked([],columns)~=obj2.y_masked([],columns)))
        return
      end
      %they are the same
      out=true;
    end
    function obj1=glue(obj1,obj2,varargin)
      %objects need to have the same time domain
      assert(obj1.isxsame(obj2),'Input objects do not share the same x-domain.')
      %make sure objects are compatible
      compatible(obj1,obj2,varargin{:})
      %augment the data, labels and units
      obj1=obj1.assign([obj1.y,obj2.y],'reset_width',true);
      obj1.labels=[obj1.labels(:);obj2.labels(:)]';
      obj1.units =[obj1.units(:); obj2.units(:) ]';
    end
    %% algebra (all operations assume the time domains are in-phase and no interpolation is done)
    function obj1=plus(obj1,obj2)
      %empty acts as zero
      if isempty(obj1)
        obj1=obj2;
      elseif isempty(obj2)
        %do nothing
      elseif isnumeric(obj2)
        obj1.y=obj1.y+obj2;
      elseif isnumeric(obj1)
        obj1=obj1+obj2.y;
      else
        %sanity
        compatible(obj1,obj2)
        %operate
        if obj1.length==1 && obj2.length==1
          assert(simpledata.isx('==',obj1.x,obj2.x,min([obj1.x_tol,obj2.x_tol])),...
            ['Cannot operate on scalar objects with different x-domains: ',...
            num2str(obj1.x),' != ',num2str(obj2.x)])
          obj1=obj1.assign_tx_mask(obj1.y+obj2.y,obj2.tx,obj1.mask&&obj2.mask);
        elseif obj1.length==1
          obj1=obj1.assign_tx_mask(ones(obj2.length,1)*obj1.y+obj2.y,obj2.tx,obj2.mask);
        elseif obj2.length==1
          obj1=obj1.assign_tx_mask(ones(obj1.length,1)*obj2.y+obj1.y,obj1.tx,obj1.mask);
        else
          %consolidate data sets
          [obj1,obj2]=obj1.merge(obj2,0);
          %operate
          obj1=obj1.assign(obj1.y+obj2.y,'mask',obj1.mask & obj2.mask);
        end
      end
    end
    function obj1=minus(obj1,obj2)
      %empty acts as zero
      if isempty(obj1)
        obj1=-obj2;
      elseif isempty(obj2)
        %do nothing
      elseif isnumeric(obj2)
        obj1.y=obj1.y-obj2;
      elseif isnumeric(obj1)
        obj1=obj1-obj2.y;
      else
        %sanity
        compatible(obj1,obj2)
        %operate
        if obj1.length==1 && obj2.length==1
          assert(simpledata.isx('==',obj1.x,obj2.x,min([obj1.x_tol,obj2.x_tol])), ...
            ['Cannot operate on scalar objects with different x-domains: ',...
            num2str(obj1.x),' != ',num2str(obj2.x)])
          if isprop(obj1,'epoch')
            assert(simpletimeseries.ist('==',obj1.epoch,obj2.epoch,min([obj1.t_tol,obj2.t_tol])), ...
              ['Cannot operate on scalar objects with different epochs: ',...
              datestr(obj1.epoch),' != ',datestr(obj2.epoch)])
          end
          obj1=obj1.assign_tx_mask(obj1.y-obj2.y,obj2.tx,obj1.mask&&obj2.mask);
        elseif obj1.length==1
          obj1=obj1.assign_tx_mask(ones(obj2.length,1)*obj1.y-obj2.y,obj2.tx,obj2.mask);
        elseif obj2.length==1
          obj1=obj1.assign_tx_mask(obj1.y-ones(obj1.length,1)*obj2.y,obj1.tx,obj1.mask);
        else
          %consolidate data sets
          [obj1,obj2]=obj1.merge(obj2,0);
          %operate
          obj1=obj1.assign(obj1.y-obj2.y,'mask',obj1.mask & obj2.mask);
        end
      end
    end
    function obj=scale(obj,scl)
      %empty acts as unit
      if isempty(scl)
        %do nothing
      elseif isscalar(scl)
        obj.y=scl*obj.y;
      elseif isvector(scl)
        if numel(scl)==obj.width
          for i=1:numel(scl)
            obj.y(:,i)=obj.y(:,i)*scl(i);
          end
        elseif numel(scl)==obj.length
          for i=1:obj.width
            obj.y(:,i)=obj.y(:,i).*scl(:);
          end
        else
          error(['if input <scale> is a vector, ',...
            'it must have the same length as either the witdh or length of <obj>.'])
        end
      else
        obj.y=obj.y.*scl;
      end
    end
    function obj=uminus(obj)
      obj=obj.scale(-1);
    end
    function obj=uplus(obj)
      obj=obj.scale(1);
    end
    function obj1=times(obj1,obj2) %this is .* (not *)
      %empty acts as unit
      if isempty(obj1)
        obj1=obj2;
      elseif isempty(obj2)
        %do nothing
      elseif isnumeric(obj2)
        obj1=obj1.scale(obj2);
      elseif isnumeric(obj1)
        obj1=obj2.scale(obj1).y;
      else
        %sanity
        compatible(obj1,obj2,'compatible_parameters',{'x_units'})
        if obj1.length==1
          obj1=obj1.assign_tx_mask(ones(obj2.length,1)*obj1.y.*obj2.y,obj2.tx,obj2.mask);
        elseif obj2.length==1
          obj1=obj1.assign_tx_mask(ones(obj1.length,1)*obj2.y.*obj1.y,obj1.tx,obj1.mask);
        else
          %consolidate data sets
          [obj1,obj2]=obj1.merge(obj2,1);
          %operate
          obj1=obj1.assign(obj1.y.*obj2.y,'mask',obj1.mask & obj2.mask);
        end
      end
    end
    function obj1=mtimes(obj1,obj2,map) %this is * (not .*)
      %shortcuts
      n=obj1.width;
      %use matlab's default column-first mapping
      if ~exist('map','var') || isempty(map)
        map=1:n*n;
      end
      %sanity
      assert(numel(map)==n*n,'Input ''map'' must have the same number of elements as obj2.width')
      assert(obj2.width==n*n,'Input ''obj2'' must represent a square matrix and have nr cols equal to obj1.width')
      %consolidate data sets
      [obj1,obj2]=obj1.merge(obj2,1);
      %get data
      y1=obj1.y;y2=obj2.y;
      %operate
      for i=1:obj1.length
        y1(i,:)=transpose(...
          reshape(y2(i,map),[n,n])*transpose(y1(i,:))...
        );
      end
      %save result
      obj1=obj1.assign(y1);
    end
    function obj1=rdivide(obj1,obj2)
      %empty is not supported
      assert(~isempty(obj1)&&~isempty(obj2),'Cannot handle empty inputs')
      if isnumeric(obj2)
        obj1=scale(obj1,1/obj2);
      elseif isnumeric(obj1)
        obj1=obj2.scale(1/obj1).y;
      else
        %sanity
        compatible(obj1,obj2,'compatible_parameters',{'x_units'})
        if obj1.length==1
          obj1=obj1.assign_tx_mask(ones(obj2.length,1)*obj1.y./obj2.y,obj2.tx,obj2.mask);
        elseif obj2.length==1
          obj1=obj1.assign_tx_mask(obj1.y./ones(obj1.length,1)*obj2.y,obj1.tx,obj1.mask);
        else
          %consolidate data sets
          [obj1,obj2]=obj1.merge(obj2,1);
          %operate
          obj1=obj1.assign(obj1.y./obj2.y,'mask',obj1.mask & obj2.mask);
        end
      end
    end
    function obj1=power(obj1,obj2)
      %empty is not supported
      assert(~isempty(obj1)&&~isempty(obj2),'Cannot handle empty inputs')
      if isnumeric(obj2)
        obj1.y=obj1.y.^obj2;
      elseif isnumeric(obj1)
        obj2.y=obj1.^obj2.y;
        obj1=obj2;
      else
        %operate
        if obj1.length==1
          obj1=obj1.assign_tx_mask(ones(obj2.length,1)*obj1.y.^obj2.y,obj2.tx,obj2.mask);
        elseif obj2.length==1
          obj1=obj1.assign_tx_mask(ones(obj1.length,1)*obj2.y.^obj1.y,obj1.tx,obj1.mask);
        else
          %consolidate data sets
          [obj1,obj2]=obj1.merge(obj2,1);
          %operate
          obj1=obj1.assign(obj1.y.^obj2.y,'mask',obj1.mask & obj2.mask);
        end
      end
    end
    %% logical ops
    function obj=and(obj,obj_new)
      %consolidate data sets
      [obj,obj_new]=obj.merge(obj_new);
      %operate
      obj=obj.mask_and(obj_new.mask);
    end
    function obj=or(obj,obj_new)
      %consolidate data sets
      [obj,obj_new]=obj.merge(obj_new);
      %operate
      obj=obj.mask_or(obj_new.mask);
    end
    function obj=not(obj)
      %build x domain of all explicit gaps
      x_now=obj.x(~obj.mask);
      %build zero data matrix
      y_now=zeros(numel(x_now),obj.width);
      %assign to output
      obj=obj.assign(y_now,'x',x_now);
    end
    %% directly-overloadable ops
    function obj=op(obj,operation)
      if ischar(operation)
        operation=str2func(operation);
      end
      assert(isa(operation,'function_handle'),['Input ''operation'' must be a function handle, not a ',class(operation),'.'])
      obj.y(obj.mask,:)=operation(obj.y_masked);
    end
    function obj=sqrt(   obj); obj=obj.op(@sqrt  );    end
    %cumulative ops
    function obj=cumsum( obj); obj=obj.op(@cumsum );    end
    function obj=cummax( obj); obj=obj.op(@cummax );    end
    function obj=cummin( obj); obj=obj.op(@cummin );    end
    function obj=cumprod(obj); obj=obj.op(@cumprod);    end
    %complex numbers
    function obj=abs(   obj); obj=obj.op(@abs   );    end
    function obj=conj(  obj); obj=obj.op(@conj  );    end
    function obj=imag(  obj); obj=obj.op(@imag  );    end
    function obj=real(  obj); obj=obj.op(@real  );    end
    function obj=sign(  obj); obj=obj.op(@sign  );    end
    %logarithms
    function obj=log(   obj); obj=obj.op(@log   );    end
    function obj=log10( obj); obj=obj.op(@log10 );    end
    function obj=log2(  obj); obj=obj.op(@log2  );    end
    %trigonometry
    function obj=sin(   obj); obj=obj.op(@sin   );    end
    function obj=cos(   obj); obj=obj.op(@cos   );    end
    function obj=tan(   obj); obj=obj.op(@tan   );    end
    function obj=cot(   obj); obj=obj.op(@cot   );    end
    function obj=sec(   obj); obj=obj.op(@sec   );    end
    function obj=csc(   obj); obj=obj.op(@csc   );    end
    function obj=asin(  obj); obj=obj.op(@asin  );    end
    function obj=acos(  obj); obj=obj.op(@acos  );    end
    function obj=atan(  obj); obj=obj.op(@atan  );    end
    function obj=acot(  obj); obj=obj.op(@acot  );    end
    function obj=asec(  obj); obj=obj.op(@asec  );    end
    function obj=acsc(  obj); obj=obj.op(@acsc  );    end
    function obj=sinh(  obj); obj=obj.op(@sinh  );    end
    function obj=cosh(  obj); obj=obj.op(@cosh  );    end
    function obj=tanh(  obj); obj=obj.op(@tanh  );    end
    function obj=coth(  obj); obj=obj.op(@coth  );    end
    function obj=sech(  obj); obj=obj.op(@sech  );    end
    function obj=csch(  obj); obj=obj.op(@csch  );    end
    function obj=asinh( obj); obj=obj.op(@asinh );    end
    function obj=acosh( obj); obj=obj.op(@acosh );    end
    function obj=atanh( obj); obj=obj.op(@atanh );    end
    function obj=acoth( obj); obj=obj.op(@acoth );    end
    function obj=asech( obj); obj=obj.op(@asech );    end
    function obj=acsch( obj); obj=obj.op(@acsch );    end
    %% decomposition/reconstruction
    function obj=parametric_decomposition(obj,varargin)
      %add some defaults
      v=varargs.wrap('sources',{....
        {...
          'epoch',    datetime('now'), @(i)isdatetime(i) || isscalar(i);...
          'timescale',      'seconds', @ischar;...
          't0',              obj.x(1), @num.isscalar;...
        },...
      },varargin{:});
      v=varargs.wrap('sources',{v,...
        {...
          'time', v.epoch+pardecomp.from_timescaled(obj.x,obj.timescale)...
        }...
      },varargin{:});
      obj=pardecomp.join(pardecomp.split(obj,v.varargin{:}),v.varargin{:});
    end
    %% differentiation
    function obj=diff(obj,varargin)
      p=machinery.inputParser;
      p.addParameter('mode','central', @ischar);
      % parse it
      p.parse(varargin{:});
      % retrieve x- and y-domain
      x_now=obj.x_masked;
      y_now=obj.y_masked;
      % branch on type of derivative
      switch p.Results.mode
      case 'forwards'
        % gather indexes
        iplus1=[false;true(obj.nr_valid-1,1)];
        iminus1=[true(obj.nr_valid-1,1);false];
        % compute forward derivative
        y_diff=[...
          (...
            y_now(iplus1,:)-y_now(iminus1,:)...
          )./(...
            (x_now(iplus1)-x_now(iminus1))*ones(1,obj.width)...
          );...
          nan(1,obj.width)...
        ];
      case 'backwards'
        % gather indexes
        iplus1=[false;true(obj.nr_valid-1,1)];
        iminus1=[true(obj.nr_valid-1,1);false];
        % compute forward derivative
        y_diff=[...
          nan(1,obj.width);...
          (...
            y_now(iplus1,:)-y_now(iminus1,:)...
          )./(...
            (x_now(iplus1)-x_now(iminus1))*ones(1,obj.width)...
          )...
        ];
      case 'central'
        % gather indexes
        iplus2=[false;false;true(obj.nr_valid-2,1)];
        iminus2=[true(obj.nr_valid-2,1);false;false];
        % compute central derivative
        y_diff=[...
          nan(1,obj.width);...
          (...
            y_now(iplus2,:)-y_now(iminus2,:)...
          )./(...
            (x_now(iplus2)-x_now(iminus2))*ones(1,obj.width)...
          );...
          nan(1,obj.width);...
        ];
      otherwise
        error(['unknown mode ''',p.Results.mode,'''.'])
      end
      %propagate data
      y_now=nan(obj.size);
      y_now(obj.mask,:)=y_diff;
      obj=obj.assign(y_now);
    end
    %% wrappers
    function obj=smooth(obj,varargin)
      x_now=obj.x_masked;
      y_now=obj.y_masked;
      mask_now=obj.mask;
      for i=1:obj.width
        obj=obj.set_cols(i,smooth(x_now,y_now(:,i),varargin{:}),mask_now);
      end
    end
    function obj=smooth_test(obj,span,varargin)
      modes={'moving','lowess','loess','sgolay','rlowess','rloess'};
      obj.plot
      for i=1:numel(modes)
        obj.smooth(span,modes{i},varargin{:}).plot
      end
      legend([{'original'},modes])
    end
    function obj=medfilt(obj,n)
      %NOTICE: this function does not decimate the data (as obj.median does), it simply applies
      %matlab's medfilt function to the data
      obj=obj.assign(medfilt1(obj.y,n)); %,'omitnan','truncate');
    end
    function obj=movmedian(obj,n)
      obj=obj.assign(movmedian(obj.y,n)); %,'omitnan','truncate');
    end
    %% vector
    function out=norm(obj,p)
      if ~exist('p','var') || isempty(p)
        p=2;
      end
      % computes the p-norm of y
      if p==2
        out=sqrt(sum(obj.y.^2,2));
      elseif p==1
        out=sum(abs(obj.y),2);
      elseif mod(p,2)==0
        out=sum(obj.y.^p,2).^(1/p);
      else
        out=sum(abs(obj.y).^p,2).^(1/p);
      end
    end
    function out=unit(obj)
      % computes the unit vector of y
      out=obj.y./(obj.norm*ones(1,obj.width));
    end
    function obj=dot(obj,obj_new)
      %sanity
      compatible(obj,obj_new,'compatible_parameters',{'x_units'})
      %consolidate data sets
      [obj,obj_new]=obj.merge(obj_new);
      %operate
      obj=obj.assign(sum(obj.y.*obj_new.y,2),'reset_width',true);
    end
    function obj=cross(obj,obj_new)
      %sanity
      compatible(obj,obj_new,'compatible_parameters',{'x_units'})
      %consolidate data sets
      [obj,obj_new]=obj.merge(obj_new);
      %operate
      obj=obj.assign(cross(obj.y,obj_new.y));
    end
    function obj=autocross(obj)
      %get autocross product
      ac=cross(obj.y(1:end-1,:),obj.y(2:end,:));
      %propagate
      obj.y=[...
        ac(1,:);...
        (ac(1:end-1,:)+ac(2:end,:))*5e-1;...
        ac(end,:)...
        ];
    end
    function obj=project(obj,x,y,z)
      %some sanity
      if obj.width~=3
        error('can only project 3D vectors')
      end
      %no need to consolidate, that is done inside dot
      px=obj.dot(x);
      py=obj.dot(y);
      pz=obj.dot(z);
      %maybe the x-domain needs to be updated
      if obj.length~=px.length || any(obj.x ~=px.x)
        obj.x=px.x;
      end
      %propagate projection
      obj.y=[px.y,py.y,pz.y];
    end
    %% calibration
    function obj1=calibrate_poly(obj1,obj2,order,varargin)
      % Determines coefficients of the polynomial expansion (y2i=obj2.y(:,i), y1i=obj1.y(:,i), k=order):
      %
      %  y2i=c0+c1*y1i+c2*y1i^2+...+ck*y1i^k
      %
      % that minimize:
      %
      %  [y2i - (c0+c1*y1i+c2*y1i^2+...+ck*y1i^k)]^2
      %
      % and applies them to y1i (i.e. the columns of the output obj1).
      if ~exist('order','var') || isempty(order)
        order=1;
      end
      %sanity
      compatible(obj1,obj2,varargin{:})
      %need to match the gaps
      [obj1,obj2]=mask_match(obj1,obj2);
      %go column-by column
      for i=1:obj1.width
        %easier names
        y1i=obj1.y_masked([],i);
        y2i=obj2.y_masked([],i);
        %build the design matrix
        A=cell2mat(arrayfun(@(o) y1i.^o,0:order,'UniformOutput',false));
        %solve coefficients
        c=transpose(A)*A\transpose(A)*y2i;
        %apply coefficients
        out=sum(cell2mat(arrayfun(@(o) c(o+1)*y1i.^o,0:order,'UniformOutput',false)),2);
        %propagate result
        obj1=obj1.set_cols(i,out,obj1.mask);
      end
    end
    %% plot methods
    function out=plot(obj,varargin)
      v=varargs.wrap('sources',{....
        {...
          'columns',      1:obj.width,@(i)isnumeric(i) || iscell(i);...
          'line'   ,      {},         @iscell;...
          'color',        {},         @iscell;...
          'zeromean',     false,      @(i)islogical(i) && isscalar(i);...
          'normalize',    false,      @(i)islogical(i) && isscalar(i);...
          'outlier_iter', 0,          @(i)isfinite(i)  && isscalar(i);...
          'title',        '',         @ischar;...
          'smooth_span',  0,          @(i)isfinite(i)  && isscalar(i) && i>=0;...
          'scale',        1,          @(i)all(isfinite(i));...
          'masked',       false,      @(i)islogical(i) && isscalar(i);...
        },...
      },varargin{:});
      % type conversions
      assert(~isempty(v.columns),'If input argument ''columns'' is given, it cannot be empty')
      if iscell(v.columns)
        columns=cell2mat(v.columns);
      else
        columns=v.columns;
      end
      % dimension conversions
      if numel(v.scale)==1
        scale=ones(size(v.columns))*v.scale;
      else
        str.sizetrap(v.scale,v.columns)
        scale=v.scale;
      end
      % plot
      for i=1:numel(columns)
        if columns(i)>obj.width
          disp(['WARNING: not plotting column ',num2str(columns(i)),' of ',obj.descriptor,...
            ' because width is limited to ',num2str(obj.width),'.'])
          continue
        end
        % get the data
        y_plot=obj.y(:,columns(i));
        % define a working mask
        out.mask{i}=true(size(obj.mask));
        % remove outliers if requested
        if v.outlier_iter>0
          for c=1:v.outlier_iter
             [y_plot,outlier_idx]=simpledata.rm_outliers(y_plot);
             %update the mask with the detected outliers
             out.mask{i}(outlier_idx)=false;
          end
        end
        % smooth if requested
        if v.smooth_span>0
          if isprop(obj,'t') && isduration(v.smooth_span)
            span=ceil(v.smooth_span/obj.step);
          else
            span=v.smooth_span;
          end
          y_plot(out.mask{i})=smooth(obj.x(out.mask{i}),y_plot(out.mask{i}),span,'moving');
          %update the mask with the to-be-deleted edges
          i1=min([           span,floor(obj.length/2)]);
          i2=max([obj.length-span, ceil(obj.length/2)]);
          out.mask{i}(1 :i1 )=false;
          out.mask{i}(i2:end)=false;
        end
        % skip showing explicit gaps
        if v.masked
          out.mask{i}=out.mask{i} & obj.mask;
        end
        %derive bias and amplitude
        if v.zeromean || v.normalize
          % get valid data
          y_valid=y_plot(out.mask{i});
          y_valid=y_valid(~isnan(y_valid));
        end
        % remove mean if requested
        if v.zeromean
          out.y_mean{i}=mean(y_valid);
        else
          out.y_mean{i}=0;
        end
        % normalize if requested
        if v.normalize
          out.y_scale{i}=diff(minmax(transpose(y_valid(:))));
        else
          out.y_scale{i}=1;
        end
        % do not plot zeros if requested
        if ~obj.plot_zeros
          out.mask{i}=out.mask{i} & ( obj.y(:,columns(i)) ~= 0 );
        end
        % get abcissae in datetime, if relevant
        x_plot=obj.tx_masked(out.mask{i});
        % apply mask to y data
        y_plot=y_plot(out.mask{i});
        % compute (de-meaned and/or normalized) ordinate
        y_plot=(y_plot-out.y_mean{i})/out.y_scale{i}*scale(i);
        % plot it
        if isempty(v.line)
          out.line_handle{i}=plot(x_plot,y_plot);hold on
        else
          if iscell(v.line{i})
            out.line_handle{i}=plot(x_plot,y_plot,v.line{i}{:});hold on
          else
            out.line_handle{i}=plot(x_plot,y_plot,v.line{i});hold on
          end
        end
      end
      % set the color, if given
      if ~isempty(v.color) && ~isempty(v.color{i})
        set(out.line_handle{i},'Color',v.color{i});
      end
      %set x axis
      if obj.length>1
        %get common axis limits (don't crop stuff)
        xl=plotting.common_lim(gca,'x');
        if isprop(obj,'t')
          try
            xlim(datetime(xl,'convertfrom','datenum'));
          catch
            xlim(xl);
          end
        else
          xlim(xl);
        end
      end
      %annotate
      if isempty(v.title)
        out.title=obj.descriptor;
      else
        out.title=v.title;
      end
      out.xlabel=['[',obj.x_units,']'];
      if numel(out.line_handle)==1
        out.ylabel=[obj.labels{columns},' [',obj.units{columns},']'];
        out.legend=obj.labels{columns};
        if v.zeromean
          out.ylabel=[out.ylabel,' ',num2str(out.y_mean{1},'%+.3g')];
          out.legend=[out.legend,' ',num2str(out.y_mean{1},'%+.3g')];
        end
        if v.normalize
          out.ylabel=[out.ylabel,' x ',num2str(out.y_scale{1},'%+.3g')];
          out.legend=[out.legend,' x ',num2str(out.y_scale{1},'%+.3g')];
        end
        out.legend={out.legend};
      else
        same_units=true;
        for i=2:numel(columns)
          if ~strcmp(obj.units{columns(1)},obj.units{columns(i)})
            same_units=false;
          end
        end
        if same_units
          out.ylabel=['[',obj.units{columns(1)},']'];
          out.legend=obj.labels(columns);
        else
          out.ylabel='';
          out.legend=arrayfun(@(i) strcat(obj.labels(i),' [',obj.units(i),']'),columns);
        end
        if v.zeromean
          for i=1:numel(columns)
            out.legend{i}=[out.legend{i},' ',num2str(out.y_mean{i})];
          end
        end
        if v.normalize
          for i=1:numel(columns)
            out.legend{i}=[out.legend{i},' x ',num2str(out.y_scale{i})];
          end
        end
      end
      %anotate
      if ~isempty(out.title);   title(str.clean(out.title, 'title')); end
      if ~isempty(out.xlabel); xlabel(str.clean(out.xlabel,'title')); end
      if ~isempty(out.ylabel); ylabel(str.clean(out.ylabel,'title')); end
      if ~isempty(out.legend); legend(str.clean(out.legend,'title')); end
      %special annotations
      if v.normalize; ylabel('[ ]'); end
      %maybe useful outside
      out.axis_handle=gca;
      out.fig_handle=gcf;
      %outputs
      if nargout == 0
        clear out
      end
    end
    %% convertion methods
    function out=xy(obj)
      %generating matrix
      out=[obj.x,obj.y];
    end
    function out=data(obj,varargin)
      % Parse inputs
      p=machinery.inputParser;
      % optional arguments
      p.addParameter('mode',  'xy',              @ischar);
      p.addParameter('masked',false,             @(i)islogical(i) && isscalar(i));
      p.addParameter('mask',  true(obj.length,1),@(i)islogical(i) && numel(i)==obj.length);
      % parse it
      p.parse(varargin{:});
      % branching on mode
      switch p.Results.mode
      case 'gaps'
        out=~(obj.mask & p.Results.mask);
      case 'mask'
        out=  obj.mask & p.Results.mask;
      otherwise
        out=obj.(p.Results.mode);
        if p.Results.masked
          out=out(obj.mask & p.Results.mask,:);
        else
          out=out(p.Results.mask,:);
        end
      end
    end
    %% export methods
    function out=pluck(obj,x_now,column)
      %retrieves one single data point, intends to be fast
      assert(isscalar(x_now) && isnumeric(x_now),'Input x_now must be a scalar number.')
      if ~exist('column','var') || isempty(column)
        column=1:obj.width;
      end
      out=obj.pluckq(x_now,column);
    end
    function out=pluckq(obj,x_now,column)
      %interpolate
      bounds=[find(obj.x<=x_now,1,'last');find(obj.x>=x_now,1,'first')];
      if numel(bounds)<2
        %out of bounds
        out=[];
      elseif bounds(1)==bounds(2)
        %data point exists, pluck it out
        out=obj.y(bounds(1),column);
      else
        %data point does not exist, interpolate
        out=interp1(obj.x(bounds),obj.y(bounds,column),x_now);
      end
    end
    function export(obj,filename,varargin)
      %determine file type
      [~,~,e]=fileparts(filename);
      %branch on extension
      switch e
      case '.mat'
        out=obj.data('masked',true);
        S.x=out(:,1);
        S.y=out(:,2:end);
        S.descriptor=obj.descriptor;
        S.units=obj.units;
        S.labels=obj.labels;
        S.x_units=obj.x_units;
        details=whos('S');
        disp(['Saving ',obj.descriptor,' to ',filename,' (',num2str(details.bytes/1024),'Kb).'])
        if nargin>1
          save(filename,'S',varargin{:})
        else
          save(filename,'S','-v7.3')
        end
      otherwise
        error(['cannot handle files of type ''',e,'''.'])
      end
    end
    %% import methods
%     function obj=import(filename)
%       %determine file type
%       [p,n,e]=fileparts(filename);
%       %branch on extension
%       switch e
%       otherwise
%         error(['cannot handle files of type ''',e,'''.'])
%       end
%     end
    end
 end