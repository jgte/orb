classdef grace
  %static
  properties(Constant)
    
%         'fields',struct(...
%           'AC0X'  ,1e-7,...
%           'AC0Z'  ,1e-6,...
%           'AC0Y1' ,1e-5,...
%           'AC0Y2' ,1e-5,...
%           'AC0Y3' ,1e-5,...
%           'AC0Y4' ,2e-5,...
%           'AC0Y5' ,2e-5,...
%           'AC0Y6' ,5e-5,...
%           'AC0Y7', 5e-5,...
%           'AC0Y8', 5e-5,...
%           'AC0XD', 5e-7,...
%           'AC0ZD', 5e-6,...
%           'AC0YD1',5e-5,...
%           'AC0YD2',5e-5,...
%           'AC0YD3',1e-4,...
%           'AC0YD4',1e-4,...
%           'AC0YD5',1e-4,...
%           'AC0YD6',5e-3,...
%           'AC0YD7',5e-3,...
%           'AC0YD8',5e-4,...
%           'AC0XQ' ,3e-7,...
%           'AC0ZQ' ,5e-6,...
%           'AC0YQ1',5e-3,...
%           'AC0YQ2',5e-3,...
%           'AC0YQ3',5e-3,...
%           'AC0YQ4',5e-3,...
%           'AC0YQ5',5e-3,...
%           'AC0YQ6',5e-3,...
%           'AC0YQ7',5e-3,...
%           'AC0YQ8',5e-4 ...
%         )...
        
    %default value of some internal parameters
    default_list=struct(...
      'par_modes',struct(...
        'stats',{{'mean','std','rms','length'}},...
        'stats2',{{'corrcoef','length'}}...
      ),...
      'par_calpar_csr',struct(...
        'longtermbias',struct(...
          'A',fullfile('input','bsA2003'),...
          'B',fullfile('input','bsB2003')...
        ),...
        'coord',{{'X','Y','Z'}},...
        'data_col_idx',1,...
        'levels',struct(...
          'aak',   1,...
          'accatt',1,...
          'estim', 1 ...
        ),...
        'fields',struct(...
          'AC0X'  ,1e-7,...
          'AC0XD', 5e-7,...
          'AC0XQ' ,3e-7,...
          'AC0Z'  ,1e-6,...
          'AC0ZD', 5e-6,...
          'AC0ZQ' ,5e-6...
        )...
      ),...
      'par_plot',struct(...
        'dir','plots',...
        'prefix','',...
        'suffix','',...
        'xlimits',[-inf,inf],...
        'ylimits',[-inf,inf],...
        'size',200+[0,0,21,9]*50,...
        'units','points',...
        'visible','on',...
        'fontsize',struct(...
          'axis', 24,...
          'title',32,...
          'label',28),...
        'line',struct(...
          'width',2)...
      )...
    );
    sats={'A','B'};
    parameter_list=struct(...
      'par_modes',      struct('default',grace.default_list.par_modes,     'validation',@(i) isstruct(i)),...
      'par_calpar_csr', struct('default',grace.default_list.par_calpar_csr,'validation',@(i) isstruct(i)),...
      'par_plot',       struct('default',grace.default_list.par_plot,      'validation',@(i) isstruct(i)),...
      'start',          struct('default',datetime([0 0 31]),               'validation',@(i) isdatetime(i)),...
      'stop',           struct('default',datetime([0 0 31]),               'validation',@(i) isdatetime(i))...
    );
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={};
  end
  %read only
  properties(SetAccess=private)
    par
    data
  end
  %private (visible only to this object)
  properties(GetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
    start
    stop
  end
  methods(Static)
    function out=parameters
      out=fieldnames(grace.parameter_list);
    end
  end
  methods
    %% constructor
    function obj=grace(varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      %declare parameters
      for j=1:numel(grace.parameters)
        %shorter names
        pn=grace.parameters{j};
        %declare parameters
        p.addParameter(pn,grace.parameter_list.(pn).default,grace.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(varargin{:});
      % reset data type list
      obj=obj.datatype_init;
      % save parameters with defaults first, they may be needed below
      for i=1:numel(grace.parameters)
        %shorter names
        pn=grace.parameters{i};
        %parameters have a dedicated member
        if strcmp(pn(1:4),'par_')
          obj.par.(pn(5:end))=collapse(p.Results.(pn));
        else
          obj.(pn)=collapse(p.Results.(pn));
        end
      end
      % data is added to this initialized object with the 'init' method (see below)
    end
    %% start/stop operations
    function out=get.start(obj)
      out=obj.sts('start',obj.cell_names);
      out=min([out{:}]);
    end
    function out=get.stop(obj)
      out=obj.sts('stop',obj.cell_names);
      out=max([out{:}]);
    end
    function obj=set.start(obj,start)
      %retrieve cell names
      names=obj.cell_names;
      %get cell array with simpletimeseries objects
      values=obj.cell_get(names);
      %loop over complete set of objects
      for i=1:numel(values)
        values{i}.start=start;
      end
      %back-propagate modified set of objects
      obj=obj.cell_set(names,values);
    end
    function obj=set.stop(obj,stop)
      %retrieve cell names
      names=obj.cell_names;
      %get cell array with simpletimeseries objects
      values=obj.cell_get(names);
      %loop over complete set of objects
      for i=1:numel(values)
        values{i}.stop=stop;
      end
      %back-propagate modified set of objects
      obj=obj.cell_set(names,values);
    end
    %% datatype operations
    function obj=datatype_init(obj)
      obj.data=struct([]);
    end
    function obj=datatype_set(obj,datatype_name,datatype_value)
      obj.data(1).(datatype_name)=datatype_value;
    end
    function out=datatype_get(obj,datatype_name)
      if isfield(obj.data,datatype_name)
        out=obj.data.(datatype_name);
      else
        out=[];
      end
    end
    function out=datatypes(obj)
      out=fieldnames(obj.data);
    end
    %% level operations
    function obj=level_init(obj,datatype)
      obj=obj.datatype_set(datatype,struct([]));
    end
    function obj=level_set(obj,datatype,level_name,level_value)
      datatype_value=obj.datatype_get(datatype);
      datatype_value.(level_name)=level_value;
      obj=obj.datatype_set(datatype,datatype_value);
    end
    function out=level_get(obj,datatype,level_name)
      datatype_value=obj.datatype_get(datatype);
      if isfield(datatype_value,level_name)
        out=datatype_value.(level_name);
      else
        out=[];
      end
    end
    function out=levels(obj,datatype)
      out=obj.datatype_get(datatype);
      if ~isempty(out)
        out=fieldnames(out);
      end
    end
    %% field operations
    function obj=field_init(obj,datatype,level)
      obj=obj.level_set(datatype,level,struct([]));
    end
    function obj=field_set(obj,datatype,level,field_name,field_value)
      level_value=obj.level_get(datatype,level);
      level_value.(field_name)=field_value;
      obj=obj.level_set(datatype,level,level_value);
    end
    function out=field_get(obj,datatype,level,field_name)
      level_value=obj.level_get(datatype,level);
      if isfield(level_value,field_name)
        out=level_value.(field_name);
      else
        out=[];
      end
    end
    function out=fields(obj,datatype,level)
      out=obj.level_get(datatype,level);
      if ~isempty(out)
        out=fieldnames(out);
      end
    end
    %% sat operations
    function obj=sat_init(obj,datatype,level,field)
      obj=obj.field_set(datatype,level,field,struct([]));
    end
    function obj=sat_set(obj,datatype,level,field,sat_name,sat_value)
      if all(cellfun(@isempty,strfind(grace.sats,sat_name)))
        error([mfilename,': sat_name can only be one of: ',strjoin(grace.sats,','),', not ''',sat_name,'''.'])
      end
      field_value=obj.field_get(datatype,level,field);
      field_value.(sat_name)=sat_value;
      obj=obj.field_set(datatype,level,field,field_value);
    end
    function out=sat_get(obj,datatype,level,field,sat_name)
      field_value=obj.field_get(datatype,level,field);
      if isfield(field_value,sat_name)
        out=field_value.(sat_name);
      else
        out=[];
      end
    end
    %% cell array operations
    function out=cell_names_sanity(~,in,name,default)
      if isempty(in)
        %assign default values
        out=default;
      else
        if ~iscellstr(in)
          if ischar(in)
            %if input is a string, convert it to cell
            out={in};
          else
            %cannot handle anything other than cell strings and strings
            error([mfilename,': input ''',name,''' must be of class ''cell'', not ''',class(in),'''.'])
          end
        else
          %we want cell strings
          out=in;
        end
      end
    end
    function out=cell_names_isvalid(obj,datatype,level,field,sat)
      try
        out=ischar(class(obj.data.(datatype).(level).(field).(sat)));
      catch
        out=false;
      end
    end
    function names=cell_names(obj,datatype,level,field,sat)
      %input arguments datatype, level, field and sat are either:
      % - a cell array of strings;
      % - a string (converted to a scalar cell array);
      % - empty (converted to all possible values of that type of argument).
      % This routine builds a list with all combinations of the given 
      % datatype, level, field and sat (after converting as mentioned above).
      if ~exist('sat',     'var'),      sat='';end
      if ~exist('field',   'var'),    field='';end
      if ~exist('level',   'var'),    level='';end
      if ~exist('datatype','var'), datatype='';end
      % some of datatype, level, field and sat do not depend on each other 
      % (those that do are inside the for loops below)
      datatype=obj.cell_names_sanity(datatype,'datatype',obj.datatypes);
           sat=obj.cell_names_sanity(sat,'sat',grace.sats);
      %make room for outputs (this is just a guess, because level and field can be empty)
      names=cell(numel(datatype)*numel(level)*numel(field)*numel(sat),1);
      c=0;
      %loop over all datatypes, levels, fields and sats
      for i=1:numel(datatype)
        level_now=obj.cell_names_sanity(level,'level',obj.levels(datatype{i}));
        for j=1:numel(level_now)
          field_now=obj.cell_names_sanity(field,'field',obj.fields(datatype{i},level_now{j}));
          for k=1:numel(field_now)
            for l=1:numel(sat)
              if obj.cell_names_isvalid(datatype{i},level_now{j},field_now{k},sat{l})
                c=c+1;
                names{c}={datatype{i},level_now{j},field_now{k},sat{l}};
              end
            end
          end
        end
      end
    end
    function values=cell_get(obj,names)
      %make room for outputs
      values=cell(size(names));
      %populate outputs
      for i=1:numel(values)
        values{i}=obj.sat_get(names{i}{:});
      end
    end
    function obj=cell_set(obj,names,values)
      %sanity
      if numel(names) ~= numel(values)
        error([mfilename,': number of elements of input ''values'' (',num2str(numel(values)),') ',...
          'must be the same as the number of data types (',num2str(numel(names)),').'])
      end
      %propagate inputs
      for i=1:numel(values)
        obj=obj.sat_set(names{i}{:},values{i});
      end
    end
    % Applies the function f to the cell array given by cell_get(datatype,level,field,sat) and
    % returns the result, so it transforms (tr) the object into something else
    function out=tr(obj,f,names,varargin)
      %collapse requested data into cell array and operate on the cell array
      out=f(obj.cell_get(names),varargin{:});
    end
    % Applies the function f to the cell array given by cell_get(datatype,level,field,sat) and
    % saves the resulting objects back to the same set, given by cell_set(datatype,level,field,sat)
    function obj=op(obj,f,names,varargin)
      %operate on the requested set
      values=obj.tr(f,names,varargin{:});
      %sanity
      if ~iscell(values)
        error([mfilename,': function ',fun2str(f),' must return a cell array, not a ',class(values),'.'])
      end
      %propagate cell array back to object
      obj=obj.cell_set(names,values);
    end
    % Retrieves the values of particular field or zero-input method from the 
    % simpletimeseries (sts) stored at the leaves defined by the input names.
    function out=sts(obj,sts_field,names)
      out=cell(size(names));
      obj_list=obj.cell_get(names);
      for i=1:numel(out)
        out{i}=obj_list{i}.(sts_field);
      end
    end
    %% utils
    function enforce_plot_par(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('fig_handle',  gcf,  @(i) ishandle(i));
      p.addParameter('axis_handle', gca,  @(i) ishandle(i));
      p.addParameter('xdate',       true, @(i) islogical(i));
      p.addParameter('ylimits',     obj.par.plot.ylimits, @(i) numel(i)==2 && isnumeric(i));
      p.addParameter('autoscale',   false, @(i)islogical(i) && isscalar(i));
      p.addParameter('automean',    false, @(i)islogical(i) && isscalar(i));
      p.addParameter('outlier',0,        @(i) isfinite(i));
      % parse it
      p.parse(varargin{:});
      % sanity
      if any(isfinite(p.Results.ylimits)) && ( p.Results.autoscale ||  p.Results.automean )
        error([mfilename,': option ''ylimits'' and ''autoscale'' or ''automean'' do not work concurrently.'])
      end
      
      % enforce fontsize and paper size
      set(    p.Results.axis_handle,          'FontSize',obj.par.plot.fontsize.axis)
      set(get(p.Results.axis_handle,'Title' ),'FontSize',obj.par.plot.fontsize.title);
      set(get(p.Results.axis_handle,'XLabel'),'FontSize',obj.par.plot.fontsize.label);
      set(get(p.Results.axis_handle,'YLabel'),'FontSize',obj.par.plot.fontsize.label);
      set(    p.Results.fig_handle, 'Position',          obj.par.plot.size,...
                                    'PaperUnits',        obj.par.plot.units,...
                                    'PaperPosition',     obj.par.plot.size);
      % enforce line properties
      line_handles=get(p.Results.axis_handle,'Children');
      for i=1:numel(line_handles)
        set(line_handles(i),'LineWidth',obj.par.plot.line.width)
      end
      % enforce x-limits
      v=axis(p.Results.axis_handle);
      if p.Results.xdate
        for i=1:2
          if isfinite(   obj.par.plot.xlimits(i))
            v(i)=datenum(obj.par.plot.xlimits(i));
          end
        end
      end
      % enforce reasonable y-limits
      for i=1:2
        if isfinite(p.Results.ylimits(i))
          v(i+2)=   p.Results.ylimits(i);
        end
      end
      % enforce auto-scale and/or auto-mean
      if p.Results.automean || p.Results.autoscale
        %gather plotted data
        dat=cell(size(line_handles));
        for i=1:numel(line_handles)
          dat{i}=get(line_handles(i),'ydata');
        end
        dat=[dat{:}];
        dat=dat(~isnan(dat(:)));
        %remove outliers (if requested)
        for c=1:p.Results.outlier
          dat=simpledata.rm_outliers(dat);
        end
        dat=dat(~isnan(dat(:)));
        %enfore data-driven mean and/or scale
        if p.Results.automean && p.Results.autoscale
          v(3:4)=mean(dat)+3*std(dat)*[-1,1];
        elseif p.Results.automean
          v(3:4)=mean(dat)+0.5*diff(v(3:4))*[-1,1];
        elseif p.Results.autoscale
          v(3:4)=mean(v(3:4))+3*std(dat)*[-1,1];
        end
        %fix non-negative data sets
        if all(dat>=0) && v(3)<0
          v(3)=0;
        end
      end
      axis(p.Results.axis_handle,v);
    end
    function out=id(obj,mode,datatype,level,field,sat,suffix)
    %Produces unique identifiers, useful for plot/data/... file names
      if ~exist('suffix','var') || isempty(suffix)
        suffix='';
      end
      switch mode
      case 'plot'
        parts={obj.par.plot.dir,    '/';...
               obj.par.plot.prefix, '.';...
               datatype,            '.';...
               level,               '.';...
               field,               '.';...
               sat,                 '.';...
               suffix,              '.';...
               obj.par.plot.suffix, '.'};
        out=cell(numel(parts)/2,1);
        for i=1:numel(parts)/2
          if ~isempty(parts{i,1})
            out{i}=parts{i,1};
            if ~strcmp(out{i}(end),parts{i,2})
              out{i}=[out{i},parts{i,2}];
            end
          else
            out{i}='';
          end
        end
        out=[strjoin(out,''),'png'];
      case 'mat'
        parts={obj.par.plot.dir,    '/';...
               datatype,            '.';...
               level,               '.';...
               field,               '.';...
               sat,                 '.';...
               suffix,              '.'};
        out=cell(numel(parts)/2,1);
        for i=1:numel(parts)/2
          if ~isempty(parts{i,1})
            out{i}=parts{i,1};
            if ~strcmp(out{i}(end),parts{i,2})
              out{i}=[out{i},parts{i,2}];
            end
          else
            out{i}='';
          end
        end
        out=[strjoin(out,''),'mat'];
      otherwise
        error([mfilename,': cannot handle mode ',mode,'.'])
      end
      if isempty(dir(fileparts(out)))
        system(['mkdir -vp ',fileparts(out)]);
      end
    end
    %% datatype initialization
    function obj=init(obj,datatype,storage,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('start',    obj.start,@(i) isdatetime(i));
      p.addParameter('stop',     obj.stop, @(i) isdatetime(i));
      p.addParameter('datafile', fullfile(storage,[datatype,'.mat']),@(i) ischar(i));
      p.addParameter('jobid_col',2,       @(i) isscalar(i) && isfinite(i));
      p.addParameter('release', '61',     @(i) ischar(i));
      p.addParameter('level',   'unknown',@(i) ischar(i));
      p.addParameter('field',   'unknown',@(i) ischar(i));
      % parse it
      p.parse(varargin{:});
      % branch according to datatype
      switch datatype
      case 'calpar_csr'
        %NOTICE: the following parameters are relevant to this datatype:
        % - datafile
        % - jobid_col
        %Everything else is ignored quietly.
        
        %check if data is already in matlab format
        if isempty(dir(p.Results.datafile))
          %get names of parameters and levels
          fields=fieldnames(obj.par.calpar_csr.fields);
          levels=fieldnames(obj.par.calpar_csr.levels);
          %need to get long-term biases
          for sat=grace.sats
            ltb.(sat{1})=flipud(transpose(...
              dlmread(['/Users/teixeira/data/csr/corral-tacc/input/bs',sat{1},'2003'])...
            ));
          end
          %load data
          for i=1:size(levels)
            for j=1:size(fields)
              tmp=struct('A',[],'B',[]);
              %read data
              for sat=grace.sats
                f=fullfile(storage,['gr',sat{1},'.',fields{j},'.',levels{i},'.GraceAccCal']);
                tmp.(sat{1})=simpletimeseries.import(f,'cut24hrs',false);
                %enforce the long-term biases
                switch fields{j}
                case 'AC0X'
                  lbt_idx=2;
                case {'AC0Y1','AC0Y2','AC0Y3','AC0Y4','AC0Y5','AC0Y6','AC0Y7','AC0Y8'}
                  lbt_idx=3;
                case 'AC0Z'
                  lbt_idx=4;
                otherwise
                  lbt_idx=0;
                end
                if lbt_idx>0
                  t=tmp.(sat{1}).mjd-ltb.(sat{1})(2,1);
                  tmp.(sat{1})=tmp.(sat{1}).assign(...
                    [tmp.(sat{1}).y(:,1)+polyval(ltb.(sat{1})(:,lbt_idx),t),tmp.(sat{1}).y(:,2:end)]...
                  );
                end
              end
              %ensure date is compatible between the satellites
              if ~tmp.A.isteq(tmp.B)
                [tmp.A,tmp.B]=tmp.A.merge(tmp.B);
              end
              %propagate data to object
              for sat=grace.sats
                obj=obj.sat_set('calpar_csr',levels{i},fields{j},sat{1},tmp.(sat{1}));
              end
            end
          end
          %loop over all sat and level
          for i=1:size(levels)
            for sat=grace.sats
              %gather names for this sat and level
              names=obj.cell_names('calpar_csr',levels{i},'',sat{1});
              %ensure the time domain is the same for all fields (in each sat and level)
              obj=obj.op(@simpledata.merge_multiple,...
                names,...
                ['GRACE-',sat{1},' ',levels{i},' calpar_csr']...
              );
              %ensure job IDs are consistent for all fields (in each sat and level)
              equal_idx=obj.tr(@simpledata.isequal_multiple,...
                names,...
                p.Results.jobid_col,['GRACE-',sat{1},' ',levels{i},' calpar_csr']...
              );
              for j=1:numel(equal_idx)
                if ~equal_idx{j}
                  error([mfilename,': Job ID inconsistency between ',...
                    '[',strjoin(names{j  },','),'] and ',...
                    '[',strjoin(names{j+1},','),']  (possibly more).'])
                end
              end
            end
          end
          %save data
          s=obj.datatype_get('calpar_csr'); %#ok<*NASGU>
          save(p.Results.datafile,'s');
          clear s
        else
          %load data
          load(p.Results.datafile,'s');
          levels=fieldnames(s); %#ok<NODEF>
          for i=1:numel(levels)
            obj=obj.level_set('calpar_csr',levels{i},s.(levels{i}));
          end
        end
      case 'acc'
        %NOTICE: the following parameters are relevant to this datatype:
        % - start
        % - stop
        % - release
        % - level
        % - field
        %Everything else is ignored quietly.

        % sanity
        if isempty(p.Results.start) || isempty(p.Results.stop)
          error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
        end
        % gather list of daily dates
        datelist=simpletimeseries.list(p.Results.start,p.Results.stop,days(1));
        for sat={'A','B'}
          switch [p.Results.level,'-',p.Results.field]
          case 'l1b-csr'
            %build required file list
            filelist=cell(size(datelist));
            for i=1:numel(datelist)
              filelist{i}=fullfile(storage,datestr(datelist(i),'yy'),'acc','asc',...
                ['ACC1B_',datestr(datelist(i),'yyyy-mm-dd'),'_',sat{1},'_',p.Results.release,'.asc']...
              );
            end
            %load and save the data
            obj=obj.sat_set('acc',p.Results.level,p.Results.field,sat{1},...
              simpletimeseries.import(filelist,'cut24hrs',false)...
            );
          case 'mod-nrtdm'
            %load the data
            tmp=nrtdm(['G',sat{1},'_Panels/Aero'],p.Results.start,p.Results.stop,'data_dir',storage);
            %save the data
            obj=obj.sat_set('acc',p.Results.level,p.Results.field,sat{1},...
              tmp.ts...
            );
          case 'mod-csr'
            %build required file list
            filelist=cell(size(datelist));
            for i=1:numel(datelist)
              filelist{i}=fullfile(storage,datestr(datelist(i),'yy'),datestr(datelist(i),'mm'),'gps_orb_l',...
                ['grc',sat{1},'_gps_orb_',datestr(datelist(i),'yyyy-mm-dd'),...
                '_RL',p.Results.release,'_GPSRL62_RL05.*.acc']...
              );
            end
            %load and save the data
            obj=obj.sat_set('acc',p.Results.level,p.Results.field,sat{1},...
              simpletimeseries.import(filelist,'cut24hrs',false)...
            );
          end
        end
      otherwise
        error([mfilename,': cannot handle datatype ''',datatype,'''.'])
      end
      %make sure start/stop options are honoured
      obj.start=p.Results.start;
      obj.stop= p.Results.stop;
    end
    %% costumized operations
    function out=data_op(obj,opname,datatype,level,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('outlier',  0, @(i) isfinite(i));
      % parse it
      p.parse(varargin{:});
      switch lower(datatype)
      case 'calpar_csr'
        switch lower(opname)
        case 'load-stats'
          stats_filename=obj.id('mat','calpar_csr',level,'','','stats');
          if isempty(dir(stats_filename))
            for sat={'A','B'};
              %retrive data for this satellite
              [dv,dn]=obj.cell_get('calpar_csr',level,'',sat{1});
              % loop over CSR cal pars
              fields=obj.fields('calpar_csr',level);
              for i=1:size(fields)
                %sanity
                if ~strcmp(fields{i},dn{i})
                  error([mfilename,': discrepancy in the field names, debug needed!'])
                end
                if i>1 && ~dv{1}.isteq(dv{i})
                  error([mfilename,': time domain discrepancy in ',fields{i}])
                end
                %compute stats per period
                s.(fields{i}).(sat{1})=dv{i}.stats(...
                  'period',years(1)/12,...
                  'overlap',seconds(0),...
                  'outlier',p.Results.outlier,...
                  'struct_fields',obj.par.modes.stats...
                );
              end
            end
            save(stats_filename,'s')
          else
            load(stats_filename)
          end
          %outputs
          out=s;
        case 'load-stats2'
          stats_filename=obj.id('mat','calpar_csr',level,'','','stats2');
          if isempty(dir(stats_filename))
            %retrive data
            A=obj.cell_get('calpar_csr',level,'','A');
            B=obj.cell_get('calpar_csr',level,'','B');
            % loop over CSR cal pars
            fields=obj.fields('calpar_csr',level);
            for i=1:size(fields)
              %compute stats per period
              s.(fields{i})=simpletimeseries.stats2(...
                A{i},...
                B{i},...
                'period',years(1)/12,...
                'overlap',seconds(0),...
                'outlier',p.Results.outlier,...
                'struct_fields',obj.par.modes.stats2...
              );
            end
            save(stats_filename,'s')
          else
            load(stats_filename)
          end
          %outputs
          out=s;
        case 'time-sanity'
          if isempty(level)
            level={'aak','accatt','estim'};
          elseif ~iscell(level)
            level={level};
          end
          %set outputs (it never gets there is there's any failed sanity check
          out=true;
          %loop over all sats
          for sat=grace.sats
            %loop over all required levels
            for i=1:numel(level)
              fields=obj.fields(datatype,level{i});
              switch level{i}
              case 'estim'
                %this check ensures that there are no arcs with lenghts longer than consecutive time stamps
                for j=1:numel(fields)
                  disp(['Checking ',datatype,'-',level{i},'-',fields{j},'-',sat{1}])
                  %save time series into dedicated var
                  ts_now=obj.sat_get(datatype,level{i},fields{j},sat{1});
                  %forget about epochs that have been artificially inserted to represent forward steps
                  idx1=find(diff(ts_now.t)>seconds(1));
                  %get arc lenths
                  al=ts_now.y(idx1,3);
                  %get consecutive time difference
                  dt=seconds(diff(ts_now.t(idx1)));
                  %find arcs that span over time stamps
                  bad_idx=find(al(1:end-1)-dt>ts_now.t_tol); %no abs here!
                  %report if any such epochs have been found
                  if ~isempty(bad_idx)
                    str=cell(1,min([numel(bad_idx),10]));
                    for k=1:numel(str)
                      idx=idx1(bad_idx(k));
                      str{k}=[...
                        datestr(ts_now.t(idx)),' ',...
                        num2str(al(bad_idx(k)),'%.5d'),' ',...
                        num2str(dt(bad_idx(k))),' ',...
                        num2str(al(bad_idx(k))-dt(bad_idx(k)))...
                        ];
                    end
                    error([mfilename,...
                      ': found ',num2str(sum(bad_idx)),' arc lengths (2nd column) longer than ',...
                      ' difference between consecutive time stamps (3rd column):',10,...
                      strjoin(str,'\n')...
                    ])
                  end
                end
              case {'aak','accatt'}
                %this check ensures that the t0 value is the same as the start of the arc
                for j=1:numel(fields)
                  %some fields do not have t0
                  if ~any(fields{j}(end)=='DQ')
                    disp(['Skipping ',datatype,'-',level{i},'-',fields{j},'-',sat{1}])
                    continue
                  end
                  disp(['Checking ',datatype,'-',level{i},'-',fields{j},'-',sat{1}])
                  %save time series into dedicated var
                  ts_now=obj.sat_get(datatype,level{i},fields{j},sat{1});
                  %forget about epochs that have been artificially inserted to represent forward steps
                  idx1=find(diff(ts_now.t)>seconds(1));
                  %get t0
                  t0=simpletimeseries.utc2gps(datetime(ts_now.y(idx1,3),'convertfrom','modifiedjuliandate'));
                  %find arcs that have (muhc) t0 different than their first epoch
                  bad_idx=find(abs(ts_now.t(idx1)-t0)>seconds(1));
                  %report if any such epochs have been found
                  if ~isempty(bad_idx)
                    str=cell(1,min([numel(bad_idx),10]));
                    for k=1:numel(str)
                      idx=idx1(bad_idx(k));
                      str{k}=[...
                        datestr(ts_now.t(idx1(bad_idx(k)))),' ',...
                        datestr(t0(bad_idx(k))),' ',...
                        char(seconds(ts_now.t(idx1(bad_idx(k))-t0(bad_idx(k)))))...
                      ];
                    end
                    error([mfilename,...
                      ': found ',num2str(sum(bad_idx)),' arc lengths (2nd column) longer than ',...
                      ' difference between consecutive time stamps (3rd column):',10,...
                      strjoin(str,'\n')...
                    ])
                  end
                end
              end
            end
          end
        case 'calmod'
          dci=obj.par.calpar_csr.data_col_idx;
          coord=obj.par.calpar_csr.coord;
          for sat=grace.sats
            %gather quantities
            acc  =obj.sat_get('acc','l1b','csr',sat{1});
            %init models container
            calmod=simpletimeseries(acc.t,zeros(acc.length,numel(coord)));
            for i=1:numel(coord)
              %skip Y coordinate for now
              if coord{i}=='Y',continue,end
              %build nice structure with the relevant calibration parameters
              cal=struct(...
                'ac0' ,obj.sat_get(datatype,level,['AC0',coord{i}    ],sat{1}).interp(acc.t),...
                'ac0d',obj.sat_get(datatype,level,['AC0',coord{i},'D'],sat{1}).interp(acc.t),...
                'ac0q',obj.sat_get(datatype,level,['AC0',coord{i},'Q'],sat{1}).interp(acc.t)...
              );
              %sanity
              if any([isempty(acc),isempty(cal.ac0),isempty(cal.ac0d),isempty(cal.ac0q)])
                error([mfilename,': not all data is available to perform this operation.'])
              end
              %retrieve time domain (it is the same for all cal pars)
              fields=fieldnames(cal);
              for l=1:numel(fields)
                t.(fields{l})=days(acc.t-simpletimeseries.ToDateTime(cal.(fields{l}).y(:,end),'modifiedjuliandate'));
              end
              %paranoid sanity check
              good_idx=~isnan(t.ac0);
              if any(t.ac0(good_idx)~=t.ac0d(good_idx)) || any(t.ac0(good_idx)~=t.ac0q(good_idx))
                error([mfilename,': calibration time domain inconsistent between parameters, debug needed!'])
              end
              %build calibration model
              calmod=calmod.set_cols(i,...
                cal.ac0.cols(dci)+...
                cal.ac0d.cols(dci).times(t.ac0d)+...
                cal.ac0q.cols(dci).times(t.ac0q.^2)...
              );
            end
            %save it
            obj=obj.sat_set('acc','calmod',level,sat{1},calmod);
          end
          %outputs
          out=obj;
        otherwise
          error([mfilename,': cannot handle operation ''',opname,''' for datatype ''',datatype,'''.'])
        end
      case 'acc'
        switch lower(opname)
        case 'load-stats'
          error([mfilename,': implementation needed.'])
        case 'load-stats2'
          error([mfilename,': implementation needed.'])
        otherwise
          error([mfilename,': cannot handle operation ''',opname,''' for datatype ''',datatype,'''.'])
        end
      otherwise
        error([mfilename,': cannot handle datatype ''',datatype,'''.'])
      end
    end
    %% costumized plotting routine
    function plot(obj,plotname,datatype,level,varargin)
      switch lower(datatype)
      case 'calpar_csr'
        switch lower(plotname)
        case ''
          % loop over CSR cal pars
          fields=obj.fields(datatype,level);
          for j=1:size(fields)
            filename=obj.id('plot',datatype,level,fields{j},'');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %retrive data
              sat=struct(...
                'A',obj.sat_get(datatype,level,fields{j},'A'),...
                'B',obj.sat_get(datatype,level,fields{j},'B')...
              );
              %plot data
              sat.A.plot('columns',obj.par.calpar_csr.data_col_idx)
              sat.B.plot('columns',obj.par.calpar_csr.data_col_idx)
              %legend
              legend({'GRACE-A','GRACE-B'});
              obj.enforce_plot_par('autoscale',true,'automean',true,'outlier',3)
              grid on
              ylabel(sat.A.y_units{1}) %could also be sat.B.y_units{1}
              title([fields{j},' ',level])
              saveas(gcf,filename)
            end
          end
        case 'stats'
          s=obj.data_op('load-stats',datatype,level,varargin{:});
          fields=fieldnames(obj.par.calpar_csr.fields);
          %plot number of data
          obj.par.plot.prefix='length.bars';
          for i=1:numel(fields)
            filename=obj.id('plot',datatype,level,fields{i},'');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              b=bar(...
                mean([...
                  datenum(s.(fields{i}).A.length.t),...
                  datenum(s.(fields{i}).B.length.t)...
                ],2),[...
                  s.(fields{i}).A.length.scale(0.5).y(:,1),...
                  s.(fields{i}).B.length.scale(0.5).y(:,1)...
                ]...
              );
              b(1).EdgeColor = 'red';
              b(1).FaceColor = 'red';
              b(2).EdgeColor = 'black';
              b(2).FaceColor = 'black';
              %plot tweaking
              obj.enforce_plot_par
              legend('GRACE-A','GRACE-B')
              datetick('x','yyyy-mm')
              xlabel('time')
              grid on
              title(['Monthly count of ',fields{i},' ',level])
              saveas(gcf,filename)
            end
          end
          %plot more stats
          for mode=obj.par.modes.stats;
            %plot the requested statistics of the data
            obj.par.plot.prefix=mode{1};
            for i=1:numel(fields)
              filename=obj.id('plot',datatype,level,fields{i},'');
              if isempty(dir(filename))
                switch mode{1}
                case 'length'
                  continue
                  mask=true(s.(fields{i}).A.(mode{1}).length,1); %#ok<UNRCH>
                  scl=0.5;
                otherwise
                  mask=s.(fields{i}).A.length.y(:,1)>10;
                  scl=1;
                end
                sat=struct(...
                  'A',s.(fields{i}).A.(mode{1}).mask_and(mask).mask_update,...
                  'B',s.(fields{i}).B.(mode{1}).mask_and(mask).mask_update...
                );
                figure('visible',obj.par.plot.visible);
                sat.A.scale(scl).plot('columns',obj.par.calpar_csr.data_col_idx);
                sat.B.scale(scl).plot('columns',obj.par.calpar_csr.data_col_idx);
                %plot tweaking
                obj.enforce_plot_par('autoscale',true,'automean',true,'outlier',1)
                legend('GRACE-A','GRACE-B')
                ylabel(['[',sat.A.y_units{obj.par.calpar_csr.data_col_idx},']'])
                grid on
                title(['monthly ',mode{1},' of ',fields{i},' ',level])
                saveas(gcf,filename)
              end
            end
          end
        case 'corr'
          s=obj.data_op('load-stats2',datatype,level,varargin{:});
          fields=fieldnames(obj.par.calpar_csr.fields);
          % loop over all required modes
          for mode=obj.par.modes.stats2;
            %plot the requested statistics of the data
            obj.par.plot.prefix=mode{1};
            % loop over CSR cal pars
            for i=1:size(fields)
              filename=obj.id('plot',datatype,level,fields{i},'');
              if isempty(dir(filename))
                switch mode{1}
                case 'length'
                  continue
                  mask=true(s.(fields{i}).(mode{1}).length,1); %#ok<UNRCH>
                  ylimits=[-Inf,Inf];
                  scl=0.5;
                otherwise
                  mask=s.(fields{i}).length.y(:,1)>=10;
                  ylimits=[-1,1];
                  scl=1;
                end    
                figure('visible',obj.par.plot.visible);
                s.(fields{i}).(mode{1}).mask_and(mask).mask_update.scale(scl).plot(...
                  'columns',obj.par.calpar_csr.data_col_idx...
                );
                obj.enforce_plot_par('ylimits',ylimits)
                datetick('x','yyyy-mm')
                xlabel('time')
                ylabel('[ ]')
                grid on
                title(['monthly ',mode{1},' of ',fields{i},' ',level])
                saveas(gcf,filename)
              end
            end
          end
        case 'overview'
          fields={...
            'AC0X',...
            'AC0Y1',...
            'AC0Z',...
            'AC0XD',...
            'AC0YD1',...
            'AC0ZD',...
            'AC0XQ',...
            'AC0YQ1',...
            'AC0ZQ'...
          };
          %single-sat stats
          s=obj.data_op('load-stats',datatype,level,varargin{:});
          stats_modes=obj.par.modes.stats;
          for j=1:numel(stats_modes);
            %skip length, not so interesting in this mode
            if strcmp(stats_modes{j},'length')
              continue
            end
            obj.par.plot.prefix=stats_modes{j};
            filename=obj.id('plot',datatype,level,'','');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %gathering data
              y=nan(numel(fields),2);
              for i=1:numel(fields)
                if ~isfield(s,fields{i})
                  continue
                end
                y(i,1)=s.(fields{i}).A.(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.calpar_csr.data_col_idx);
                y(i,2)=s.(fields{i}).B.(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.calpar_csr.data_col_idx);
              end
              pos_idx=(y>0);
              %plot negative parameters downwards
              y_now=abs(y); y_now(pos_idx)=nan;
              b{1}=bar(log10(y_now)); hold on
              %plot positive parameters upwards
              y_now=abs(y); y_now(~pos_idx)=nan;
              b{2}=bar(-log10(y_now)); hold on
              %annotating
              obj.enforce_plot_par
              for k=1:2
                b{k}(1).EdgeColor = 'red';
                b{k}(1).FaceColor = 'red';
                b{k}(2).EdgeColor = 'black';
                b{k}(2).FaceColor = 'black';
              end
              if all(pos_idx(~isnan(y)))
                axis([0 numel(fields)+1 0 10])
              else
                axis([0 numel(fields)+1 -10 10])
              end
              set(gca,'XTick',1:numel(fields))
              set(gca,'XTickLabel',fields)
              set(gca,'XTickLabelRotation',45)
              yticks=get(gca,'YTickLabel');
              for k=1:numel(yticks)
                ytickval=str2double(yticks{k});
                if ytickval==0
                  yticks{k}='1';
                elseif ytickval>0
                  yticks{k}=['y.10^{-',yticks{k},'}'];
                elseif ytickval<0
                  yticks{k}=['-y.10^{',yticks{k},'}'];
                end
              end
              set(gca,'YTickLabel',yticks);
              legend('GRACE-A','GRACE-B')
              grid on
              title([...
                stats_modes{j},' of ',level,' for ',...
                datestr(obj.par.plot.xlimits(1)),' to ',...
                datestr(obj.par.plot.xlimits(2))...
              ]);
              saveas(gcf,filename)
            end
          end
          %all-sat stats
          s=obj.data_op('load-stats2',datatype,level,varargin{:});
          stats_modes=obj.par.modes.stats2;
          for j=1:numel(stats_modes);
            %skip length, not so interesting in this mode
            if strcmp(stats_modes{j},'length')
              continue
            end
            obj.par.plot.prefix=stats_modes{j};
            filename=obj.id('plot',datatype,level,'','');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %gathering data
              y=nan(numel(fields),1);
              for i=1:numel(fields)
                if ~isfield(s,fields{i})
                  continue
                end
                y(i,1)=s.(fields{i}).(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(data_col_idx);
              end
              bar(y,'k')
              %annotating
              obj.enforce_plot_par
              set(gca,'XTick',1:numel(fields))
              set(gca,'XTickLabel',fields)
              set(gca,'XTickLabelRotation',45)
              axis([0 numel(fields)+1 -1 1])
              grid on
              title([...
                stats_modes{j},' of ',level,' for ',...
                datestr(obj.par.plot.xlimits(1)),' to ',...
                datestr(obj.par.plot.xlimits(2))...
              ]);
              saveas(gcf,filename)
            end
          end  
%           case 'calpar_csr-tables'
%             fields=fieldnames(obj.par.calpar_csr.fields);
%             tab=12;
%             %single-sat stats
%             s=obj.data_op('load-stats',datatype,level,varargin{:});
%             stats_modes=obj.par.modes.stats;
%             for sat={'a','b'};
%               disp(['GRACE-',upper(sat{1})])
%               for i=0:numel(fields)
%                 out=cell(1,numel(stats_modes)+1);
%                 if i==0
%                   calpar_clean='cal par';
%                 else
%                   calpar_clean=str.clean(fields{i},'_');
%                 end
%                 out{1}=str.tabbed(calpar_clean,tab);
%                 for j=1:numel(stats_modes);
%                   if i==0
%                     out{j+1}=str.tabbed(stats_modes{j},tab,true);
%                   else
%                     d=s{i}.(sat{1}).(stats_modes{j}).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
%                     out{j+1}=str.tabbed(num2str(d(obj.par.calpar_csr.data_col_idx),'% .3g'),tab,true);
%                   end
%                 end
%                 disp(strjoin(out,' '))
%               end
%             end
%             %all-sat stats
%             s=obj.data_op('load-stats2',datatype,level,varargin{:});
%             stats_modes=obj.par.modes.stats2;
%             for i=0:numel(fields)
%               out=cell(1,numel(stats_modes)+1);
%               if i==0
%                 calpar_clean='cal par';
%               else
%                 calpar_clean=str.clean(fields{i},'_');
%               end
%               out{1}=str.tabbed(calpar_clean,tab);
%               for j=1:numel(stats_modes);
%                 if i==0
%                   out{j+1}=str.tabbed(stats_modes{j},tab,true);
%                 else
%                   d=s{i}.(stats_modes{j})(:,1).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
%                   out{j+1}=str.tabbed(num2str(mean(d(obj.par.calpar_csr.data_col_idx)),'% .3f'),tab,true);
%                 end
%               end
%               disp(strjoin(out,' '))
%             end
        otherwise
          error([mfilename,': cannot handle plotname ''',plotname,''' for datatype ''',datatype,'''.'])
        end
      case 'acc'
        switch lower(plotname)
        case ''
        otherwise
          error([mfilename,': cannot handle plotname ''',plotname,''' for datatype ''',datatype,'''.'])
        end
      otherwise
        error([mfilename,': unknown data of type ''',datatype,'''.'])
      end
    end
  end
end

function out=collapse(in)
  if ~isscalar(in) && ~iscell(in)
    %vectors are always lines (easier to handle strings)
    out=transpose(in(:));
  else
    out=in;
  end
end
