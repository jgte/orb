classdef gswarm
  methods(Static)
    function obj=load_models(obj,product,varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{...
        {...
          'import_dir',        '.',@(i) ischar(i) && exist(i, 'dir')~=0;...
          'model_format','unknonw',@(i) ischar(i);...
          'date_parser', 'unknonw',@(i) ischar(i);...
          'consistent_GM',   false,@(i) islogical(i);...
          'consistent_R',    false,@(i) islogical(i);...
          'max_degree'           0,@(i) isnumeric(i) && isscalar(i);...
          'use_GRACE_C20',   false,@(i) islogical(i);...
          'delete_C00',      false,@(i) islogical(i);...
          'delete_C20',      false,@(i) islogical(i);...
          'model_span', {time.zero_date,time.inf_date},@(i) isnumeric(i) && numel(i)==2;...
          'static_model',   'none',@(i) ischar(i);...
          'overwrite_common_t',false,@(i) islogical(i);...
        },...
        product.args...
      },varargin{:});
      %load all available data 
      [s,e]=gravity.load_dir(v.import_dir,v.model_format,str2func(v.date_parser),...
        'wilcarded_filename',v.wilcarded_filename,...
        'start',obj.start,...
        'stop',obj.stop,...
        'descriptor',product.name,...
        'overwrite_common_t',v.overwrite_common_t...
      );
      %apply model processing options
      [s,e]=gswarm.load_models_op('all',v,product,s,e);
      %make sure we got a gravity object
      assert(isa(s,'gravity'),['failed to load product ',product.codename])
      %propagate relevant data
      for i=1:numel(v.model_types)
        switch lower(v.model_types{i})
        case {'signal','sig','s'}
          obj=obj.data_set(product.dataname.append_field_leaf(v.model_types{i}),s);
        case {'error','err','e'}
          obj=obj.data_set(product.dataname.append_field_leaf(v.model_types{i}),e);
        otherwise
          error([mfilename,': unknown model type ''',v.model_types{i},'''.'])
        end
      end
    end
    function [mod,err]=load_models_op(mode,v,product,mod,err)
      if iscellstr(mode)
        for i=1:numel(mode)
          if exist('err','var')
            [mod,err]=gswarm.load_models_op(mode{i},v,product,mod,err);
          else
            mod=gswarm.load_models_op(mode{i},v,product,mod);
          end
        end
      else
        %sanity
        if exist('mod','var') && ~isa(mod,'gravity')
          str.say(['Ignoring model ',mod.descriptor,' since it is of class ''',class(mod),''' and not of class ''gravity''.'])
          return
        end
        %inits user feedback vars
        msg='';show_msg=false;
        switch mode
        case 'mode_list'
          mod={'consistent_GM','consistent_R','max_degree','use_GRACE_C20','delete_C00','delete_C20','start','stop','static_model'};
        case 'all'
          if exist('err','var')
            [mod,err]=gswarm.load_models_op(gswarm.load_models_op('mode_list'),v,product,mod,err);
          else
            mod=gswarm.load_models_op(gswarm.load_models_op('mode_list'),v,product,mod);
          end
        case 'consistent_GM'
          mod=mod.setGM(gravity.parameters('GM'));
          if exist('err','var')
            err=err.setGM(gravity.parameters('GM'));
          end
          show_msg=true;
        case 'consistent_R'
          mod=mod.setR( gravity.parameters('R'));
          if exist('err','var')
            err=err.setR( gravity.parameters('R'));
          end
          show_msg=true;
        case 'max_degree'
          %set maximum degree (if requested)
          if v.isparameter('max_degree') && v.max_degree>0
            mod.lmax=v.max_degree;
            if exist('err','var')
              err.lmax=v.max_degree;
            end
            msg=['set to ',num2str(v.max_degree)];
            show_msg=true;
          end
        case 'use_GRACE_C20'
          %set C20 coefficient
          if v.isparameter('use_GRACE_C20') && ~str.none(v.use_GRACE_C20)
            %some sanity
            if strcmpi(v.date_parser,'static')
              error([mfilename,': there''s no point in replacing GRACE C20 coefficients in a static model.'])
            end
            %legacy support
            if islogical(v.use_GRACE_C20) && v.use_GRACE_C20
              v.use_GRACE_C20='TN-07';
            end
            %get C20 timeseries, interpolated to current time domain
            c20=gravity.graceC20('version',v.use_GRACE_C20);
            %use parametric decomposition if outside of time domain
            if max(c20.t_masked)<max(mod.t_masked) || min(c20.t_masked)>min(mod.t_masked)
              switch  v.use_GRACE_C20
              case 'TN-11'
                %NOTICE: these periods were derived with:
                %[~,pd]=c20.parametric_decomposition_search('np',2,'T',[365.2426,182.6213],'timescale','days');pd.T
                T=[365.2426 182.6213 7936 992 566.85714 1984 91.218391 881.77778 721.45455 ...
                  214.48649 1587.2 529.06667 377.90476 122.09231 661.33333 273.65517 330.66667];
                np=2;
              case 'TN-07'
                error('needs implementation')
              end
              c20=c20.parametric_decomposition('np',np,'T',T,'timescale','days','time',mod.t);
            else
              c20=c20.interp(mod.t);
            end
  %         figure
  %         plot(c20.x_masked,c20.y_masked([],1),'x-','MarkerSize',10,'LineWidth',4), hold on
  %         plot(c20.x,spline(c20.x_masked,c20.y_masked([],1),c20.x),'o-','MarkerSize',10,'LineWidth',2)
  %         plot(c20.x,pchip(c20.x_masked,c20.y_masked([],1),c20.x),'+-','MarkerSize',10,'LineWidth',2)
  %         legend('data','spline','pchip')
            %interpolate in case there are gaps
            if c20.nr_gaps>0
              c20=c20.assign([...
                interp1(c20.x_masked,c20.y_masked([],1),c20.x),...
                interp1(c20.x_masked,c20.y_masked([],2),c20.x)...
              ],'t',c20.t,'mask',true(size(c20.x)));
            end
            mod=mod.setC(2,0,c20.y(:,1),mod.t);
            err=err.setC(2,0,c20.y(:,2),mod.t);
  %           figure
  %           plot(c20.t,c20.y(:,1),'o-'), hold on
  %           plot(m.t,m.C(2,0),'x-')
  %           legend('GRACE',m.descriptor)
            msg=['set to ',v.use_GRACE_C20];
            show_msg=true;
          end
        case 'delete_C00'
          %remove C00 bias
          if v.isparameter('delete_C00') && v.delete_C00
            mod=mod.setC(0,0,0);
            if exist('err','var')
              err=err.setC(0,0,0);
            end
            show_msg=true;
          end
        case 'delete_C20'
          %remove C20 bias
          if v.isparameter('delete_C20') &&  v.delete_C20
            mod=mod.setC(2,0,0);
            if exist('err','var')
              err=err.setC(2,0,0);
            end
            show_msg=true;
          end
        case 'model_span'
          error('no longer supported')
          msg={};
          %append extremeties, if requested
          if v.isparameter('model_span') && v.model_span{1}~=time.zero_date && v.model_span{1}<mod.start
            mod=mod.append(gravity.nan(mod.lmax,'t',v.model_span{1},'R',mod.R,'GM',mod.GM));
            if exist('err','var')
              err=err.append(gravity.nan(err.lmax,'t',v.model_span{1},'R',err.R,'GM',err.GM));
            end
            msg{end+1}=['from ',datestr(v.model_span{1})];
            show_msg=true;
          end
          if v.isparameter('model_span') && v.model_span{2}~=time.inf_date && v.model_span{2}>mod.stop
            mod=mod.append(gravity.nan(mod.lmax,'t',v.model_span{2},'R',mod.R,'GM',mod.GM));
            if exist('err','var')
              err=err.append(gravity.nan(err.lmax,'t',v.model_span{2},'R',err.R,'GM',err.GM));
            end
            msg{end+1}=['to ',datestr(v.model_span{2})];
            show_msg=true;
          end
          msg=strjoin(msg,' ');
        case 'start' %NOTICE: this is only done when loading the raw data (afterwards the matlab data is read directly, bypassing this routine altogher)
          if v.isparameter('start') && v.start~=time.zero_date
            if v.start<mod.start
              %append extremeties 
              mod=mod.append(gravity.nan(mod.lmax,'t',v.start,'R',mod.R,'GM',mod.GM));
              if exist('err','var')
                err=err.append(gravity.nan(err.lmax,'t',v.start,'R',err.R,'GM',err.GM));
              end
            elseif v.start>mod.start
              %trim extremeties (this is redundant unless data is saved before the start metadata is increased)
              mod=mod.trim(v.start,mod.stop);
              if exist('err','var')
                err=err.trim(v.start,mod.stop);
              end
            end
            msg=['at ',datestr(v.start)];
            show_msg=true;
          end
        case 'stop' %NOTICE: this is only done when loading the raw data (afterwards the matlab data is read directly, bypassing this routine altogher)
          if v.isparameter('stop') && v.stop~=time.inf_date
            if v.stop>mod.stop
              %append extremeties
              mod=mod.append(gravity.nan(mod.lmax,'t',v.stop,'R',mod.R,'GM',mod.GM));
              if exist('err','var')
                err=err.append(gravity.nan(err.lmax,'t',v.stop,'R',err.R,'GM',err.GM));
              end
            elseif v.stop<mod.stop
              %trim extremeties (this is redundant unless data is saved before the stop metadata is decreased)
              mod=mod.trim(mod.start,v.stop);
              if exist('err','var')
                err=err.trim(mod.start,v.stop);
              end
            end
            msg=['at ',datestr(v.stop)];
            show_msg=true;
          end
        case 'static_model'
          %remove static field (if requested)
          if v.isparameter('static_model') && ~strcmpi(v.static_model,'none')
            static=datastorage().init(v.static_model).data_get_scalar(datanames(v.static_model).append_field_leaf('signal'));
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
            show_msg=true;
          end
        otherwise
          error(['Cannot handle operantion ''',mode,'''.'])
        end
        if show_msg;str.say(product.str,':',mode,msg);end
      end
    end
    function obj=load_project_models(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{product.args},varargin{:});
      %loop over all requested scenarios
      for i=1:numel(v.scenarios)
        p=product;
        p.dataname=p.dataname.append_field_leaf(v.scenarios{i});
        obj=gswarm.load_models(obj,p,'import_dir',fullfile(v.import_dir,[v.project_name,'-',v.scenarios{i}]),varargin{:});
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=smooth_models(obj,product,varargin)
      %retrieve relevant parameters
      smoothing_degree  =product.mdget('smoothing_degree');
      smoothing_method  =product.mdget('smoothing_method');
      model_types       =product.mdget('model_types','always_cell_array',true);
      %sanity
      assert(product.nr_sources==1,['Can only handle one source model, not ',num2str(product.nr_sources),'.'])
      %loop over all models
      for i=1:numel(model_types)
        %gather model
        m=obj.data_get_scalar(product.sources(1).append_field_leaf(model_types{i}));
        %branch on the type of model, do not smooth error models
        switch lower(model_types{i})
        case {'error','err','e'}
          %do nothing
        otherwise
          if smoothing_degree>0
            %apply smoothing
            m=m.scale(smoothing_degree,smoothing_method);
          end
        end
        %save the smoothed model
        obj=obj.data_set(product.dataname.append_field_leaf(model_types{i}),m);
      end
    end
    function obj=combine_models(obj,product,varargin)
      %retrieve relevant parameters
      combination_type =product.mdget('combination_type');
      model_type       =product.dataname.field_path{end};
      %collect the models
      m=cell(product.nr_sources,1);
      for i=1:product.nr_sources
        m{i}=obj.data_get_scalar(product.sources(i));
      end
      %propagate relevant data
      switch lower(model_type)
      case {'error','err','e'}
        obj=obj.data_set(product,gravity.combine(m,'mode',combination_type,'type','error'));
      otherwise
        obj=obj.data_set(product,gravity.combine(m,'mode',combination_type,'type','signal'));
      end
    end
    function obj=stats(obj,product,varargin)
      %retrieve relevant parameters
      derived_quantity=product.mdget('plot_derived_quantity');
      functional      =product.mdget('functional');
      stats           =product.mdget('stats');
      stats_outlier   =product.mdget('stats_outlier');
      stats_period    =product.mdget('stats_period');
      stats_overlap   =product.mdget('stats_overlap');
      %sanity
      assert(product.nr_sources==1,['Can only handle one source model, not ',num2str(product.nr_sources),'.'])
      %get the model
      g=obj.data_get_scalar(product.sources(1)).scale(functional,'functional');
      %compute cumulative amplitude time series
      ts=g.derived(derived_quantity);
      %derive statistics
      s=ts.stats(...
        'struct_fields',stats,...
        'outlier',stats_outlier,...
        'period', stats_period,...
        'overlap',stats_overlap...
      );
      %add field with degree
      s.degree=0:g.lmax;
      %save data
      obj=obj.data_set(product,s);
    end
    %% plots
    function out=plot_ops(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_min_degree',        2 ,@(i) isnumeric(i) && isscalar(i);... %NOTICE: this needs to be handled externally!
          'plot_max_degree',      inf ,@(i) isnumeric(i) && isscalar(i);...
          'plot_smoothing_degree', [] ,@(i) isnumeric(i) && isscalar(i);...
          'plot_smoothing_method', '' ,@(i) ischar(i);...
          'plot_spatial_mask', 'none' ,@(i) iscellstr(i);...
          'plot_save_data',      'no' ,@(i) islogical(i) || ischar(i);...
          'stats_relative_to', 'none' ,@(i) ischar(i);...
          'model_types',   {'signal'} ,@(i) iscellstr(i);...
        },...
        product.args({'stats_relative_to','model_types'}),...
        product.plot_args...
      },varargin{:});
      %don't save data by default
      savedata=false;
      %check if loading the data is possible
      try
        loaddata=str.logical(v.plot_save_data);
      catch
        switch lower(v.plot_save_data)
        case 'force'
          loaddata=false;
          savedata=true;
        otherwise
          error(['Cannot handle parameter ''plot_save_data'' with value ''',v.plot_save_data,'''.'])
        end
      end
      %first build data filename
      datafilename=obj.plotdatafilename(product);
      %check if data is to be loaded
      if loaddata
        %check if plot data is saved
        if ~isempty(datafilename) && file.exist(datafilename)
          str.say('Loading plot data from ',datafilename)
          load(datafilename,'out')
          %we're done
          return
        else
          savedata=true;
        end
      end
      %collect the models 
      out.source.n=0; %this acts as a counter inside the following for loop but becomes a constant thereafter
      out.source.signal=cell(0);out.source.datanames=cell(0);
      for i=1:product.nr_sources
        dn_sources=obj.data_list(product.sources(i));
        %dn_sources is usually 'signal' (and 'error' sometimes)
        for j=1:numel(dn_sources)
          %check if there's any value in the field path to define the type of model
          if ~isempty(dn_sources{j}.field_path)
            %only plot the relevant model types (the metadata 'model_types' defines the relevant model types)
            %NOTICE: this function has not been tested for model_types other than 'signal'
            if ~cells.isincluded(v.model_types,dn_sources{j}.field_path{end}); continue;end
            %get this model type
            model_type=dn_sources{j}.field_path{end};
          else
            %assume signal if no field path to specify the type of model
            model_type='signal';
          end
          %this is legacy but still needed, otherwise the errors associated with the signal will also iterate the counter
          if strcmp(model_type,'signal')
            %only iterate on the signal
            out.source.n=out.source.n+1;
            %save general dataname
            out.source.datanames(out.source.n)=dn_sources(j);
          end
          %save this model_type (most likely out.source.signal{out.source.n})
          out.source.(model_type){out.source.n}=obj.data_get_scalar(dn_sources(j));
          %handle default value of plot_max_degree
          %TODO: fix this, it makes it impossible to plot degree ranges above inclusive maximum
          v.plot_max_degree=min([out.source.(model_type){out.source.n}.lmax,v.plot_max_degree]);
        end
      end
      %enforce maximum degree
      out.source.signal=cellfun(@(i) i.set_lmax(v.plot_max_degree),out.source.signal,'UniformOutput',false);
      %enforce smoothing
      if ~isempty(v.plot_smoothing_degree) && ~isempty(v.plot_smoothing_method)
        v.plot_smoothing_degree=cells.c2m(v.plot_smoothing_degree);
        if isscalar(v.plot_smoothing_degree)
          v.plot_smoothing_degree=v.plot_smoothing_degree*ones(size(out.source.signal));
        else
          assert(numel(v.plot_smoothing_degree)==out.source.n,...
            ['If plot_smoothing_degree is a vector (now with length ',num2str(numel(v.plot_smoothing_degree)),...
            '), it needs to have the same number of elemets as sources (now equal to ',num2str(out.source.n),').'])
        end
        for i=1:out.source.n
          %smooth the data
          out.source.signal{i}=out.source.signal{i}.scale(...
            v.plot_smoothing_degree(i),...
            v.plot_smoothing_method...
          );
          %save smoothing annotations
          smoothing_name=gravity.gauss_smoothing_name(v.plot_smoothing_degree(i));
          out.source.file_smooth{i}=strjoin({'smooth',v.plot_smoothing_method,smoothing_name},'_');
          out.source.title_smooth{i}=str.show({smoothing_name,...
                                  gravity.smoothing_name(v.plot_smoothing_method),'smoothing'});
        end
        %NOTICE: these are used in plots where all sources are shown together
        out.file_smooth=out.source.file_smooth{end};
        out.title_smooth=out.source.title_smooth{end};
      else
        for i=1:out.source.n
          out.source.file_smooth{i}='';
          out.source.title_smooth{i}='';
        end
        %NOTICE: these are used in plots where all sources are shown together
        out.file_smooth='';
        out.title_smooth='';
      end
      %find reference product
      switch v.stats_relative_to
        case 'none'
          out.source.names=cellfun(...
            @(i) strtrim(strrep(str.clean(i,v.plot_title_suppress),'.',' ')),...
            cellfun(@(i) i.name,out.source.datanames,'UniformOutput',false),...
            'UniformOutput',false...
          );
          out.mod.names=out.source.names;
          %alias sources to mods and res
          out.mod.dat=out.source.signal;
          out.mod.res=out.source.signal;
          %patch remaining details
          out.source.ref_idx=[];
          out.source.mod_idx=1:out.source.n;
          out.title_wrt='';
          %get time domain common to all sources
          out.t=simpletimeseries.t_mergev(out.source.signal);
      otherwise
        for i=1:out.source.n
          str.say(out.source.datanames{i}.str,'?=',[v.stats_relative_to,'/'])
          if strfind(out.source.datanames{i}.str,[v.stats_relative_to,'/'])
            out.mod.ref=out.source.signal{i};
            out.source.ref_idx=i;
            out.mod.ref_name=out.source.datanames{out.source.ref_idx};
            break
          end
        end
        assert(isfield(out.source,'ref_idx'),['None of the sources of product ',product.str,...
          ' match ''stats_relative_to'' equal to ''',v.stats_relative_to,'''.'])
        %get the indexes of the sources to plot (i.e. not the reference)
        out.source.mod_idx=[1:out.source.ref_idx-1,out.source.ref_idx+1:out.source.n];
        %get product name difference between products to derive statistics from
        if numel(out.source.mod_idx)<2
          out.mod.names={strjoin(datanames.unique(out.source.datanames(out.source.mod_idx)),' ')};
          out.source.names(out.source.mod_idx)=strrep(out.mod.names,'_',' ');
        else
          out.mod.names=cellfun(@(i) strjoin(i,' '),datanames.unique(out.source.datanames(out.source.mod_idx)),'UniformOutput',false);
          out.source.names(out.source.mod_idx)=cellfun(@(i) strrep(i,'_',' '),out.mod.names,'UniformOutput',false);
        end
        out.source.names{out.source.ref_idx}=upper(strtrim(strrep(str.clean(out.mod.ref_name.name,v.plot_title_suppress),'.',' ')));
        %reduce the models to plot
        out.mod.dat=out.source.signal(out.source.mod_idx);
        %get time domain common to all sources
        out.t=simpletimeseries.t_mergev(out.source.signal);
        %compute residual
        out.mod.res=cellfun(@(i) out.mod.ref.interp(out.t)-i.interp(out.t),out.mod.dat,'UniformOutput',false);
        %title
        out.title_wrt=str.show({'wrt',out.source.names{out.source.ref_idx}});
      end
      %easier names
      out.source.names_str=strjoin(str.rep(out.source.names,' ','_'),'-');
      %enforce spatial mask
      switch v.plot_spatial_mask
        case 'none'
          %do nothing
          out.title_masking='';
        otherwise
          switch v.plot_spatial_mask
          case {'ocean','land'}
            out.title_masking=str.show(v.plot_spatial_mask,'areas');
          otherwise
            out.title_masking=v.plot_spatial_mask;
          end
        %enforce masking
        for i=1:numel(out.mod.dat)
          str.say('Applying',v.plot_spatial_mask,'mask to product',out.mod.names{i})
          out.mod.dat{i}=out.mod.dat{i}.spatial_mask(v.plot_spatial_mask);
          out.mod.res{i}=out.mod.res{i}.spatial_mask(v.plot_spatial_mask);
        end
      end
      %filename particles
      %NOTICE: this is needed so that only one name is returned, otherwise 
      %         a cell array with numerous filenames may be returned, 
      %         depending on the storage period
      out.file_root=product.mdset('storage_period','direct').file('plot',...
        'start',obj.start,...
        'stop',obj.stop,...
        'start_timestamp_only',false,...
        'no_extension',true,...
        'sub_dirs','single',...
        'use_storage_period',true... 
      );
      out.file_deg=['deg',num2str(v.plot_min_degree),'-',num2str(v.plot_max_degree)];
      %legend
      if numel(out.source.datanames)>1
        out.source.legend_str=cellfun(@(i) strjoin(i,' '),datanames.unique(out.source.datanames),'UniformOutput',false);
      else
        out.source.legend_str=str.rep(out.source.datanames{1}.str,'.',' ','/',' ');
      end
      %patch empty legend entries (this expects there to be only one empty legend entry) 
      if any(cells.isempty(out.source.legend_str))
        out.source.legend_str(cells.isempty(out.source.legend_str))={strjoin(datanames.common(out.source.datanames),' ')};
      end
      %time info in the title
      for i=1:out.source.n
        out.source.title_startstop{i}=['(',...
          datestr(out.source.signal{i}.t_masked([],'start'),'yyyy-mm'),' to ',...
          datestr(out.source.signal{i}.t_masked([],'stop' ),'yyyy-mm'),')'...
        ];
      end
      %NOTICE: this is used in plots where all sources are shown together
      out.title_startstop=['(',...
        datestr(out.source.signal{1}.t_masked([],'start'),'yyyy-mm'),' to ',...
        datestr(out.source.signal{1}.t_masked([],'stop' ),'yyyy-mm'),')'...
      ];
      %save data if requested
      if savedata
        str.say('Saving plot data to ',datafilename)
        save(datafilename,'out');  
      end
      %start/stop
      [~,out.startlist,out.stoplist]=product.file('data',v.varargin{:},'start',obj.start,'stop',obj.stop);
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_low_degrees(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_pause_on_save',false,@islogical;...
          'plot_max_lines',inf,@(i) isnumeric(i) && isscalar(i);...
          'show_legend_stats',  'yes' ,@(i) islogical(i) || ischar(i);...
        },...
        product.args...
      },varargin{:});
      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      %build title suffix and legend prefix
      legend_str_prefix=cell(1,out.source.n);
      for k=1:out.source.n
        legend_str_prefix{k}=upper(out.source.legend_str{k});
      end
      if v.plot_legend_include_smoothing
        for k=1:out.source.n
          legend_str_prefix{k}=strjoin([legend_str_prefix(k),out.source.title_smooth(out.source.mod_idx(k))],' ');
        end
        title_suffix=strjoin({out.title_masking},' ');
      else
        title_suffix=strjoin({out.title_masking,out.title_smooth},' ');
      end

      %get collection of degrees/orders (this looks at arguments 'degrees' and 'orders')
      [degrees,orders]=gravity.resolve_degrees_orders(v.varargin{:});
      %loop over all requested degrees and orders
      for i=1:numel(degrees)
        d=degrees(i);
        o=orders(i);
        filename=file.build(out.file_root,['C',num2str(d),',',num2str(o)],'png');
        %plot only if not done yet
        if ~file.exist(filename)
          %build legend string
          legend_str=cell(1,out.source.n);
          trivial_idx=true(size(legend_str));
          %make room for loop records
          ts_now=cell(1,out.source.n); stats=cell(size(ts_now));
          plotting.figure(v.varargin{:});
          %loop over all models
          for k=1:out.source.n
            %plot the time series for this degree/order and model
            ts_now{k}=out.source.signal{k}.ts_C(d,o);
            %don't plot trivial data
            if isempty(ts_now{k}) || ts_now{k}.iszero
              trivial_idx(k)=false;
              stats{k}.corrcoef=-1;stats{k}.rms=inf;
              continue
            end
            %add statistics to the legend (unless this is the produce from which the stats are derived)
            if k~=out.source.ref_idx
              mod_ref_now=out.mod.ref.ts_C(d,o).interp(ts_now{k}.t);
              if all(isnan(mod_ref_now.y(:)))
                stats{k}.corrcoef=NaN;stats{k}.rms=NaN;
              else
                stats{k}=ts_now{k}.stats2(mod_ref_now,'mode','struct','struct_fields',{'corrcoef','rms'},'period',seconds(inf));
              end
              if v.show_legend_stats
                legend_str{k}=[legend_str_prefix{k},...
                  ' corr=',num2str(stats{k}.corrcoef,'%.2f'),...
                  ', RMS{\Delta}=',num2str(stats{k}.rms,'%.2g')];
              else
                legend_str{k}=legend_str_prefix{k};
              end
            else
              if ~isempty(out.source.ref_idx)
                legend_str{k}=[legend_str_prefix{k},' (reference)'];
              else
                legend_str{k}=legend_str_prefix{k};
              end
              stats{k}.corrcoef=1;stats{k}.rms=0;
            end
          end
          %get rid of trivial data
          ts_now=ts_now(trivial_idx);
          legend_str=legend_str(trivial_idx);
          %sort it
          [~,idx]=sort(cellfun(@(i) i.rms,stats),'ascend');
          %truncate it
          idx=idx(1:min(out.source.n,v.plot_max_lines));
          ts_now=ts_now(idx);
          legend_str=legend_str(idx);
          %plot it
          for k=1:numel(ts_now)
            %tweak markers
            switch ts_now{k}.nr_valid
              case 1;    line_fmt='o';
              case 2;    line_fmt='+-';
              otherwise; line_fmt='-';
            end
            ts_now{k}.plot('line',{line_fmt});
          end
          %annotate it
          product.enforce_plot(v.varargin{:},...
            'plot_ylabel','[ ]',...
            'plot_legend',legend_str,...
            'plot_line_color_order',idx,...
            'plot_title',['C',num2str(d),',',num2str(o),' ',title_suffix]...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        else
          disp(['NOTICE: plot already available: ',filename])
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_spatial_stats(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_min_degree',           2, @(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',         inf, @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',     'geoid', @gravity.isfunctional;...
          'plot_type',            'line', @(i) cells.isincluded(i,{'line','bar'});...
          'plot_monthly_error',    false, @islogical;...
          'plot_pause_on_save',    false, @islogical;...
          'plot_show_legend_stats',false, @islogical;...
          'plot_max_nr_lines',        20, @(i) isnumeric(i) && isscalar(i);...
          'plot_derived_quantity',          'cumdrms', @(i) ismethod(gravity.unit(1),i);...
          'plot_spatial_stat_list',{'diff','monthly'}, @iscellstr; ... TODO: corr
          'plot_legend_include_smoothing', false, @islogical;...
        },...
        product.plot_args...
      },varargin{:});
      %make sure this model type is relevant
      if ~isempty(product.dataname.field_path) && ~ismember(product.mdget('model_types'),product.dataname.field_path{end})
        disp(['WARNING: ignoring ',product.codename,' because model type is none of ',...
          strjoin(product.mdget('model_types'),', '),'.'])
        return
      end
      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      
      %check if this plot is requested
      if cells.isincluded(v.plot_spatial_stat_list,'diff')
        
        %% prepare data for (cumulative) diff rms plots
        
        %build filenames
        filenames={...
          file.build(out.file_root,'cumdrms',        out.file_smooth,out.source.names_str,'png');...
          file.build(out.file_root,'cumdrms_summary',out.file_smooth,out.source.names_str,'png');...
        };
        %prepare plot data only if needed
        if any(~file.exist(filenames))
          %build data array
          y=zeros(numel(out.t),numel(out.mod.res));
          for di=1:numel(out.mod.res)
            tmp=out.mod.res{di}.interp(out.t).scale(v.plot_functional,'functional').(v.plot_derived_quantity);
            y(:,di)=tmp.y;
          end
          %average it
          yc=sum(y,1,'omitnan')./sum(~isnan(y));
          %sort it
          [yc_sorted,idx]=sort(yc,'ascend');
          y_sorted=y(:,idx);
          %build legend
          legend_str=upper(out.mod.names);
          if v.plot_legend_include_smoothing
            for i=1:numel(legend_str)
              legend_str{i}=strjoin([legend_str(i),out.source.title_smooth(out.source.mod_idx(i))],' ');
            end
            title_suffix=strjoin({out.title_masking},' ');
          else
            title_suffix=strjoin({out.title_masking,out.title_smooth},' ');
          end
          %sort the legend
          legend_str=legend_str(idx);
          %truncate it
          if v.plot_max_nr_lines < numel(idx)
              y_sorted=  y_sorted(:,1:v.plot_max_nr_lines);
             yc_sorted= yc_sorted(  1:v.plot_max_nr_lines);
            legend_str=legend_str(  1:v.plot_max_nr_lines);
                   idx=idx(         1:v.plot_max_nr_lines);
          end
          %trim nans for plotting
          plot_idx=any(num.trim_NaN(y_sorted),2);
          y_sorted=y_sorted(plot_idx,:);
          t_sorted=out.t(plot_idx);
        end
        
        %% plot diff rms (only if not done yet)
        
        fn_idx=1;
        if ~file.exist(filenames{fn_idx})
          %plot it
          plotting.figure(v.varargin{:});
          switch v.plot_type
          case 'bar';   bar(t_sorted,y_sorted,'EdgeColor','none');
          case 'line'; plot(t_sorted,y_sorted,'Marker','o');
          end
          %deal with legend stats
          if v.plot_show_legend_stats
            for i=1:numel(legend_str)
              legend_str{i}=[legend_str{i},...
                   ' ',num2str(mean(y_sorted(~isnan(y_sorted(:,i)),i)),2),...
               ' +/- ',num2str( std(y_sorted(~isnan(y_sorted(:,i)),i)),2),...
            ' \Sigma=',num2str( sum(y_sorted(~isnan(y_sorted(:,i)),i)),2)];
            end
          end
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',legend_str,...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_xdate',true,...
            'plot_xlimits',[t_sorted(1),t_sorted(end)+days(1)],...
            'plot_line_color_order',idx,...
            'plot_title',str.show({'Residual',out.title_wrt,title_suffix})...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filenames{fn_idx})
          str.say('Plotted',filenames{fn_idx})
        else
          disp(['NOTICE: plot already available: ',filenames{fn_idx}])
        end

        %% plot cumulative diff rms (only if not done yet)
        
        fn_idx=2;
        if ~file.exist(filenames{fn_idx}) && numel(yc_sorted)>1
          %plot it
          plotting.figure(v.varargin{:});
          grey=[0.5 0.5 0.5];
          barh(flipud(yc_sorted(:)),'EdgeColor',grey,'FaceColor',grey);
          set(gca,'YTick',1:numel(legend_str),'yticklabels',str.clean(legend_str,'title'));
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend_location','none',...
            'plot_ylabel','none',...
            'plot_xlabel',gravity.functional_label(v.plot_functional),...
            'plot_ylimits',[0 numel(legend_str)+1],...
            'plot_title',str.show({'Cum. residual',out.title_wrt,out.title_startstop,title_suffix})...
          );
          set(gca,...
            'YTick',1:numel(legend_str),...
            'yticklabels',plotting.legend_replace_clean(v.varargin{:},'plot_legend',flipud(legend_str(:)))...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filenames{fn_idx})
          str.say('Plotted',filenames{fn_idx})
        else
          disp(['NOTICE: plot already available: ',filenames{fn_idx}])
        end
      end
      
      %% epoch-wise degree amplitude plots 
            
      %check if this plot is requested
      if cells.isincluded(v.plot_spatial_stat_list,'monthly')
        for i=1:numel(out.t)
          %gather models valid now
          dat    =cellfun(@(j) j.interp(out.t(i)),out.source.signal,'UniformOutput',false);
          dat_idx=cellfun(@(j) j.nr_valid>0,dat);
          %gather error if requested
          if v.plot_monthly_error
            dat_error=cellfun(@(j) j.interp(out.t(i)),out.source.error,'UniformOutput',false);
          end
          %check if there's any data to plot
          if all(~dat_idx); continue;end
          %reduce data
          dat=dat(dat_idx);
          dat_error=dat_error(dat_idx);
          legend_str=out.source.names(dat_idx);
          %build filename
          filename=cells.scalar(product.file('plot',...
            'start',out.t(i),'stop',out.t(i),...
            'suffix',file.build('drms',out.file_smooth,out.source.names_str)...
          ),'get');
          if isempty(dir(filename))
            plotting.figure(v.varargin{:});
            for j=1:numel(dat)
              dat{j}.plot('mode','drms','functional',v.plot_functional);
            end
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_legend',legend_str,...
              'plot_ylabel',gravity.functional_label(v.plot_functional),...
              'plot_colormap','jet',...
              'plot_title',str.show({datestr(out.t(i),'yyyy-mm'),'degree-RMS',title_suffix})...
            );
            if v.plot_monthly_error
              %get previous lines
              lines_before=findobj(gca,'Type','line');
              %plot errors
              for j=1:numel(dat)
                dat_error{j}.plot('mode','drms','functional',v.plot_functional,'line','--');
              end
              %get all lines
              lines_all=findobj(gca,'Type','line');
              %set consistent line clours
              for j=1:numel(dat)
                set(lines_all(j),'Color',get(lines_before(j),'Color'),'LineWidth',v.plot_line_width)
              end
              axis auto
              title(str.show({datestr(out.t(i),'yyyy-mm'),'degree-RMS ',title_suffix}))
            end
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          else
            disp(['NOTICE: plot already available: ',filename])
          end
        end
      end
      %check if this plot is requested
      if cells.isincluded(v.plot_spatial_stat_list,'corr')
        error('not yet implemented)')
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_temporal_stats(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          'plot_min_degree',     2,          @(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',     inf,        @(i) isnumeric(i) && isscalar(i);...
          'plot_rms_caxis',  [-inf inf], @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',     'geoid',    @gravity.isfunctional;...
          'plot_type',           'line',     @(i) cells.isincluded(i,{'line','bar'});...
          'plot_pause_on_save',  false,      @islogical;...
          'plot_max_nr_lines',   20,         @(i) isnumeric(i) && isscalar(i);...
          'plot_temp_stat_list', {            'corrcoeff',           'rms',           'std'},@(i) iscellstr(i); ...
          'plot_temp_stat_title',{'temporal corr. coeff.','temporal RMS\Delta','temporal STD\Delta'},@(i) iscellstr(i); ...
          'plot_legend_include_smoothing', false, @islogical;...
        },...
        product.plot_args...
      },varargin{:});
      %legacy checking
      assert(~v.isparameter('plot_temp_stat_func'),...
        'the ''plot_temp_stat_func'' metadata entry is no longer supported, use ''plot_functional'' instead.')
      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      
      %% triangular plots
      
      %loop over all statistics
      for s=1:numel(v.plot_temp_stat_list)
        %maybe nothing is requested to be plotted
        if strcmp(v.plot_temp_stat_list{s},'none'); continue; end
        %loop over all sources
        %NOTICE: out.mod.* have one less element than out.source.*, so keep that in mind and don't mix them!
        for i=1:numel(out.mod.dat)
          %build filename
          filename=file.build(out.file_root,...
            {v.plot_temp_stat_list{s},'triang'},out.source.file_smooth{i},strsplit(out.mod.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            plotting.figure(v.varargin{:});
            if isfield(out.mod,'ref')
              %compute the requested stat between this model and mod_ref
              d=out.mod.ref.scale(v.plot_functional,'functional').interp(out.mod.dat{i}.t).stats2(...
                out.mod.dat{i}.scale(v.plot_functional,'functional'),...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'period',seconds(inf)...
              );
            else
              d=out.mod.dat{i}.scale(v.plot_functional,'functional').stats(...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'period',seconds(inf)...
              );
            end
            %set y-units
            d.y_units(:)={gravity.functional_units(v.plot_functional)};
            %plot it
            h=d.plot('method','triang');
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_caxis',cells.c2m(v.(['plot_',v.plot_temp_stat_list{s},'_caxis'])),...
              'plot_title',str.show({...
                upper(out.mod.names{i}),v.plot_temp_stat_title{s},out.title_wrt,...
                out.title_masking,out.source.title_startstop{out.source.mod_idx(i)},...
                strjoin(cells.rm_empty({' ',out.source.title_smooth{out.source.mod_idx(i)}}),newline)...
            }));cb.nan(h.axis_handle);
            %need to adapt caxis label for correlation coefficients
            switch v.plot_temp_stat_list{s}
            case 'corrcoeff'
              cb.label('Corr. Coeff. []');
            end
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
        
        %% prepare data for (cumulative) degree-mean/mean corr coeff plots
        
        %plot degree-mean stat
        filenames={...
          file.build(out.file_root,{v.plot_temp_stat_list{s},'dmean'        },out.file_smooth,out.source.names_str,'png');...
          file.build(out.file_root,{v.plot_temp_stat_list{s},'dmean_summary'},out.file_smooth,out.source.names_str,'png');...
        };
        %prepare plot data only if needed
        if any(~file.exist(filenames))
          %init plot counter and data container
          d=zeros(numel(out.mod.dat),out.mod.dat{1}.lmax+1);
          %loop over all sources
          for i=1:numel(out.mod.dat)
            if isfield(out.mod,'ref')
              %compute it
              d(i,:)=out.mod.ref.scale(v.plot_functional,'functional').interp(out.mod.dat{i}.t).stats2(...
                out.mod.dat{i}.scale(v.plot_functional,'functional'),...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'period',seconds(inf)...
              ).dmean;
            else
              d(i,:)=out.mod.dat{i}.scale(v.plot_functional,'functional').stats(...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'period',seconds(inf)...
              ).dmean;
            end
          end
          %enforce minimum degree
          d=d(:,v.plot_min_degree+1:end);
          %filter out nans
          good_idx=~all(isnan(d),2);
          d=d(good_idx,:);
          %build legend
          legend_str=upper(out.mod.names);
          if v.plot_legend_include_smoothing
            for i=1:numel(legend_str)
              legend_str{i}=strjoin([legend_str(i),out.source.title_smooth(out.source.mod_idx(i))],' ');
            end
            title_suffix=strjoin({...
              out.title_wrt,out.title_masking,out.title_startstop...
            },' ');
          else
            title_suffix=strjoin({...
              out.title_wrt,out.title_masking,out.title_startstop,...
              strjoin(cells.rm_empty({' ',out.title_smooth}),newline)...
            },' ');
          end
          %sort the legend         
          legend_str=legend_str(good_idx);
          %accumulate it
          dc=sum(d,2);
          %need to adapt some plotting aspects to correlation coefficients
          switch v.plot_temp_stat_list{s}
          case 'corrcoeff'
            sort_mode='descend';
          otherwise
            sort_mode='ascend';
          end
          %sort it
          [dc_sorted,idx]=sort(dc,sort_mode);
          d_sorted=d(idx,:);
          legend_str=legend_str(idx);
          %truncate it
          if v.plot_max_nr_lines < numel(idx)
              d_sorted=  d_sorted(1:v.plot_max_nr_lines,:);
             dc_sorted= dc_sorted(1:v.plot_max_nr_lines);
            legend_str=legend_str(1:v.plot_max_nr_lines);
                   idx=idx(       1:v.plot_max_nr_lines);
          end
        end
        
        %% plot degree-mean/mean corr coeff (only if not done yet)
        
        fn_idx=1;
        if ~file.exist(filenames{fn_idx})
          %need to adapt y label for correlation coefficients
          switch v.plot_temp_stat_list{s}
          case 'corrcoeff'
            y_label_str='Corr. Coeff. []';
          otherwise
            y_label_str=gravity.functional_label(v.plot_functional);
          end
          %plot it
          plotting.figure(v.varargin{:});
          switch v.plot_type
          case 'bar';   bar(v.plot_min_degree:out.mod.dat{1}.lmax,d_sorted','EdgeColor','none');
          case 'line'; plot(v.plot_min_degree:out.mod.dat{1}.lmax,d_sorted','Marker','o');
          end
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',legend_str,...
            'plot_ylabel',y_label_str,...
            'plot_xlabel','SH degree',...
            'plot_xlimits',[v.plot_min_degree-1,v.plot_max_degree+1],...
            'plot_line_color_order',idx,...
            'plot_title',str.show({'degree-mean',v.plot_temp_stat_title{s},newline,title_suffix})...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filenames{fn_idx})
          str.say('Plotted',filenames{fn_idx})
        else
          disp(['NOTICE: plot already available: ',filenames{fn_idx}])
        end
        
        %% plot cumulative degree-mean/mean corr coeff stat  (only if not done yet)
        
        fn_idx=2;
        if ~file.exist(filenames{fn_idx}) && numel(dc_sorted)>1
          %need to adapt some plotting aspects to correlation coefficients
          switch v.plot_temp_stat_list{s}
          case 'corrcoeff'
            x_label_str='[]';
          otherwise
            x_label_str=['[',gravity.functional_units(v.plot_functional),']'];
          end
          %plot it
          plotting.figure(v.varargin{:});
          grey=[0.5 0.5 0.5];
          barh(flipud(dc_sorted(:)),'EdgeColor',grey,'FaceColor',grey);
          set(gca,'YTick',1:numel(legend_str),'yticklabels',str.clean(legend_str,'title'));
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend_location','none',...
            'plot_ylabel','none',...
            'plot_xlabel',x_label_str,...
            'plot_ylimits',[0 numel(legend_str)+1],...
            'plot_title',str.show({...
              'Degrees',v.plot_min_degree,'-',v.plot_max_degree,'cum. degree-mean',v.plot_temp_stat_title{s},newline,...
              title_suffix...
             })...
          );
          set(gca,...
            'YTick',1:numel(legend_str),...
            'yticklabels',plotting.legend_replace_clean(v.varargin{:},'plot_legend',flipud(legend_str(:)))...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filenames{fn_idx})
          str.say('Plotted',filenames{fn_idx})
        else
          disp(['NOTICE: plot already available: ',filenames{fn_idx}])
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_stats(obj,product,varargin)
      obj=gswarm.plot_spatial_stats( obj,product,varargin{:});
      obj=gswarm.plot_temporal_stats(obj,product,varargin{:});
      obj=gswarm.plot_low_degrees(obj,product,varargin{:});
    end
    function obj=plot_parametric_decomposition(obj,product,varargin)
      
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          't_mod_f',                                        1, @(i) isnumeric(i) && isscalar(i);... %TODO: implement this 
          'polyorder',                                      2, @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',                           'eqwh', @gravity.isfunctional;...
          'polynames',      {'constant','linear','quadratic'}, @iscellstr;...
          'sin_period',              [      12,            6], @isnumeric;...
          'sin_names',               {'yearly','semi-yearly'}, @iscellstr;...
          'sin_period_unit',                        'months' , @ischar;...
          'timescale',                               'years' , @ischar;...
          'plot_spatial_step',                              1, @(i) isnumeric(i) && isscalar(i);...
          'plot_pause_on_save',                         false, @islogical;...
          'plot_save_data',                             'yes', @(i) islogical(i) || ischar(i);...
        },...
        product.args...
      },varargin{:});
      %some options are derived from others
      v=varargs.wrap('sources',{...
        {...
          'plot_amplitude_range'     inf(size(v.sin_period)), @isnumeric;...
          'plot_poly_range',         inf(1,v.polyorder+1),    @isnumeric;...
        },...
        v...
      });
      %collect the models; NOTICE: all operations in plot_ops are done now, so (e.g.) smoothing can be done in one go
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      
      %check if loading the data is possible
      try
        loaddata=str.logical(v.plot_save_data,'logical');
      catch
        switch lower(v.plot_save_data)
        case 'force';
          loaddata=false;
        otherwise
          error(['Cannot handle parameter ''plot_save_data'' with value ''',v.plot_save_data,'''.'])
        end
      end
      %init data container  
      out.pd=cell(size(out.source.signal));
      
      %loop over all models
      for i=1:out.source.n

        %TODO: check if all datafilenames here are generated consistently
        
        %build data file name
        datafilename=cells.scalar(product.mdset('storage_period','direct').file('plot',...
          'start',obj.start,'stop',obj.stop,...
          'ext','mat',...
          'sub_dirs','single',...
          'start_timestamp_only',false,...
          'suffix',strjoin(cells.rm_empty({out.source.datanames{i}.name,out.source.file_smooth{i},'pardecomp-data'}),'.')...
        ),'get');
        %check if plot data is saved
        if loaddata && ~isempty(datafilename) && file.exist(datafilename)
          str.say('Loading plot data from ',datafilename)
          load(datafilename,'pd')
          out.pd{i}=pd;
        else
          %compute parametric decompositions
          out.pd{i}=pardecomp.split(out.source.signal{i},...
            'np',v.polyorder+1,...
            'T',time.duration2num(time.num2duration(v.sin_period,v.sin_period_unit),v.timescale),...
            'epoch',out.source.signal{i}.epoch,...
            'timescale',v.timescale,...
            't0',time.duration2num(obj.start-out.source.signal{i}.epoch,v.timescale)...
          );
          %split useful bit from the data
          pd=out.pd{i};
          save(datafilename,'pd'); 
          str.say('Saved plot data to ',datafilename)
        end
      
        for j=0:v.polyorder
          par_name=['p',num2str(j)];
          %build filename
          filename=file.build(out.file_root,...
            v.plot_functional,par_name,out.source.file_smooth{i},strsplit(out.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            a=out.pd{i}.(par_name).scale(...
              v.plot_functional,'functional'...
            ).grid('spatial_step',v.plot_spatial_step);
            %build cb titles
            if j==0
              cb_title=gravity.functional_label(v.plot_functional);
            else
              cb_title=[gravity.functional_names(v.plot_functional),' [',gravity.functional_units(v.plot_functional),...
                '/',v.timescale,'^',num2str(j),']'];
            end
            %plot it
            plotting.figure(v.varargin{:});
            a.imagesc(...
              'cb_title',cb_title...
            );
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_legend_location','none',...
              'plot_caxis',[-v.plot_poly_range(j+1),v.plot_poly_range(j+1)],...
              'plot_title',str.show({...
                v.polynames{j+1},'term for',out.source.names{i},...
                out.title_masking,out.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',out.source.title_smooth{i}}),newline)...
            }));
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
        for j=1:numel(v.sin_period)
          %build filename for amplitude
          filename=file.build(out.file_root,...
            v.plot_functional,{'a',j},out.source.file_smooth{i},strsplit(out.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            %get the sine and cosine terms in the form of grids
            s=out.pd{i}.(['s',num2str(j)]).scale(v.plot_functional,'functional').grid('spatial_step',v.plot_spatial_step);
            c=out.pd{i}.(['c',num2str(j)]).scale(v.plot_functional,'functional').grid('spatial_step',v.plot_spatial_step);
            %compute amplitudes
            a=c.power(2)+s.power(2);
            a=a.sqrt;
            %plot it
            plotting.figure(v.varargin{:});
            a.imagesc(...
              'cb_title',gravity.functional_label(v.plot_functional)...
            );
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_legend_location','none',...
              'plot_caxis',[0,v.plot_amplitude_range(j)],...
              'plot_title',str.show({...
                v.sin_names{j},'amplitude for',out.source.names{i},...
                out.title_masking,out.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',out.source.title_smooth{i}}),newline)...
            }));
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
          %build filename for phase
          filename=file.build(out.file_root,...
            v.plot_functional,{'f',j},out.source.file_smooth{i},strsplit(out.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            %get the sine and cosine terms in the form of grids
            s=out.pd{i}.(['s',num2str(j)]).grid('spatial_step',v.plot_spatial_step);
            c=out.pd{i}.(['c',num2str(j)]).grid('spatial_step',v.plot_spatial_step);
            %compute phase
            a=s./c;
            a=a.atan.scale(v.sin_period(j)/pi);
            %plot it
            plotting.figure(v.varargin{:});
            a.imagesc(...
              'cb_title',v.sin_period_unit...
            );
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_legend_location','none',...
              'plot_caxis',[-v.sin_period(j),v.sin_period(j)]/2,...
              'plot_title',str.show({...
                v.sin_names{j},'phase for',out.source.names{i},...
                out.title_masking,out.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',out.source.title_smooth{i}}),newline)...
            }));
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
      end
  
      for i=1:out.source.n
        %build filename        
        filename=file.build(out.file_root,...
          v.plot_functional,'std',out.source.file_smooth{i},strsplit(out.source.names{i},' '),'png'...
        );
        %plot only if not done yet
        if ~file.exist(filename)
          plotting.figure(v.varargin{:});
          out.source.signal{i}.scale(v.plot_functional,'functional').grid('spatial_step',v.plot_spatial_step).stats('mode','std').imagesc(...
            'cb_title',gravity.functional_label(v.plot_functional)...
          );
          plotting.enforce(v.varargin{:},...
            'plot_legend_location','none',...
            'plot_caxis',[0,v.plot_std_range],...
            'plot_colormap','parula',... %this colormap goes well with the red boxes
            'plot_title',str.show({...
              'temporal STD of ',out.source.names{i},...
              out.title_masking,out.source.title_startstop{i},...
              strjoin(cells.rm_empty({' ',out.source.title_smooth{i}}),newline)...
          }));
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      
    end
    function obj=plot_catchments(obj,product,varargin)

      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          'polyorder',                                      2, @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',                           'eqwh', @gravity.isfunctional;...
          'polynames',      {'constant','linear','quadratic'}, @iscellstr;...
          'sin_period',              [      12,            6], @isnumeric;...
          'sin_names',               {'yearly','semi-yearly'}, @iscellstr;...
          'sin_period_unit',                        'months' , @ischar;...
          'timescale',                               'years' , @ischar;...
          'plot_spatial_step',                              1, @(i) isnumeric(i) && isscalar(i);...
          'plot_pause_on_save',                         false, @islogical;...
          'plot_save_data',                             'yes', @(i) islogical(i) || ischar(i);...
          'catchment_list',    simplegrid.catchment_list(:,1), @iscellstr; ....
          'pardecomp_catchments',                        true, @islogical;... %do pardecomp on the catchment timeseries
          'plot_pardecomp_catchments',            {'p0','p1'}, @iscellstr;... %plot these components
          'plot_legend_include_smoothing',              false, @islogical;...
        },...
        product.args...
      },varargin{:});
    
      %collect the models; NOTICE: all operations in plot_ops are done now, so (e.g.) smoothing can be done in one go
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      
      %check if loading the data is possible
      try
        loaddata=str.logical(v.plot_save_data,'logical');
      catch
        switch lower(v.plot_save_data)
        case 'force';
          loaddata=false;
        otherwise
          error(['Cannot handle parameter ''plot_save_data'' with value ''',v.plot_save_data,'''.'])
        end
      end
      %init data container  
      out.catch=cell(size(out.source.signal));
%       %try to load all data (otherwise it's way slower)
%       alldatafilename=cells.scalar(product.mdset('storage_period','direct').file('plot',...
%             'start',obj.start,'stop',obj.stop,...
%             'ext','mat',...
%             'sub_dirs','single',...
%             'start_timestamp_only',false,...
%             'suffix',strjoin(cells.rm_empty({out.file_smooth,'all-catchment-data'}),'.')...
%           ),'get');
%       if loaddata && ~isempty(alldatafilename) && file.exist(alldatafilename)
%         str.say('Loading all plot data from ',alldatafilename)
%         load(alldatafilename,'allcatchments')
%         out.catch=allcatchments;
%       else
        %loop over all catchments
        for j=1:numel(v.catchment_list)
          %loop over all models
          for i=1:out.source.n

            %TODO: check if all datafilenames here are generated consistently

            %build data file name
            datafilename=cells.scalar(product.mdset('storage_period','direct').file('plot',...
              'start',obj.start,'stop',obj.stop,...
              'ext','mat',...
              'sub_dirs','single',...
              'start_timestamp_only',false,...
              'suffix',strjoin(cells.rm_empty({out.source.datanames{i}.name,out.source.file_smooth{i},v.catchment_list{j},'catchment-data'}),'.')...
            ),'get');
            %check if plot data is saved
            if loaddata && ~isempty(datafilename) && file.exist(datafilename)
              str.say('Loading plot data from ',datafilename)
              load(datafilename,'catchment')
              out.catch{i,j}=catchment;
            else
              %compute parametric decompositions
              out.catch{i,j}=out.source.signal{i}.scale(...
                v.plot_functional,'functional'...
              ).grid(...
                'spatial_step',v.plot_spatial_step...
              ).catchment_get(...
                v.catchment_list{j},...
                'parametric_decomposition',v.pardecomp_catchments,...
                'np',v.polyorder+1,...
                'T',time.duration2num(time.num2duration(v.sin_period,v.sin_period_unit),v.timescale),...
                'epoch',out.source.signal{i}.epoch,...
                'timescale',v.timescale,...
                't0',time.duration2num(obj.start-out.source.signal{i}.epoch,v.timescale)...
              );
              %split useful bit from the data
              catchment=out.catch{i,j};
              save(datafilename,'catchment'); 
              str.say('Saved plot data to ',datafilename)
            end
          end
        end
%         allcatchments=out.catch;
%         save(alldatafilename,'allcatchments','-v7.3'); 
%         str.say('Saved all plot data to ',alldatafilename)        
%       end
    
      %show stats header (if relevant)
      msg=cell(out.source.n*(numel(v.catchment_list)+1),1);mc=0;
      tab_len=20;
      mc=mc+1;msg{mc}=str.tablify(tab_len,'catch','solution','bias [cm]','bias diff [cm]','trend [cm/y]','trend diff [cm/y]','corr coeff');
      norm_stats=zeros(out.source.n,5,numel(v.catchment_list));
      for j=1:numel(v.catchment_list)
        %build filename
        filename=file.build(out.file_root,...
          v.plot_functional,'catch',out.file_smooth,strsplit(v.catchment_list{j},' '),'png'...
        );
        %plot only if not done yet
        if ~file.exist(filename)
          ref_idx=out.source.ref_idx; if isempty(ref_idx),ref_idx=1;end
          plotting.figure(v.varargin{:});
          legend_str=cell(1,out.source.n);
          for i=1:out.source.n
            out.catch{i,j}=simplegrid.catchment_plot(out.catch{i,j},...
              'plot_parametric_components',v.plot_pardecomp_catchments...
            );
            legend_str{i}=upper(out.source.names{i});
            if v.plot_legend_include_smoothing
              legend_str{i}=strjoin([legend_str(i),out.source.title_smooth(i)]);
            end
            if i==ref_idx
              l=plotting.line_handles(gca);
              set(l(1),'Color','k')
              if v.pardecomp_catchments; set(l(2),'Color','k'); end
            end
          end
          plotting.enforce(...
            v.varargin{:},...
            'plot_line_style','none',... %these 2 nones are needed so that the formating defined in simplegrid.catchment is not over written
            'plot_line_color','none',...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_title',v.catchment_list{j},...
            'plot_legend',legend_str...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
          
          %show some stats, if parametric decomposition was made
          if isfield(out.catch{i,j},'pws')
            latex_data=cell(out.source.n,6);
            latex_name=cell(out.source.n,1);
            latex_scale=[100 100 years(1)/simpletimeseries.timescale(1)*100 years(1)/simpletimeseries.timescale(1)*100 1];
            %loop over all sources
            for i=1:out.source.n
              %get the names of this solutions and enforce string replacement and cleaning
              latex_name{i}=out.source.names{i};
              latex_name{i}=str.rep(latex_name{i},v.plot_legend_replace{:});
              latex_name{i}=str.clean(latex_name{i},[v.plot_legend_suppress(:);{'title'}]);
              latex_data(i,:)={...
                latex_name{i},...                                                          %solution name
                latex_scale(1)* out.catch{i,j}.pws.p0.y,...                                %bias
                latex_scale(2)*(out.catch{i,j}.pws.p0.y-out.catch{ref_idx,j}.pws.p0.y),... %biasdiff
                latex_scale(3)* out.catch{i,j}.pws.p1.y,...                                %trend cm/year
                latex_scale(4)*(out.catch{i,j}.pws.p1.y-out.catch{ref_idx,j}.pws.p1.y),... %trenddiff
                latex_scale(5)* out.catch{i,j}.ws.interp(out.catch{ref_idx,j}.ws.t).stats2(out.catch{ref_idx,j}.ws,'mode','corrcoeff')...
              };
              mc=mc+1;msg{mc}=str.tablify(tab_len,out.catch{i,j}.name,latex_data{i,:});
              %save some stats for later
              norm_stats(i,:,j)=[latex_data{i,2:end}];
            end
            file.strsave(...
              fullfile(fileparts(filename),str.rep([out.catch{i,j}.name,'.tex'],' ','_')...
            ),str.latex_table(latex_data,{'','%5.1f','%5.1f','%5.1f','%5.1f','%5.2f'}));
          end
        end
      end
      %show stats
      if numel(v.catchment_list)>0 && mc>1; disp(msg);end
      %compute stats for all catchments, if any was computed
      if ~all(norm_stats(i)==0)
        latex_data=cell(out.source.n,4);
        latex_idx=[NaN,3,5,6]-1;
        for i=1:out.source.n
          j=1; latex_data{i,j}=latex_name{i};
          for j=2:3
            latex_data{i,j}=rms(norm_stats(i,latex_idx(j),:));
          end
          j=4;latex_data{i,j}=mean(norm_stats(i,latex_idx(j),:));
        end
        file.strsave('all_catchments.tex',str.latex_table(latex_data,{'','%5.2f','%5.2f','%5.2f'}));
      end
     
    end
    %% utils
    function d=quality(varargin)
      p =dataproduct('swarm.sh.gswarm.rl01.err');
      ps=dataproduct(p.sources(1));
      m='deep ocean';
      debug=true;
      [~,f]=fileparts(strrep(ps.metadata.wilcarded_filename,'*_',''));
      str.say(' --- Input combined models are:',f)
      f=fullfile(getenv('HOME'),'data','gswarm','dissemination',[f,'.quality']);
      str.say(' --- Getting datastorage:')
      d=datastorage(...
        'start',datetime('2013-12-01'),...
        'stop', datetime('2018-09-30'),...
        'inclusive',true,...
        'debug',debug...
      ).init(p);
      str.say(' --- Getting gravity:')
      g =d.data_get_scalar( p.dataname.append_field_leaf('signal'));
      gs=d.data_get_scalar(ps.dataname.append_field_leaf('signal'));
      str.say(' --- Time domain is:')
      disp([g.t_masked,gs.t_masked]);t=g.t_masked;
      f=[f,'.',datestr(t(1),'yyyymm'),'-',datestr(t(end),'yyyymm')];
      str.say(' --- Applying',m,'spatial mask:')
      gm=g.set_lmax(20).spatial_mask('deep ocean');
      str.say(' --- Smoothing, computing geoid and cumulative RMS:')
      c=gm.scale(750e3,'gauss').scale('geoid','functional').cumdrms;
      str.say(' --- Plotting the time series:')
      plotting.figure(varargin{:});
      c.plot;
      plotting.enforce(....
        'plot_legend_location','none'...
      );
      saveas(gcf,[f,'.png'])
      str.say(' --- Exporting  the time series:')
      c.export([f,'.dat'],'ascii')
    end
    function production
      %parameters
      recompute=false;
      
      
      %definitions
%       datafilename=file.unresolve('~/data/gswarm/analyses/2018-06-24/d.mat');
%       p=dataproduct('gswarm.swarm.all.res.plots','plot_dir',fileparts(datafilename));

%       datafilename=file.unresolve('~/data/gswarm/analyses/2018-07-04/d.mat');
%       p=dataproduct('gswarm.swarm.all.res-unsmoothed.plots','plot_dir',fileparts(datafilename));

%       datafilename=file.unresolve('~/data/gswarm/analyses/2018-08-16/d.mat');
%       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...

%       datafilename=file.unresolve('~/data/gswarm/analyses/2018-10-04/d.mat');
%       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...   
%         'gswarm.swarm.all.TN-03_2.catchments',...
%         'gswarm.swarm.all.TN-03_2.unsmoothed',...
%         'gswarm.swarm.all.TN-03_2.smoothed'...
%       },'UniformOutput',false);

%       datafilename=file.unresolve('~/data/gswarm/analyses/2018-11-19/d.mat');
%       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...   
%         'gswarm.swarm.all.TN-03_2.ocean',...
%         'gswarm.swarm.all.TN-03_2.smoothed',...
%       },'UniformOutput',false);
%       plot_stop=datetime('2017-07-31');
     
%         'gswarm.swarm.all.TN-03_2.land',...
%         'gswarm.swarm.all.TN-03_2.ocean',...
% 
%       datadir=file.unresolve('~/data/gswarm/analyses/2019-02-17/');
      datadir=file.unresolve('~/data/gswarm/analyses/2019-04-15/');
      p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datadir)),{...  
        'gswarm.swarm.all.TN-03.catchments',...
        'gswarm.swarm.all.TN-03.smoothed',...
        'gswarm.swarm.all.TN-03.ocean',...
        'gswarm.swarm.all.TN-03.land',...
        'gswarm.swarm.all.TN-03.unsmoothed',...
        'gswarm.swarm.all.TN-03.pardecomp',...
      },'UniformOutput',false);
      start=datetime('2014-07-01');
      stop =datetime('2017-07-31');
      %save version numbers into latex table
      file.strsave(fullfile(datadir,'versions.tex'),gswarm.sourcelatextable(p{1}));

%         'comp.grace.swarm.catchments',...
%         'swarm.sh.gswarm.rl01.land',...
%         'swarm.sh.gswarm.rl01.lowdeg',...
%         'swarm.sh.gswarm.rl01.pardecomp',...
%       datadir=file.unresolve('~/data/gswarm/analyses/2019-04-05/');
%       p=cellfun(@(i) dataproduct(i,'plot_dir',datadir),{...  
%         'swarm.sh.gswarm.rl01.quality',...
%       },'UniformOutput',false);
%       start=datetime('2002-04-01');
%       stop =datetime('2018-12-31');      
     
      if numel(p)==1
        datafilename=fullfile(datadir,[p{1}.str,'.mat']);
      else
        datafilename=fullfile(datadir,'d.mat');
      end
    


      
      %TODO: There's problem with the handling of model errors!!!
      %TODO: can't have GRACE data and Swarm data all the way until Sep 2018
      
      %load data if already available
      if file.exist(datafilename) && ~recompute
        str.say('Loading analysis data from',datafilename)
        load(datafilename,'d')
      else
        d=datastorage('debug',true,'start',start,'stop',stop,'inclusive',true);
        for i=1:p{1}.nr_sources
          d=d.init(p{1}.sources(i),'force',recompute);
        end
        file.ensuredir(datafilename,true);
        save(datafilename,'d')
      end
      %plot it
      for i=1:numel(p)
        d.init(p{i});
      end
      
    end
    function out=sourcelatextable(p)
      sourcelist=cellfun(@(i) str.rep(i,'gswarm.',''),p.source_list_str,'UniformOutput',false);
      out={'\acl{GFM}','version','\acl{KO}'};c=1;
      %count the number of valid entries
      for i=1:numel(sourcelist)
        s=strsplit(sourcelist{i},'.');
        if ~strcmp(s{1},'swarm'); continue; end
        c=c+1;
        if mod(c+1,3)==0; out{c,1}='\rowcolor{Gray}';  c=c+1; end
        switch numel(s)
        case 3
          out{c,1}=s{2}; %gravity field 
          out{c,2}=str.rep(s{3},'v','');     %version
          out{c,3}='\ac{N/A}';               %orbit
        case 4
          out{c,1}=['\ac{',str.rep(upper(s{2}),'IFG','IfG'),'}']; %gravity field 
          out{c,2}=str.rep(s{4},'v','');     %version
          out{c,3}=['\ac{',str.rep(upper(s{3}),'IFG','IfG'),'}']; %orbit
        otherwise
          error(['Cannot handle source ''',sourcelist{i},'''.'])
        end
      end
      out=str.latex_table(out);
    end
  end
end
