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
      % resolve dataname to save the data to: sometimes, the 'signal' field path is given explicitly and
      % we don't want to save data to product.signal.signal and product.signal.error
      dn=product.dataname;
      if any(dn.isanyfield_path(v.model_types))
        for i=1:numel(dn.field_path)
          if any(strcmp(dn.field_path{i},v.model_types))
            if i==1
              dn=dn.set_field_path({});
            else
              dn=dn.set_field_path(dn.field_path(1:i-1));
            end
          end
        end  
      end
      %propagate relevant data
      for i=1:numel(v.model_types)
        switch lower(v.model_types{i})
        case {'signal','sig','s'}
          obj=obj.data_set(dn.append_field_leaf(v.model_types{i}),s);
        case {'error','err','e'}
          obj=obj.data_set(dn.append_field_leaf(v.model_types{i}),e);
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
            %parse using a model or the original data
            if contains(v.use_GRACE_C20,'-model')
              mode='model';
              v.use_GRACE_C20=strrep(v.use_GRACE_C20,'-model','');
            else
              mode='interp';
            end
            %get C20 timeseries, interpolated to current time domain
            c20=gravity.graceC20('version',v.use_GRACE_C20,'mode',mode,'time',mod.t);
  %         figure
  %         plot(c20.x_masked,c20.y_masked([],1),'x-','MarkerSize',10,'LineWidth',4), hold on
  %         plot(c20.x,spline(c20.x_masked,c20.y_masked([],1),c20.x),'o-','MarkerSize',10,'LineWidth',2)
  %         plot(c20.x,pchip(c20.x_masked,c20.y_masked([],1),c20.x),'+-','MarkerSize',10,'LineWidth',2)
  %         legend('data','spline','pchip')
%             %TODO: this is probably not a good idea because any gap in the C20 timeseries is blindly interpolated
%                    better use the model version of it
%             %interpolate in case there are gaps
%             if c20.nr_gaps>0
%               c20=c20.assign([...
%                 interp1(c20.x_masked,c20.y_masked([],1),c20.x),...
%                 interp1(c20.x_masked,c20.y_masked([],2),c20.x)...
%               ],'t',c20.t,'mask',true(size(c20.x)));
%             end
            mod=mod.setC(2,0,c20.y(:,1),mod.t);
            err=err.setC(2,0,c20.y(:,2),mod.t);
%             figure
%             plot(c20.t,c20.y(:,1),'o-'), hold on
%             plot(mod.t,mod.C(2,0),'x-')
%             legend('TN-11',mod.descriptor)
%             msg=['set to ',v.use_GRACE_C20];
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
    function pod=plot_ops(obj,product,varargin)
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
          'model_types',        {'*'} ,@(i) iscellstr(i);...
          'plot_lines_over_gaps_narrower_than',days(120),@isduration;...
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
      %first build data filename; TODO: this datafilename is messed up, the smoothing only
      %appears as 'gauss'
      datafilename=obj.plotdatafilename(product);
      %check if data is to be loaded
      if loaddata
        %check if plot data is saved
        if ~isempty(datafilename) && file.exist(datafilename)
          str.say('Loading plot data from ',datafilename)
          load(datafilename,'out')
          pod=out;
          %we're done
          return
        else
          savedata=true;
        end
      end
      %collect the models 
      pod.source.n=0; %this acts as a counter inside the following for loop but becomes a constant thereafter
      pod.source.dat=cell(0);pod.source.datanames=cell(0);
      for i=1:product.nr_sources
        dn_sources=obj.data_list(product.sources(i));
        %dn_sources is usually 'signal' (and 'error' sometimes) but can be anything
        %NOTICE: all different dn_sources (more specifically the model_type) are clumped together under pod.source.dat
        for j=1:numel(dn_sources)
          %check if there's any value in the field path to define the type of model
          if ~isempty(dn_sources{j}.field_path)
            %model_types {'*'} gathers everything
            if ~cells.isincluded(v.model_types,'*')
              %only plot the relevant model types (the metadata 'model_types' defines the relevant model types)
              if ~cells.isincluded(dn_sources{j}.field_path(end),v.model_types); continue;end
            end
            %get this model type
            model_type=dn_sources{j}.field_path{end};
          else
            %assume signal if no field path to specify the type of model
            model_type='signal';
          end
          %only iterate on the requested model type
          if cells.isincluded({model_type},v.model_types)
            pod.source.n=pod.source.n+1;
            %save general dataname
            pod.source.datanames(pod.source.n)=dn_sources(j);
          end
          %save this model_type (most likely out.source.dat{out.source.n})
          pod.source.dat{pod.source.n}=obj.data_get_scalar(dn_sources(j));
          %handle default value of plot_max_degree
          %TODO: fix this, it makes it impossible to plot degree ranges above inclusive maximum
          v.plot_max_degree=min([pod.source.dat{pod.source.n}.lmax,v.plot_max_degree]);
        end
      end
      %enforce maximum degree
      pod.source.dat=cellfun(@(i) i.set_lmax(v.plot_max_degree),pod.source.dat,'UniformOutput',false);
      %enforce smoothing
      if ~isempty(v.plot_smoothing_degree) && ~isempty(v.plot_smoothing_method) && ~str.none(v.plot_smoothing_method)
        v.plot_smoothing_degree=cells.c2m(v.plot_smoothing_degree);
        if isscalar(v.plot_smoothing_degree)
          v.plot_smoothing_degree=v.plot_smoothing_degree*ones(size(pod.source.dat));
        else
          assert(numel(v.plot_smoothing_degree)==pod.source.n,...
            ['If plot_smoothing_degree is a vector (now with length ',num2str(numel(v.plot_smoothing_degree)),...
            '), it needs to have the same number of elemets as sources (now equal to ',num2str(pod.source.n),').'])
        end
        for i=1:pod.source.n
          %smooth the data
          pod.source.dat{i}=pod.source.dat{i}.scale(...
            v.plot_smoothing_degree(i),...
            v.plot_smoothing_method...
          );
          %save smoothing annotations
          smoothing_name=gravity.gauss_smoothing_name(v.plot_smoothing_degree(i));
          pod.source.file_smooth{i}=strjoin({'smooth',v.plot_smoothing_method,smoothing_name},'_');
          pod.source.title_smooth{i}=str.show({smoothing_name,...
                                  gravity.smoothing_name(v.plot_smoothing_method),'smoothing'});
        end
        %NOTICE: these are used in plots where all sources are shown together
        pod.file_smooth=pod.source.file_smooth{end};
        pod.title_smooth=pod.source.title_smooth{end};
      else
        for i=1:pod.source.n
          pod.source.file_smooth{i}='';
          pod.source.title_smooth{i}='';
        end
        %NOTICE: these are used in plots where all sources are shown together
        pod.file_smooth='';
        pod.title_smooth='';
      end
      %define the time domain to plot
      idx_t=product.mdget('plot_time_domain_source','default',0);
      if idx_t==0
        %set time domain common to all sources
        pod.t=simpletimeseries.t_mergev(pod.source.dat);
      else
        %set time domain equal to the requested source
        pod.t=pod.source.dat{idx_t}.t;
      end
      %find reference product
      switch v.stats_relative_to
        case 'none'
          pod.source.names=cellfun(...
            @(i) strtrim(strrep(str.clean(i,v.plot_title_suppress),'.',' ')),...
            cellfun(@(i) i.name,pod.source.datanames,'UniformOutput',false),...
            'UniformOutput',false...
          );
          pod.mod.names=pod.source.names;
          %alias sources to mods and res
          pod.mod.dat=pod.source.dat;
          pod.mod.res=pod.source.dat;
          %patch remaining details
          pod.source.ref_idx=[];
          pod.source.mod_idx=1:pod.source.n;
          pod.title_wrt='';
      otherwise
        for i=1:pod.source.n
          str.say(pod.source.datanames{i}.filename,'?=',v.stats_relative_to)
          if str.iseq(pod.source.datanames{i}.filename,v.stats_relative_to)
            pod.mod.ref=pod.source.dat{i};
            pod.source.ref_idx=i;
            pod.mod.ref_name=pod.source.datanames{pod.source.ref_idx};
            break
          end
        end
        assert(isfield(pod.source,'ref_idx'),['None of the sources of product ',product.str,...
          ' match ''stats_relative_to'' equal to ''',v.stats_relative_to,'''.'])
        %get the indexes of the sources to plot (i.e. not the reference)
        pod.source.mod_idx=[1:pod.source.ref_idx-1,pod.source.ref_idx+1:pod.source.n];
        %get product name difference between products to derive statistics from
        if numel(pod.source.mod_idx)<2
          pod.mod.names={strjoin(datanames.unique(pod.source.datanames(pod.source.mod_idx)),' ')};
          pod.source.names(pod.source.mod_idx)=strrep(pod.mod.names,'_',' ');
        else
          pod.mod.names=cellfun(@(i) strjoin(i,' '),datanames.unique(pod.source.datanames(pod.source.mod_idx)),'UniformOutput',false);
          pod.source.names(pod.source.mod_idx)=cellfun(@(i) strrep(i,'_',' '),pod.mod.names,'UniformOutput',false);
        end
        pod.source.names{pod.source.ref_idx}=upper(strtrim(strrep(str.clean(pod.mod.ref_name.name,v.plot_title_suppress),'.',' ')));
        %reduce the models to plot
        pod.mod.dat=pod.source.dat(pod.source.mod_idx);
        %compute residual, need to interpolate twice to force explicit gaps
        ref=pod.mod.ref.interp(pod.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than);
        pod.mod.res=cellfun(...
          @(i) ref-i.interp(pod.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than),...
        pod.mod.dat,'UniformOutput',false);
%         pod.mod.res=cell(size(pod.mod.dat));
%         for i=1:numel(pod.mod.dat)
%           dat=pod.mod.dat{i}.interp(pod.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than);
%           pod.mod.res{i}=ref-dat.interp(pod.t);
%         end
%         %NOTICE: probably it's not a good idea to use 'interp_over_gaps_narrower_than' here because that can 
%         %        produce different time domains and break the subtraction.
%         pod.mod.res=cellfun(@(i) pod.mod.ref.interp(pod.t)-i.interp(pod.t),pod.mod.dat,'UniformOutput',false);
        %title
        pod.title_wrt=str.show({'wrt',pod.source.names{pod.source.ref_idx}});
      end
      %easier names
      pod.source.names_str=strjoin(str.rep(str.clean(pod.source.names,'succ_blanks'),' ','_'),'-');
      %enforce spatial mask
      switch v.plot_spatial_mask
        case 'none'
          %do nothing
          pod.title_masking='';
        otherwise
          switch v.plot_spatial_mask
          case {'ocean','land'}
            pod.title_masking=str.show(v.plot_spatial_mask,'areas');
          otherwise
            pod.title_masking=v.plot_spatial_mask;
          end
        %enforce masking
        for i=1:numel(pod.mod.dat)
          str.say('Applying',v.plot_spatial_mask,'mask to product',pod.mod.names{i})
          pod.mod.dat{i}=pod.mod.dat{i}.spatial_mask(v.plot_spatial_mask);
          pod.mod.res{i}=pod.mod.res{i}.spatial_mask(v.plot_spatial_mask);
        end
      end
      %filename particles
      %NOTICE: this is needed so that only one name is returned, otherwise 
      %         a cell array with numerous filenames may be returned, 
      %         depending on the storage period
      pod.file_root=product.mdset('storage_period','direct').file('plot',...
        'start',obj.start,...
        'stop',obj.stop,...
        'start_timestamp_only',false,...
        'no_extension',true,...
        'sub_dirs','single',...
        'use_storage_period',true... 
      );
      pod.file_deg=['deg',num2str(v.plot_min_degree),'-',num2str(v.plot_max_degree)];
      %legend
      if numel(pod.source.datanames)>1
        pod.source.legend_str=cellfun(@(i) strjoin(i,' '),datanames.unique(pod.source.datanames),'UniformOutput',false);
      else
        pod.source.legend_str=str.rep(pod.source.datanames{1}.str,'.',' ','/',' ');
      end
      %patch empty legend entries (this expects there to be only one empty legend entry) 
      if any(cells.isempty(pod.source.legend_str))
        pod.source.legend_str(cells.isempty(pod.source.legend_str))={strjoin(datanames.common(pod.source.datanames),' ')};
      end
      %time info in the title
      for i=1:pod.source.n
        pod.source.title_startstop{i}=['(',...
          datestr(pod.source.dat{i}.t_masked([],'start'),'yyyy-mm'),' to ',...
          datestr(pod.source.dat{i}.t_masked([],'stop' ),'yyyy-mm'),')'...
        ];
      end
      %NOTICE: this is used in plots where all sources are shown together
      pod.title_startstop=['(',...
        datestr(pod.source.dat{1}.t_masked([],'start'),'yyyy-mm'),' to ',...
        datestr(pod.source.dat{1}.t_masked([],'stop' ),'yyyy-mm'),')'...
      ];
      %save data if requested
      if savedata
        str.say('Saving plot data to ',datafilename)
        out=pod; %#ok<NASGU>
        save(datafilename,'out');  
      end
      %start/stop
      [~,pod.startlist,pod.stoplist]=product.file('data',v.varargin{:},'start',obj.start,'stop',obj.stop);
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
          'plot_lines_over_gaps_narrower_than',days(120),@isduration;...
          'plot_legend_sorting','ascend', @ischar;...
        },...
        product.args...
      },varargin{:});
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
      
      %build title suffix and legend prefix
      legend_str_prefix=cell(1,v.pod.source.n);
      for k=1:v.pod.source.n
        legend_str_prefix{k}=upper(v.pod.source.legend_str{k});
      end
      if v.plot_legend_include_smoothing
        for k=1:v.pod.source.n
          legend_str_prefix{k}=strjoin([legend_str_prefix(k),v.pod.source.title_smooth(k)],' ');
        end
        title_suffix=strjoin({v.pod.title_masking},' ');
      else
        title_suffix=strjoin({v.pod.title_masking,v.pod.title_smooth},' ');
      end

      %get collection of degrees/orders (this looks at arguments 'degrees' and 'orders')
      [degrees,orders]=gravity.resolve_degrees_orders(v.varargin{:});
      %loop over all requested degrees and orders
      for i=1:numel(degrees)
        d=degrees(i);
        o=orders(i);
        filename=file.build(v.pod.file_root,['C',num2str(d),',',num2str(o)],'png');
        %plot only if not done yet
        if ~file.exist(filename)
          %build legend string
          legend_str=cell(1,v.pod.source.n);
          trivial_idx=true(size(legend_str));
          %make room for loop records
          ts_now=cell(1,v.pod.source.n); stats=cell(size(ts_now));
          plotting.figure(v.varargin{:});
          %loop over all models
          for k=1:v.pod.source.n
            %plot the time series for this degree/order and model
            ts_now{k}=v.pod.source.dat{k}.ts_C(d,o).addgaps(v.plot_lines_over_gaps_narrower_than);
            %don't plot trivial data
            if isempty(ts_now{k}) || ts_now{k}.iszero
              trivial_idx(k)=false;
              stats{k}.corrcoef=-1;stats{k}.rms=inf;
              continue
            end
            %add statistics to the legend (unless this is the produce from which the stats are derived)
            if k~=v.pod.source.ref_idx
              mod_ref_now=v.pod.mod.ref.ts_C(d,o).interp(ts_now{k}.t);
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
              if ~isempty(v.pod.source.ref_idx)
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
          %sort it (if requested)
          switch v.plot_legend_sorting
          case {'ascend','descend'};   [~,idx]=sort(cellfun(@(i) i.rms,stats),'ascend');
          case {'none','no'};             idx=1:numel(stats);
          otherwise
            error(['Cannot handle ''plot_legend_sorting'' with value ''',v.plot_legend_sorting,'''.'])
          end
          %truncate it
          idx=idx(1:min(v.pod.source.n,v.plot_max_lines));
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
          'plot_derived_quantity',     {'cumdas','gridmean'}, @(i) all(cellfun(@(j) ismethod(gravity.unit(1),j),i));...
          'plot_spatial_stat_list',    {'diff','monthly'}, @iscellstr; ... TODO: corr
          'plot_legend_include_smoothing',          false, @islogical;...
          'plot_lines_over_gaps_narrower_than', days(120), @isduration;...
          'plot_signal',           false, @islogical;...
          'plot_summary',           true, @islogical;...
          'plot_logy',             false, @islogical;...
          'plot_legend_sorting','ascend', @ischar;...
        },...
        product.plot_args...
      },varargin{:});
      %make sure this model type is relevant
      if ~isempty(product.dataname.field_path) && ~ismember(product.mdget('model_types'),product.dataname.field_path{end})
        disp(['WARNING: ignoring ',product.codename,' because model type is none of ',...
          strjoin(product.mdget('model_types'),', '),'.'])
        return
      end
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
      
      %check if this plot is requested
      if cells.isincluded(v.plot_spatial_stat_list,'diff')
        
        %loop over all requested derived quantities
        for qi=1:numel(v.plot_derived_quantity)
          
          %% prepare data for (cumulative) diff rms plots

          %build filenames
          filenames={...
            file.build(v.pod.file_root, v.plot_derived_quantity{qi},            v.pod.file_smooth,v.pod.source.names_str,'png');...
            file.build(v.pod.file_root,[v.plot_derived_quantity{qi},'_summary'],v.pod.file_smooth,v.pod.source.names_str,'png');...
          };
          %prepare plot data only if needed
          if any(~file.exist(filenames))
            %build data array
            y=zeros(numel(v.pod.t),numel(v.pod.mod.res));
            for di=1:numel(v.pod.mod.res)
              [tmp,~,title_name]=v.pod.mod.res{di}.interp(...
                v.pod.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than...
              ).scale(v.plot_functional,'functional').(v.plot_derived_quantity{qi});
              y(:,di)=tmp(:,end);
            end
            if v.plot_signal && ~str.none(product.mdget('stats_relative_to','default','none'))
              [tmp,~,title_name]=v.pod.source.dat{v.pod.source.ref_idx}.interp(...
                v.pod.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than...
              ).scale(v.plot_functional,'functional').(v.plot_derived_quantity{qi});
              y(:,di+1)=tmp(:,end);
            end
            %average it
            yc=sum(y,1,'omitnan')./sum(~isnan(y));
            %sort it (if requested)
            switch v.plot_legend_sorting
            case {'ascend','descend'};   [yc_sorted,idx]=sort(yc,'ascend');
            case {'none','no'};           yc_sorted=yc;idx=1:numel(yc);
            otherwise
              error(['Cannot handle ''plot_legend_sorting'' with value ''',v.plot_legend_sorting,'''.'])
            end
            y_sorted=y(:,idx);
            %build legend
            legend_str=upper(v.pod.mod.names);
            if v.plot_legend_include_smoothing
              for i=1:numel(legend_str)
                legend_str{i}=strjoin([legend_str(i),v.pod.source.title_smooth(v.pod.source.mod_idx(i))],' ');
              end
              title_suffix=strjoin({v.pod.title_masking},' ');
            else
              title_suffix=strjoin({v.pod.title_masking,v.pod.title_smooth},' ');
            end
            %add legend signal
            if v.plot_signal && ~str.none(product.mdget('stats_relative_to','default','none'))
              n=numel(legend_str);
              legend_str{n+1}=upper(v.pod.source.legend_str{v.pod.source.ref_idx});
              for i=1:n
                legend_str{i}=strjoin([legend_str(i),'diff w.r.t.',legend_str(n+1)],' ');
              end
              if v.plot_legend_include_smoothing
                legend_str{n+1}=strjoin([legend_str(n+1),v.pod.source.title_smooth(v.pod.source.ref_idx)],' ');
              end
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
            t_sorted=v.pod.t(plot_idx);
            %compute abs value in case of log scale
            if v.plot_logy
              y_sorted=abs(y_sorted);
            end
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
              'plot_title',str.show({title_name,v.pod.title_wrt,title_suffix})...
            );
            if v.plot_logy; set(gca, 'YScale', 'log'); end
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filenames{fn_idx})
            str.say('Plotted',filenames{fn_idx})
          else
            disp(['NOTICE: plot already available: ',filenames{fn_idx}])
          end

          %% plot cumulative diff rms (only if not done yet)

          fn_idx=2;
          if ~file.exist(filenames{fn_idx}) && numel(yc_sorted)>1 && v.plot_summary
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
              'plot_title',str.show({'Cum. residual',v.pod.title_wrt,v.pod.title_startstop,title_suffix})...
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
      end
      
      %% epoch-wise degree amplitude plots 
            
      %check if this plot is requested
      if cells.isincluded(v.plot_spatial_stat_list,'monthly')
        for i=1:numel(v.pod.t)
          %gather models valid now
          dat    =cellfun(@(j) j.interp(v.pod.t(i)),v.pod.source.dat,'UniformOutput',false);
          dat_idx=cellfun(@(j) j.nr_valid>0,dat);
          %gather error if requested
          if v.plot_monthly_error
            dat_error=cellfun(@(j) j.interp(v.pod.t(i)),v.pod.source.error,'UniformOutput',false);
          end
          %check if there's any data to plot
          if all(~dat_idx); continue;end
          %reduce data
          dat=dat(dat_idx);
          dat_error=dat_error(dat_idx);
          legend_str=v.pod.source.names(dat_idx);
          %build filename
          filename=cells.scalar(product.file('plot',...
            'start',v.pod.t(i),'stop',v.pod.t(i),...
            'suffix',file.build('drms',v.pod.file_smooth,v.pod.source.names_str)...
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
              'plot_title',str.show({datestr(v.pod.t(i),'yyyy-mm'),'degree-RMS',title_suffix})...
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
              title(str.show({datestr(v.pod.t(i),'yyyy-mm'),'degree-RMS ',title_suffix}))
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
          'plot_rms_caxis',      [-inf inf], @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',     'geoid',    @gravity.isfunctional;...
          'plot_type',           'line',     @(i) cells.isincluded(i,{'line','bar'});...
          'plot_pause_on_save',  false,      @islogical;...
          'plot_max_nr_lines',   20,         @(i) isnumeric(i) && isscalar(i);...
          'plot_temp_stat_list', {...
                                    'corrcoeff',...
                                    'rms',...
                                    'std'...
                                  },         @(i) iscellstr(i); ...
          'plot_temp_stat_title',{...
                                    'temporal corr. coeff.',...
                                    'temporal RMS\Delta',...
                                    'temporal STD\Delta'....
                                 },          @(i) iscellstr(i); ...
          'plot_legend_include_smoothing', false, @islogical;...
          'plot_logy',                     false, @islogical;...
          'plot_legend_sorting',        'ascend', @ischar;...
          'plot_corrcoeff_caxis',         [-1,1], @(i) isnumeric(i) && numel(i)==2 && all(abs(i)<=1)
          'plot_summary',                  false, @islogical;...
        },...
        product.plot_args...
      },varargin{:});
      %legacy checking
      assert(~v.isparameter('plot_temp_stat_func'),...
        'the ''plot_temp_stat_func'' metadata entry is no longer supported, use ''plot_functional'' instead.')
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
      
      %% triangular plots
      
      %loop over all statistics
      for s=1:numel(v.plot_temp_stat_list)
        %maybe nothing is requested to be plotted
        if strcmp(v.plot_temp_stat_list{s},'none'); continue; end
        %loop over all sources
        %NOTICE: out.mod.* have one less element than out.source.*, so keep that in mind and don't mix them!
        for i=1:numel(v.pod.mod.dat)
          %build filename
          filename=file.build(v.pod.file_root,...
            {v.plot_temp_stat_list{s},'triang'},v.pod.source.file_smooth{i},strsplit(v.pod.mod.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            plotting.figure(v.varargin{:});
            if isfield(v.pod.mod,'ref')
              %compute the requested stat between this model and mod_ref
              d=v.pod.mod.ref.scale(v.plot_functional,'functional').interp(v.pod.mod.dat{i}.t).stats2(...
                v.pod.mod.dat{i}.scale(v.plot_functional,'functional'),...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'period',seconds(inf)...
              );
            else
              d=v.pod.mod.dat{i}.scale(v.plot_functional,'functional').stats(...
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
                upper(v.pod.mod.names{i}),v.plot_temp_stat_title{s},v.pod.title_wrt,...
                v.pod.title_masking,v.pod.source.title_startstop{v.pod.source.mod_idx(i)},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{v.pod.source.mod_idx(i)}}),newline)...
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
          file.build(v.pod.file_root,{v.plot_temp_stat_list{s},'dmean'        },v.pod.file_smooth,v.pod.source.names_str,'png');...
          file.build(v.pod.file_root,{v.plot_temp_stat_list{s},'dmean_summary'},v.pod.file_smooth,v.pod.source.names_str,'png');...
        };
        %prepare plot data only if needed
        if any(~file.exist(filenames))
          %init plot counter and data container
          d=zeros(numel(v.pod.mod.dat),v.pod.mod.dat{1}.lmax+1);
          %loop over all sources
          for i=1:numel(v.pod.mod.dat)
            if isfield(v.pod.mod,'ref')
              %compute it
              d(i,:)=v.pod.mod.ref.scale(v.plot_functional,'functional').interp(v.pod.mod.dat{i}.t).stats2(...
                v.pod.mod.dat{i}.scale(v.plot_functional,'functional'),...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'period',seconds(inf)...
              ).dmean;
            else
              d(i,:)=v.pod.mod.dat{i}.scale(v.plot_functional,'functional').stats(...
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
          legend_str=upper(v.pod.mod.names);
          if v.plot_legend_include_smoothing
            for i=1:numel(legend_str)
              legend_str{i}=strjoin([legend_str(i),v.pod.source.title_smooth(v.pod.source.mod_idx(i))],' ');
            end
            title_suffix=strjoin({...
              v.pod.title_wrt,v.pod.title_masking,v.pod.title_startstop...
            },' ');
          else
            title_suffix=strjoin({...
              v.pod.title_wrt,v.pod.title_masking,v.pod.title_startstop,...
              strjoin(cells.rm_empty({' ',v.pod.title_smooth}),newline)...
            },' ');
          end
          %sort the legend         
          legend_str=legend_str(good_idx);
          %accumulate it
          dc=sum(d,2);
          %need to adapt some plotting aspects to correlation coefficients
          switch v.plot_temp_stat_list{s}
          case 'corrcoeff'
            switch v.plot_legend_sorting
              case 'ascend';  sort_mode='descend';
              case 'descend'; sort_mode='ascend';
              otherwise;      sort_mode='none';
            end
            plot_logy=false;
          otherwise
            sort_mode=v.plot_legend_sorting;
            plot_logy=v.plot_logy;
          end
          %sort it (if requested)
          switch sort_mode
          case {'ascend','descend'};   [dc_sorted,idx]=sort(dc,sort_mode);
          case {'none','no'};           dc_sorted=dc;idx=1:numel(dc);
          otherwise
            error(['Cannot handle ''plot_legend_sorting'' with value ''',v.plot_legend_sorting,'''.'])
          end
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
          case 'bar';   bar(v.plot_min_degree:v.pod.mod.dat{1}.lmax,d_sorted','EdgeColor','none');
          case 'line'; plot(v.plot_min_degree:v.pod.mod.dat{1}.lmax,d_sorted','Marker','o');
          end
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',legend_str,...
            'plot_ylabel',y_label_str,...
            'plot_xlabel','SH degree',...
            'plot_xlimits',[v.plot_min_degree-1,v.plot_max_degree+1],...
            'plot_line_color_order',idx,...
            'plot_title',str.show({'degree mean',v.plot_temp_stat_title{s},newline,title_suffix})...
          );
          if plot_logy; set(gca, 'YScale', 'log'); end
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filenames{fn_idx})
          str.say('Plotted',filenames{fn_idx})
        else
          disp(['NOTICE: plot already available: ',filenames{fn_idx}])
        end
        
        %% plot cumulative degree-mean/mean corr coeff stat  (only if not done yet)
        
        fn_idx=2;
        if ~file.exist(filenames{fn_idx}) && numel(dc_sorted)>1 && v.plot_summary
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
              'Degrees',v.plot_min_degree,'-',v.plot_max_degree,'cum. degree mean',v.plot_temp_stat_title{s},newline,...
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
    function obj=plot_temporal_std(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          'plot_functional',                            'eqwh', @gravity.isfunctional;...
          'plot_spatial_step',                               1, @(i) isnumeric(i) && isscalar(i);...
          'plot_pause_on_save',                          false, @islogical;...
          'catchment_list',     simplegrid.catchment_list(:,1), @iscellstr; ....
          'plot_std_range'                                 inf, @(i) isnumeric(i) && isscalar(i);...
          'plot_std_colormap'                     'parulazero', @ischar;...
        },...
        product.args...
      },varargin{:});
    
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
    
      for i=1:v.pod.source.n
        %build filename        
        filename=file.build(v.pod.file_root,...
          v.plot_functional,'std',v.pod.source.file_smooth{i},strsplit(v.pod.source.names{i},' '),'png'...
        );
        %plot only if not done yet
        if ~file.exist(filename)
          plotting.figure(v.varargin{:});
          v.pod.source.dat{i}.scale(v.plot_functional,'functional').grid('spatial_step',v.plot_spatial_step).stats('mode','std').imagesc(...
            'boxes',v.catchment_list,...
            'cb_title',gravity.functional_label(v.plot_functional)...
          );
          plotting.enforce(v.varargin{:},...
            'plot_legend_location','none',...
            'plot_caxis',[0,str.num(v.plot_std_range)],...
            'plot_colormap',v.plot_std_colormap,... %the "parula" colormap goes well with the red boxes
            'plot_xlabel','none','plot_ylabel','none',...
            'plot_title',str.show({...
              'temporal STD of ',upper(v.pod.source.names{i}),...
              v.pod.title_masking,v.pod.source.title_startstop{i},...
              strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
          }));
          plotting.no_ticks
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
    end
    function obj=plot_stats(obj,product,varargin)
      %collect the models
      pod=gswarm.plot_ops(obj,product,varargin{:});
      %call the lower tier plotting routines, re-use pod
      obj=gswarm.plot_spatial_stats( obj,product,varargin{:},'pod',pod);
      obj=gswarm.plot_temporal_stats(obj,product,varargin{:},'pod',pod);
      obj=gswarm.plot_low_degrees(obj,product,varargin{:},   'pod',pod);
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
          't0',                        datetime('2000-01-01'), @isdatetime;...
          'plot_spatial_step',                              1, @(i) isnumeric(i) && isscalar(i);...
          'plot_pause_on_save',                         false, @islogical;...
          'plot_save_data',                             'yes', @(i) islogical(i) || ischar(i);...
          'plot_catchment_list',simplegrid.catchment_list(:,1),@iscellstr; ....
          'plot_catchment_boxes',                       false, @islogical; ....
          'plot_std_colormap',                       'parula', @ischar;...
          'plot_std_range'                                inf, @(i) isnumeric(i) && isscalar(i);...
        },...
        product.args...
      },varargin{:});
      %some options are derived from others
      v=varargs.wrap('sources',{...
        {...
          'plot_amplitude_range'     inf(size(v.sin_period)), @isnumeric;...
          'plot_poly_range',            inf(1,v.polyorder+1), @isnumeric;...
          'pod',                                          [], @isstruct;...
        },...
        v...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
      
      %check if loading the data is possible
      try
        loaddata=str.logical(v.plot_save_data,'logical');
      catch
        switch lower(v.plot_save_data)
        case 'force'
          loaddata=false;
        otherwise
          error(['Cannot handle parameter ''plot_save_data'' with value ''',v.plot_save_data,'''.'])
        end
      end
      %init data container  
      v.pod.pd=cell(size(v.pod.source.dat));
      
      %loop over all models
      for i=1:v.pod.source.n

        %TODO: check if all datafilenames here are generated consistently
        
        %build data file name
        datafilename=cells.scalar(product.mdset('storage_period','direct').file('plot',...
          'start',obj.start,'stop',obj.stop,...
          'ext','mat',...
          'sub_dirs','single',...
          'start_timestamp_only',false,...
          'suffix',strjoin(cells.rm_empty({v.pod.source.datanames{i}.name,v.pod.source.file_smooth{i},'pardecomp-data'}),'.')...
        ),'get');
        %check if plot data is saved
        if loaddata && ~isempty(datafilename) && file.exist(datafilename)
          str.say('Loading plot data from ',datafilename)
          load(datafilename,'pd')
          v.pod.pd{i}=pd;
        else
          %compute parametric decompositions
          v.pod.pd{i}=pardecomp.split(v.pod.source.dat{i},...
            'np',v.polyorder+1,...
            'T',time.duration2num(time.num2duration(v.sin_period,v.sin_period_unit),v.timescale),...
            'epoch',v.pod.source.dat{i}.epoch,...
            'timescale',v.timescale,...
            't0',time.duration2num(obj.start-v.t0,v.timescale)...
          );
          %split useful bit from the data
          pd=v.pod.pd{i};
          save(datafilename,'pd'); 
          str.say('Saved plot data to ',datafilename)
        end
      
        for j=0:v.polyorder
          par_name=['p',num2str(j)];
          %build filename
          filename=file.build(v.pod.file_root,...
            v.plot_functional,par_name,v.pod.source.file_smooth{i},strsplit(v.pod.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            a=v.pod.pd{i}.(par_name).scale(...
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
              'plot_ylabel','none','plot_xlabel','none',...
              'plot_caxis',[-v.plot_poly_range(j+1),v.plot_poly_range(j+1)],...
              'plot_title',str.show({...
                v.polynames{j+1},'term for',v.pod.source.names{i},...
                v.pod.title_masking,v.pod.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
            }));
            plotting.no_ticks
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
        for j=1:numel(v.sin_period)
          %build filename for amplitude
          filename=file.build(v.pod.file_root,...
            v.plot_functional,{'a',j},v.pod.source.file_smooth{i},strsplit(v.pod.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            %get the sine and cosine terms in the form of grids
            s=v.pod.pd{i}.(['s',num2str(j)]).scale(v.plot_functional,'functional').grid('spatial_step',v.plot_spatial_step);
            c=v.pod.pd{i}.(['c',num2str(j)]).scale(v.plot_functional,'functional').grid('spatial_step',v.plot_spatial_step);
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
              'plot_ylabel','none','plot_xlabel','none',...
              'plot_caxis',[0,v.plot_amplitude_range(j)],...
              'plot_title',str.show({...
                v.sin_names{j},'amplitude for',v.pod.source.names{i},...
                v.pod.title_masking,v.pod.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
            }));
            plotting.no_ticks
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
          %build filename for phase
          filename=file.build(v.pod.file_root,...
            v.plot_functional,{'f',j},v.pod.source.file_smooth{i},strsplit(v.pod.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename)
            %get the sine and cosine terms in the form of grids
            s=v.pod.pd{i}.(['s',num2str(j)]).grid('spatial_step',v.plot_spatial_step);
            c=v.pod.pd{i}.(['c',num2str(j)]).grid('spatial_step',v.plot_spatial_step);
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
              'plot_ylabel','none','plot_xlabel','none',...
              'plot_caxis',[-v.sin_period(j),v.sin_period(j)]/2,...
              'plot_title',str.show({...
                v.sin_names{j},'phase for',v.pod.source.names{i},...
                v.pod.title_masking,v.pod.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
            }));
            plotting.no_ticks
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
      end
  
      if ~v.plot_catchment_boxes
        v.plot_catchment_list={};
      end
      
    end
    function obj=plot_catchments(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          'plot_functional',                           'eqwh', @gravity.isfunctional;...
          'plot_spatial_step',                              1, @(i) isnumeric(i) && isscalar(i);...
          'plot_pause_on_save',                         false, @islogical;...
          'plot_save_data',                             'yes', @(i) islogical(i) || ischar(i);...
          'catchment_list',    simplegrid.catchment_list(:,1), @iscellstr; ....
          'parametric_decomposition',                    true, @islogical;... %do pardecomp on the catchment timeseries
          'polyorder',                                      2, @(i) isnumeric(i) && isscalar(i);...
          'polynames',          {'bias','linear','quadratic'}, @iscellstr;...
          'sin_period',              [      12,            6], @isnumeric;...
          'sin_period_unit',                        'months' , @ischar;...
          'sin_names',               {'yearly','semi-yearly'}, @iscellstr;...
          'timescale',                               'years' , @ischar;...
          't0',                        datetime('2000-01-01'), @isdatetime;...
          'plot_parametric_components',           {'p0','p1'}, @iscellstr;... %plot these components
          'plot_parametric_timestep_value',                 7, @(i) isnumeric(i) && isscalar(i);...
          'plot_parametric_timestep_units',            'days', @ischar;...
          'plot_legend_include_smoothing',              false, @islogical;...
          'plot_line_colormap',                      'spiral', @ischar;...
          'plot_lines_over_gaps_narrower_than',     days(120), @isduration;...
        },...
        product.args...
      },varargin{:});
    
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
  
      %check if loading the data is possible
      try
        loaddata=str.logical(v.plot_save_data,'logical');
      catch
        switch lower(v.plot_save_data)
        case 'force'
          loaddata=false;
        otherwise
          error(['Cannot handle parameter ''plot_save_data'' with value ''',v.plot_save_data,'''.'])
        end
      end
      %support legacy
      if isfield(v.pod.source,'dat')
        sfn='dat';
      else
        sfn='signal';
      end
      %init data container  
      v.pod.catch=cell(size(v.pod.source.(sfn)));
      %init user feedback
      s.msg='Loading/computing catchments';s.n=numel(v.catchment_list)*v.pod.source.n;s.stationary=false;
      pd_args={...
        'T',time.duration2num(time.num2duration(v.sin_period,v.sin_period_unit),v.timescale),...
        'np',v.polyorder+1,...
        't0',time.duration2num(obj.start-v.t0,v.timescale)...
        'epoch',v.pod.source.(sfn){1}.epoch,...
        'timescale',v.timescale,...
        'descriptor',v.pod.source.(sfn){1}.descriptor,...
      };
      if v.parametric_decomposition
        disp(['Parametric decomposition will consider the following parameters:',newline,pardecomp.str(pd_args{:})])
      end
      %loop over all catchments
      for j=1:numel(v.catchment_list)
        %check what needs to be loaded/computed
        for i=1:v.pod.source.n
          %build data file name
          datafilename=cells.scalar(product.mdset('storage_period','direct').file('plot',...
            'start',obj.start,'stop',obj.stop,...
            'ext','mat',...
            'sub_dirs','single',...
            'start_timestamp_only',false,...
            'suffix',strjoin(cells.rm_empty({v.pod.source.datanames{i}.name,v.pod.source.file_smooth{i},v.catchment_list{j},'catchment-data'}),'.')...
          ),'get');
          if loaddata && ~isempty(datafilename) && file.exist(datafilename)
            str.say('Loading catchment data from ',datafilename)
            load(datafilename,'catchment')
            v.pod.catch{i,j}=catchment;
          else
            %compute parametric decompositions
            v.pod.catch{i,j}=v.pod.source.(sfn){i}.scale(...
              v.plot_functional,'functional'...
            ).grid(...
              'spatial_step',v.plot_spatial_step...
            ).catchment_get(...
              v.catchment_list{j},...
              pd_args{:},...
              'parametric_decomposition',v.parametric_decomposition,...
              'epoch',v.pod.source.(sfn){i}.epoch...
            );
            %split useful bit from the data
            catchment=v.pod.catch{i,j};
            save(datafilename,'catchment'); 
            str.say('Saved catchment data to ',datafilename)
          end
          %user feedback
          s=time.progress(s);
        end
      end
    
      %parameters
      ref_idx=v.pod.source.ref_idx; if isempty(ref_idx),ref_idx=1;end

      for j=1:numel(v.catchment_list)
        %build filename
        filename=file.build(v.pod.file_root,...
          v.plot_functional,'catch',v.pod.file_smooth,strsplit(v.catchment_list{j},' '),'png'...
        );
        %plot only if not done yet
        if ~file.exist(filename)
          plotting.figure(v.varargin{:});
          legend_str=cell(1,v.pod.source.n);
          %the default color scheme for lines only has a limited number of entries, so to
          %be safe against that, just build a colormap from scratch
          plot_line_colormap=plotting.line_color_map(v.plot_line_colormap,v.pod.source.n);
          %loop over all sources
          for i=1:v.pod.source.n
            v.pod.catch{i,j}=simplegrid.catchment_plot(v.pod.catch{i,j},...
              'time',...
                v.pod.catch{i,j}.pws.time(1):...
                time.num2duration(v.plot_parametric_timestep_value,v.plot_parametric_timestep_units):...
                v.pod.catch{i,j}.pws.time(end),...
              'color',plot_line_colormap(i),...
              'plot_lines_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than,...
              'plot_parametric_components',v.plot_parametric_components...
            );
            legend_str{i}=upper(v.pod.source.names{i});
            if v.plot_legend_include_smoothing
              legend_str{i}=strjoin([legend_str(i),v.pod.source.title_smooth(i)]);
            end
            if i==ref_idx
              l=plotting.line_handles(gca);
              set(l(1),'Color','k')
              if v.parametric_decomposition; set(l(2),'Color','k'); end
            end
          end
          plotting.enforce(...
            v.varargin{:},...
            'plot_line_style','none',... %these 2 nones are needed so that the formating defined in simplegrid.catchment_plot is not over written
            'plot_line_color','none',...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_xlabel','none',...
            'plot_title',v.catchment_list{j},...
            'plot_legend',legend_str...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      
      %show stats header (if relevant)
      msg=cell(v.pod.source.n*(numel(v.catchment_list)+1),1);mc=0;
      tab_len=20;
      mc=mc+1;msg{mc}=str.tablify(tab_len,'catch','solution','bias [cm]','bias diff [cm]','trend [cm/y]','trend diff [cm/y]','corr coeff');
      norm_stats=zeros(v.pod.source.n,5,numel(v.catchment_list));
      for j=1:numel(v.catchment_list)
        %build filename
        filename=file.build(v.pod.file_root,...
          v.plot_functional,'catch',v.pod.file_smooth,strsplit(v.catchment_list{j},' '),'tex'...
        );
        %plot only if not done yet
        if ~file.exist(filename)
          %show some stats, if parametric decomposition was made
          if isfield(v.pod.catch{i,j},'pws')
            latex_data=cell(v.pod.source.n,6);
            latex_name=cell(v.pod.source.n,1);
            latex_scale=[100 100 100 100 1]; %converts m to cm
            %compute bias of the reference (the p0 value is only relevant at t0 and it doesn't really give a good idea of the bias)
            v.pod.catch{ref_idx,j}.mean=mean(v.pod.catch{ref_idx,j}.ws.y_masked);
            v.pod.catch{      i,j}.mean=mean(v.pod.catch{      i,j}.ws.y_masked);
            %loop over all sources
            for i=1:v.pod.source.n
              %get the names of this solutions and enforce string replacement and cleaning
              latex_name{i}=v.pod.source.names{i};
              latex_name{i}=str.rep(latex_name{i},v.plot_legend_replace{:});
              latex_name{i}=str.clean(latex_name{i},[v.plot_legend_suppress(:);{'title'}]);
              %compute bias of this source
              v.pod.catch{i,j}.bias=mean(v.pod.catch{i,j}.ws.minus(v.pod.catch{i,j}.pws.ts_p1).minus(v.pod.catch{i,j}.pws.ts_p0).y_masked);
              latex_data(i,:)={...
                latex_name{i},...                                                              %solution name
                latex_scale(1)* v.pod.catch{i,j}.mean,...                                      %bias
                latex_scale(2)*(v.pod.catch{i,j}.mean-v.pod.catch{ref_idx,j}.mean),...         %biasdiff
                latex_scale(3)* v.pod.catch{i,j}.pws.p1.y,...                                  %trend cm/year
                latex_scale(4)*(v.pod.catch{i,j}.pws.p1.y-v.pod.catch{ref_idx,j}.pws.p1.y),... %trenddiff
                latex_scale(5)* v.pod.catch{i,j}.ws.interp(v.pod.catch{ref_idx,j}.ws.t).stats2(v.pod.catch{ref_idx,j}.ws,'mode','corrcoeff')...
              };
              mc=mc+1;msg{mc}=str.tablify(tab_len,v.pod.catch{i,j}.name,latex_data{i,:});
              %save some stats for later
              norm_stats(i,:,j)=[latex_data{i,2:end}];
            end
            %enforce legend editing
            latex_data(:,1)=plotting.legend_replace_clean(v.varargin{:},'plot_legend',upper(latex_data(:,1)));
            %save table
            file.strsave(filename,str.latex_table(latex_data,{'','%6.2f','%6.2f','%6.2f','%6.2f','%5.2f'}));
          end
        end
      end
      %show stats
      if numel(v.catchment_list)>0 && mc>1; disp(msg);end
      %build filename
      filename=file.build(v.pod.file_root,'all_catchments.tex');
      %compute stats for all catchments, if any was computed
      if ~all(norm_stats(i)==0)
        latex_data=cell(v.pod.source.n,4);
        latex_idx=[NaN,3,5,6]-1;
        for i=1:v.pod.source.n
          j=1; latex_data{i,j}=latex_name{i};
          for j=2:3
            latex_data{i,j}=rms(norm_stats(i,latex_idx(j),:));
          end
          j=4;latex_data{i,j}=mean(norm_stats(i,latex_idx(j),:));
        end
        file.strsave(filename,str.latex_table(latex_data,{'','%5.2f','%5.2f','%5.2f'}));
      end
     
      %plot temporal std to show catchments
      obj=gswarm.plot_temporal_std(obj,product,varargin{:},'pod',v.pod);
      
    end
    %% utils
    function d=grace_model(varargin)
      v=varargs.wrap('sources',{....
        {...
          'debug',     true,                      @islogical;...
          'force',     false,                     @islogical;...
          'start',     datetime('2002-04-01'),    @isdatetime;...
          'stop',      datetime('2017-06-30'),    @isdatetime;...
          'functional','eqwh',                    @ischar;...
          'smoothing', 350e3,                     @numeric;...
          'mode',      'all',                     @(i)ischar(i)||iscellstr(i);...
          'products',  {...
            'grace.sh.rl06.csr';...
            'grace.sh.rl06.csr.pd';...
            'grace.sh.rl06.csr.pd.ts';...
          },                                      @iscellstr;...
        },... 
      },varargin{:});
      %create vector with relevant data product
      p=cellfun(@(i) dataproduct(i),v.products,'UniformOutput',false);
      %one does not usually want to ignore an explicit force flag, so be sure that failsafes are off
      p=cellfun(@(i) i.mdset('never_force',false,'always_force',false),p,'UniformOutput',false);
      %may need to init grace climatological model if force is true (this needs to be done separately from the Swarm processing
      %because otherwise the Swarm period is what is used in the regression: no bueno)
      d=datastorage(...
        'inclusive',true,...
        'start',v.start,...
        'stop', v.stop,...
        'debug',v.debug...
      );
      for i=1:numel(p)
        d=d.init(p{i},'force',v.force);
      end
      %derive data
      a={
        d.data_get_scalar(p{1}.dataname.append_field_leaf('signal')),'original' ;... 
        pardecomp.join(d.data.grace_sh_rl06_csr_pd),              're-modelled' ;...
        d.data_get_scalar(p{3}),                                        'model' ;...
      };
      %plot C20
      if cells.isincluded(v.mode,{'C20','all'})
        stmn=dataproduct(p{1}).mdget('static_model');
        stc20=datastorage().init(stmn).data_get_scalar(datanames(stmn).append_field_leaf('signal')).C(2,0);
        tn=gravity.graceC20('mode','read' ); tn=tn.trim(v.start,v.stop).addgaps(days(45))-stc20(1);            
        plotting.figure
        for i=1:size(a,1)
          a{i,1}.ts_C(2,0).addgaps(days(45)).plot;
        end
        tn.plot('columns',1)
        plotting.enforce('plot_legend',[a(:,2);{'TN-11'}],'plot_title','C20',...
          'plot_ylabel',gravity.functional_label(v.functional))
      end
      if cells.isincluded(v.mode,{'freq-stats','all'})
        for i=1:size(a,1)
          a{i,3}=a{i,1}.scale(v.functional,'functional');
        end
        %plot cumulative degree stats
        for s={'mean','std','rms'}
          plotting.figure
          c=0;legend_str=cell(0);
          i=1;
          [~,ts]=a{i,3}.(['cumd',s{1}]); ts.addgaps(days(45)).plot;
          c=c+1; legend_str{c}='with C20 (TN-11)';
          [~,ts]=a{i,3}.scale(v.smoothing,'gauss').(['cumd',s{1}]); ts.addgaps(days(45)).plot;
          c=c+1; legend_str{c}=str.show({'with C20',v.smoothing/1e3,'km gauss'});
          [~,ts]=a{i,3}.setC(2,0,0).(['cumd',s{1}]); ts.addgaps(days(45)).plot;
          c=c+1; legend_str{c}='without C20';
          i=2;
          [~,ts]=a{i,3}.(['cumd',s{1}]); ts.addgaps(days(45)).plot;
          c=c+1; legend_str{c}='p12 with C20';
          plotting.enforce('plot_legend',legend_str,'plot_title',['cumulative degree ',s{1}],...
            'plot_ylabel',gravity.functional_label(v.functional))
        end
      end
      if cells.isincluded(v.mode,{'spatial-stats','all'})
        for i=1:size(a,1)
          a{i,4}=a{i,1}.scale(v.functional,'functional').detrend.grid;
        end
        %plot spatial stats (TODO: the mean should be zero!)
        for s={'mean','std','rms'}
          plotting.figure;
          for i=1:size(a,1)
            [~,~,ts]=a{i,4}.(s{1})('tot'); ts.addgaps(days(45)).plot;
          end
          plotting.enforce('plot_legend',a(:,2),'plot_title',['spatial ',s{1}],...
            'plot_ylabel',gravity.functional_label(v.functional))
        end
      end
    end
    function d=quality(varargin)
      %NOTICE: this is used to produce the plots in ~/data/gswarm/dissemination
      %TODO: merge this with the validation method (and augment the report there)
      %parse inputs
      v=varargs.wrap('sources',{....
        {...
          'force_grace_pd',         false, @islogical;...
          'force',                  false, @islogical;...
          'mask',            'deep ocean', @ischar;...
          'debug',                   true, @ischar;...
          'start', datetime('2013-11-01'), @isdatestime;...
          'stop',  datetime('2019-06-30'), @isdatestime;...
          'smoothing',              750e3, @isnumeric;...
          'functional',           'geoid', @ischar;...
        },... 
      },varargin{:});
      %define Swarm products
      p =dataproduct('swarm.sh.gswarm.rl01.err');
      ps=dataproduct(p.sources(1));
      %fileame stem
      [~,f]=fileparts(strrep(ps.metadata.wilcarded_filename,'*_',''));
      str.say(' --- Input combined models are:',f)
      f=fullfile(getenv('HOME'),'data','gswarm','dissemination',...
        strjoin({...
          f,datestr(now,'yyyy-mm-dd'),[datestr(v.start,'yyyymm'),'-',datestr(v.stop,'yyyymm')],...
          num2str(v.smoothing/1e3)...
        },'.')...
      );
      %re-process GRACE data if requested
      if v.force_grace_pd
        p_now='grace.sh.rl06.csr.pd.ts';
        str.say(' --- re-computing :',p_now)
        d=datastorage(...
          'start',datetime('2002-04-01'),...
          'stop', datetime('2017-06-30'),...
          'debug',v.debug...
        ).init(p_now,'force',true);
%         ).init(p_now,'force',false).init('grace.sh.rl06.csr');
        %make pretty plots
        plotting.figure(varargin{:});
        grs=d.data_get_scalar('grace_sh_rl06_csr/signal').ts_C(2,0).plot('zeromean',true);
        grm=d.data_get_scalar('grace_sh_rl06_csr_pd_ts' ).ts_C(2,0).plot('zeromean',true);
        tn=gravity.graceC20.plot('columns',1,'zeromean',true);
        plotting.enforce('plot_legend',{...
          str.rep(grs.legend{1},'C2,0','GRACE data');...
          str.rep(grm.legend{1},'C2,0','GRACE model');...
          str.rep(tn.legend{1},'C20','TN-11 (w/ static)');...
        },'plot_title','C20','plot_ylabel','non-dim');
        saveas(gcf,[f,'.grace-C20.png'])
        plotting.figure(varargin{:});
        grs=d.data_get_scalar('grace_sh_rl06_csr/signal').ts_C(3,0).plot('zeromean',true);
        grm=d.data_get_scalar('grace_sh_rl06_csr_pd_ts' ).ts_C(3,0).plot('zeromean',true);
        plotting.enforce('plot_legend',{...
          str.rep(grs.legend{1},'C3,0','GRACE data');...
          str.rep(grm.legend{1},'C3,0','GRACE model');...
        },'plot_title','C30','plot_ylabel','non-dim');
        saveas(gcf,[f,'.grace-C30.png'])
      end
      %reprocess Swarm products if needed
      if v.force
%         str.say(' --- re-computing grace.sh.rl06.csr for the GRACE-FO period')
%         datastorage(...
%           'start',datetime('2018-06-01'),...
%           'stop', v.stop,...
%           'debug',v.debug...
%         ).init('grace.sh.rl06.csr','force',true);
        str.say(' --- re-computing :',p)
        d=datastorage(...
          'inclusive',false,...
          'start',datetime('2013-10-01'),...
          'stop', v.stop,...
          'debug',v.debug...
        ).init(dataproduct('grace.sh.rl06.csr').mdset('never_force',true)...
        ).init(dataproduct('grace.sh.rl06.csr.pd.ts').mdset('never_force',true)...
        ).init(dataproduct('swarm.sh.gswarm.rl01')...
        ).init(p,'force',false);
        %need ggm05c (make sure this is the static model used in all products)
        stmn=dataproduct('swarm.sh.gswarm.rl01').mdget('static_model');
        stc20=datastorage().init(stmn).data_get_scalar(datanames(stmn).append_field_leaf('signal')).C(2,0);
        %make pretty plots
        plotting.figure(varargin{:});
        sw=d.data_get_scalar('swarm_sh_gswarm_rl01/signal' ).ts_C(2,0).addgaps(days(45)); swp=sw.plot;
        gs=d.data_get_scalar('grace_sh_rl06_csr/signal'    ).ts_C(2,0).addgaps(days(45)); gsp=gs.plot;
        gm=d.data_get_scalar('grace_sh_rl06_csr_pd_ts'     ).ts_C(2,0).addgaps(days(45)); gmp=gm.plot;
        tn=gravity.graceC20('mode','read' ); tn=tn.trim(v.start,v.stop).addgaps(days(45))-stc20(1);            
        tnp=tn.plot('columns',1);
        tm=gravity.graceC20('mode','model','time',gm.t);tm=tm.addgaps(days(45))-stc20(1);
        tmp=tm.plot('columns',1);
        plotting.enforce('plot_legend',{...
          str.rep(swp.legend{1},'C2,0','Swarm');...
          str.rep(gsp.legend{1},'C2,0','GRACE data');...
          str.rep(gmp.legend{1},'C2,0','GRACE model (p12)');...
          str.rep(tnp.legend{1},'C20','TN-11');...
          ['TN-11 model (p22)',tmp.legend{1}];...
        },'plot_title','C20','plot_ylabel','non-dim');
        saveas(gcf,[f,'.swarm-C20.png'])
        plotting.figure(varargin{:});
        sw=d.data_get_scalar('swarm_sh_gswarm_rl01/signal').ts_C(3,0).addgaps(days(45)).plot;
        gs=d.data_get_scalar('grace_sh_rl06_csr/signal'   ).ts_C(3,0).addgaps(days(45)).plot;
        gm=d.data_get_scalar('grace_sh_rl06_csr_pd_ts'    ).ts_C(3,0).addgaps(days(45)).plot;
        plotting.enforce('plot_legend',{...
          str.rep(sw.legend{1},'C3,0','Swarm');...
          str.rep(gs.legend{1},'C3,0','GRACE data');...
          str.rep(gm.legend{1},'C3,0','GRACE model');...
        },'plot_title','C30','plot_ylabel','non-dim');
        saveas(gcf,[f,'.swarm-C30.png']);
      end
      str.say(' --- Getting datastorage:')
      d=datastorage(...
        'start',v.start,...
        'stop', v.stop,...
        'debug',v.debug...
      ).init(p).init(ps);
      str.say(' --- Getting gravity:')
      g =d.data_get_scalar( p.dataname.append_field_leaf('signal'));
      gs=d.data_get_scalar(ps.dataname.append_field_leaf('signal'));
      str.say(' --- Time domain is:')
      disp([g.t_masked,gs.t_masked]);
      str.say(' --- Applying',v.mask,'spatial mask:')
      gm=g.set_lmax(40).spatial_mask(v.mask);
      str.say(' --- Smoothing, computing geoid and cumulative RMS:')
      [~,cts]=gm.scale(v.smoothing,'gauss').scale(v.functional,'functional').cumdrms;
      str.say(' --- Plotting the time series:')
      plotting.figure(varargin{:});
      cts.plot;
      plotting.enforce(....
        'plot_title',cts.labels{1},...
        'plot_ylabel',gravity.functional_label(v.functional),...
        'plot_legend_location','none'...
      );
      saveas(gcf,strjoin({f,v.functional,'quality','png'},'.'))
      str.say(' --- Exporting  the time series:')
      cts.export(strjoin({f,v.functional,'quality','dat'},'.'),'ascii')
    end
    function d=paper(varargin)
%             'swarm.sh.gswarm.rl01.catchments';...
%             'swarm.sh.gswarm.rl01.land';...
%             'swarm.sh.gswarm.rl01.lowdeg';...
%             'swarm.sh.gswarm.rl01.quality';...
%             'swarm.sh.gswarm.rl01.individual';...
%             'swarm.sh.gswarm.rl01.quality';...
%             'swarm.sh.gswarm.rl01.quality-land';...
%             'swarm.sh.gswarm.rl01.quality-land-deg40';...
%             'swarm.sh.gswarm.rl01.smooth0';...
%             'swarm.sh.gswarm.rl01.smooth300';...
%             'swarm.sh.gswarm.rl01.smooth750';...
%             'swarm.sh.gswarm.rl01.smooth1500';...
%             'swarm.sh.gswarm.rl01.smooth-ocean1000';...
%             'swarm.sh.gswarm.rl01.smooth-ocean1500';...
%             'swarm.sh.gswarm.rl01.smooth-ocean3000';...
%             'swarm.sh.gswarm.rl01.smooth-ocean5000';...
      v=varargs.wrap('sources',{....
        {...
          'debug',                   true, @logical;...
          'start', datetime('2013-09-01'), @isdatestime;...
          'stop',  datetime('2019-06-30'), @isdatestime;...
          'force',                  false, @islogical;...
          'plot_save_data',          true, @(i) islogical(i) || ischar(i);... %can be 'force'
          'plot_dir',fullfile(getenv('HOME'),'cloud','Work','articles','2019-04.gswarm'), @ischar;...
          'plot_type',                 '', @ischar;...
          'products',{...
            'swarm.sh.gswarm.rl01.individual';...
          }, @iscellstr;
        },... 
      },varargin{:});
      switch v.plot_type
      case ''
        args={'plot_dir',v.plot_dir,'force',v.force,'plot_save_data',v.plot_save_data};
        d=datastorage(...
          'start',v.start,...
          'stop', v.stop,...
          'inclusive',true,...
          'debug',v.debug...
        );
        for i=1:numel(v.products)
          d=d.init(v.products{i},args{:});
        end
      case 'TN11'
        gravity.graceC20('mode','model-plot')
        saveas(gcf,fullfile(v.plot_dir,'figures','TN11.png'))
      case 'TN11-parameters'
        gravity.graceC20('mode','model-list')
      case 'deepoceanmask'
        l=4;
        a=simplegrid.unit(l*90,l*45);
        a=a.spatial_mask('deep ocean');
        c=gray(l*16);
        c=flipud(c(l*4:l*12,:));
        plotting.figure;
        a.imagesc;
        colormap(c)
        colorbar off
        plotting.enforce('plot_legend_location','none')
        saveas(gcf,fullfile(v.plot_dir,'figures','deepoceanmask.png'))
      case 'graceC30'
        d=datastorage(...
          'start',datetime('2002-04-01'),...
          'stop', datetime('2017-06-30'),...
          'debug',v.debug...
        ).init('grace.sh.rl06.csr'...
        ).init('grace.sh.rl06.csr.pd.ts'...
        );
        %make pretty plot
        plotting.figure(varargin{:});
        grs=d.data_get_scalar('grace_sh_rl06_csr/signal').ts_C(3,0).add_expl_gaps(days(45)).plot('line',{'-o'});
        grm=d.data_get_scalar('grace_sh_rl06_csr_pd_ts' ).ts_C(3,0).plot;
        plotting.enforce('plot_legend',{...
          str.rep(grs.legend{1},'C3,0','data');...
          str.rep(grm.legend{1},'C3,0','model');...
        },'plot_title','','plot_ylabel','C_{3,0} [ ]');
        saveas(gcf,fullfile(v.plot_dir,'figures','graceC30.png'))
      end
    end
    function d=precombval(varargin)
      %NOTICE: this method produced the plots in ~/data/gswarm/analyses/<date>-precombval
      %        needed to produce the pre-combination report in ~/data/gswarm/analyses/<date>-precombval/report
      v=varargs.wrap('sources',{....
        {...
          'date',      datestr(now,'yyyy-mm-dd'), @isdatetime;...
          'force',     false,                     @islogical;... %this affects datastorage.init
        },... 
      },varargin{:});
      d=gswarm.production(...
        'start',     datetime('2016-01-01'),...
        'stop',      dateshift(datetime('now'),'end','month'),...
        'inclusive', true,...     %this can be false, because the GRACE data is only used to derive grace.sh.rl06.csr.ld.ts, separately 
        'force',     v.force,...  %NOTICE: grace.sh.rl06.csr.ld.ts has metadata never_force set as true (usually!) so this as false will only reload the Swarm individual models
        'products',  {...
          'swarm.sh.gswarm.rl01.individual';...
        },...
        'datadir',file.unresolve_home(['~/data/gswarm/analyses/',v.date,'-precombval/'])...
      );
    end
    function d=validation(varargin)
      %NOTICE: this method produced the plots in ~/data/gswarm/analyses/<date>
      %        needed to produce the signal content report on ~/data/gswarm/analyses/<date>/report (needs to be revised)
      v=varargs.wrap('sources',{....
        {...
          'date',      datestr(now,'yyyy-mm-dd'), @isdatetime;...
          'force',     false,                     @islogical;... %this affects datastorage.init
        },... 
      },varargin{:});
      d=gswarm.production(...
        'start',     datetime('2013-11-01'),...
        'stop',      dateshift(datetime('now'),'end','month'),...
        'inclusive', false,...    %this can be false, because the GRACE data is only used to derive grace.sh.rl06.csr.ld.ts, separately 
        'products',  {...
          'gswarm.swarm.validation.catchments';...
          'gswarm.swarm.validation.pardecomp';...
          'gswarm.swarm.validation.land';...
          'gswarm.swarm.validation.ocean';...
          'gswarm.swarm.validation.deepocean';...
          'gswarm.swarm.validation.smoothed';...
          'gswarm.swarm.validation.unsmoothed';...
        },...
        'datadir',file.unresolve_home(['~/data/gswarm/analyses/',v.date,'/'])...
      );
    end
    function d=SDQW9(varargin)
      v=varargs.wrap('sources',{....
        {...
          'date',      datestr(now,'yyyy-mm-dd'), @isdatetime;...
          'force',     false,                     @islogical;... %this affects datastorage.init
        },... 
      },varargin{:});
      d=gswarm.production(...
        'start',     datetime('2016-01-01'),...
        'stop',      datetime(v.date),...
        'inclusive', false,...    %this can be false, because the GRACE data is only used to derive grace.sh.rl06.csr.ld.ts, separately 
        'products',  {...
          'gswarm.swarm.validation.catchments';...
          'gswarm.swarm.validation.pardecomp';...
          'gswarm.swarm.validation.land';...
          'gswarm.swarm.validation.ocean';...
          'gswarm.swarm.validation.deepocean';...
          'gswarm.swarm.validation.smoothed';...
          'gswarm.swarm.validation.unsmoothed';...
        },...
        'datadir',file.unresolve_home(['~/data/gswarm/analyses/',v.date,'/'])...
      );
    end
    function d=production(varargin)
      %NOTICE: this method expects some input arguments, notably:
      % - products (cellstr)
      % - datadir (char)
      v=varargs.wrap('sources',{....
        {...
          'stop',      now,                       @isdatetime;...
          'start',     datetime('2013-11-01'),    @isdatetime;...
          'debug',     true,                      @islogical;...
          'inclusive', true,                      @islogical;...
          'force',     false,                     @islogical;... %this affects datastorage.init
          'force_d',   false,                     @islogical;... %this affects load(datafilename,'d') if force is false
          'products',  {},                        @iscellstr;...
          'datadir',   '',                        @ischar;...
        },... 
      },varargin{:});
      %check if all needed arguments are available
      assert(~isempty(v.products),'need argument "products"')
      assert(~isempty(v.datadir ),'need argument "datadir"' )
      %create vector with relevant data product
      p=cellfun(@(i) dataproduct(i,'plot_dir',v.datadir),v.products,'UniformOutput',false);
      %define data filename
      datafilename=fullfile(v.datadir,'d.mat');
      %load data if already available
      if file.exist(datafilename) && ~v.force_d && ~v.force
        str.say('Loading analysis data from',datafilename)
        load(datafilename,'d')
      else
        d=datastorage('debug',v.debug,'start',v.start,'stop',v.stop,'inclusive',v.inclusive);
        for i=1:p{1}.nr_sources
          d=d.init(p{1}.sources(i),'force',v.force);
        end
        file.ensuredir(datafilename,true);
        save(datafilename,'d')
      end
      %plot it
      tstart=tic;
      for i=1:numel(p)
        d.init(p{i});
        disp(['Finished plotting product ',p{i}.str,' after ',time.str(toc(tstart)),' elapsed'])
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

%NOTICE: this was an old (and rapidly-evolving) version that is now depracated
%     function production
%       %parameters
%       recompute=false;
%       
%       
%       %definitions
% %       datafilename=file.unresolve_home('~/data/gswarm/analyses/2018-06-24/d.mat');
% %       p=dataproduct('gswarm.swarm.all.res.plots','plot_dir',fileparts(datafilename));
% 
% %       datafilename=file.unresolve_home('~/data/gswarm/analyses/2018-07-04/d.mat');
% %       p=dataproduct('gswarm.swarm.all.res-unsmoothed.plots','plot_dir',fileparts(datafilename));
% 
% %       datafilename=file.unresolve_home('~/data/gswarm/analyses/2018-08-16/d.mat');
% %       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...
% 
% %       datafilename=file.unresolve_home('~/data/gswarm/analyses/2018-10-04/d.mat');
% %       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...   
% %         'gswarm.swarm.all.TN-03_2.catchments',...
% %         'gswarm.swarm.all.TN-03_2.unsmoothed',...
% %         'gswarm.swarm.all.TN-03_2.smoothed'...
% %       },'UniformOutput',false);
% 
% %       datafilename=file.unresolve_home('~/data/gswarm/analyses/2018-11-19/d.mat');
% %       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...   
% %         'gswarm.swarm.all.TN-03_2.ocean',...
% %         'gswarm.swarm.all.TN-03_2.smoothed',...
% %         'gswarm.swarm.all.TN-03_2.land',...
% %         'gswarm.swarm.all.TN-03_2.ocean',...
% %       },'UniformOutput',false);
% %       plot_stop=datetime('2017-07-31');
% %       datadir=file.unresolve_home('~/data/gswarm/analyses/2019-02-17/');
% 
% %       datadir=file.unresolve_home('~/data/gswarm/analyses/2019-04-05/');
% %       p=cellfun(@(i) dataproduct(i,'plot_dir',datadir),{...  
% %         'comp.grace.swarm.catchments',...
% %         'swarm.sh.gswarm.rl01.land',...
% %         'swarm.sh.gswarm.rl01.lowdeg',...
% %         'swarm.sh.gswarm.rl01.pardecomp',...
% %         'swarm.sh.gswarm.rl01.quality',...
% %       },'UniformOutput',false);
% %       start=datetime('2002-04-01');
% %       stop =datetime('2018-12-31');
% %       inclusive=true;
% 
% %       datadir=file.unresolve_home('~/data/gswarm/analyses/2019-04-15/');
% %       p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datadir)),{...  
% %         'gswarm.swarm.all.TN-03.unsmoothed',...
% %         'gswarm.swarm.all.TN-03.ocean',...
% %         'gswarm.swarm.all.TN-03.smoothed',...
% %         'gswarm.swarm.all.TN-03.land',...
% %         'gswarm.swarm.all.TN-03.catchments',...
% %         'gswarm.swarm.all.TN-03.pardecomp',...
% %       },'UniformOutput',false);
% %       start=datetime('2014-06-01');
% %       stop =datetime('2017-07-31');
% %       inclusive=false;
% %       %save version numbers into latex table
% %       file.strsave(fullfile(datadir,'versions.tex'),gswarm.sourcelatextable(p{1}));
% 
% %       datadir=file.unresolve_home('~/data/gswarm/analyses/2019-05-06/');
% %       p=cellfun(@(i) dataproduct(i,'plot_dir',datadir),{...  
% %         'swarm.sh.gswarm.rl01.catchments',...
% %         'swarm.sh.gswarm.rl01.land',...
% %         'swarm.sh.gswarm.rl01.quality',...
% %       },'UniformOutput',false);
% %       start=datetime('2012-05-01');
% %       stop =datetime('2018-12-31');
% %       inclusive=true;
% 
%       datadir=file.unresolve_home('~/data/gswarm/analyses/2019-05-26/');
%       p=cellfun(@(i) dataproduct(i,'plot_dir',datadir),{...  
%         'swarm.sh.gswarm.rl01.quality-check',...
%       },'UniformOutput',false);
%       start=datetime('2013-11-01');
%       stop =datetime('2019-03-31');
%       inclusive=true;
%     
%       if numel(p)==1
%         datafilename=fullfile(datadir,[p{1}.str,'.mat']);
%       else
%         datafilename=fullfile(datadir,'d.mat');
%       end
%       
%       %TODO: There's problem with the handling of model errors!!!
%       %TODO: can't have GRACE data and Swarm data all the way until Sep 2018
%       
%       %load data if already available
%       if file.exist(datafilename) && ~recompute
%         str.say('Loading analysis data from',datafilename)
%         load(datafilename,'d')
%       else
%         d=datastorage('debug',true,'start',start,'stop',stop,'inclusive',inclusive);
%         for i=1:p{1}.nr_sources
%           d=d.init(p{1}.sources(i),'force',recompute);
%         end
%         file.ensuredir(datafilename,true);
%         save(datafilename,'d')
%       end
%       %plot it
%       for i=1:numel(p)
%         d.init(p{i});
%       end
%       
%     end
