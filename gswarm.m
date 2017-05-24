classdef gswarm
  methods(Static)
    function obj=load_models(obj,product,varargin)
      %retrieve relevant parameters
      model_types       =product.mdget('model_types');
      indir             =product.mdget('import_dir');
      wilcarded_filename=product.mdget('wilcarded_filename');
      model_format      =product.mdget('model_format');
      date_parser       =str2func(product.mdget('date_parser'));
      max_degree        =product.mdget('max_degree');
      use_GRACE_C20     =product.mdget('use_GRACE_C20');
      delete_C00        =product.mdget('delete_C00');
      delete_C20        =product.mdget('delete_C20');
      static_field      =product.mdget('static_model');
      %load all available data 
      [s,e]=gravity.load_dir(indir,model_format,date_parser,...
        'wilcarded_filename',wilcarded_filename,...
        'start',obj.start,...
        'stop',obj.stop,...
        'descriptor',product.name...
      );
      %enforce consistent GM and R
      s=s.setGM(gravity.default_list.GM);
      s=s.setR( gravity.default_list.R );
      e=e.setGM(gravity.default_list.GM);
      e=e.setR( gravity.default_list.R );
      %set maximum degree (if requested)
      if max_degree>0
        s.lmax=max_degree;
        e.lmax=max_degree;
      end
      %set C20 coefficient
      if use_GRACE_C20
        disp([product.str,': use_GRACE_C20'])
        %some sanity
        if strcmpi(func2str(date_parser),'static')
          error([mfilename,': there''s no point in replacing GRACE C20 coefficients in a static model.'])
        end
        %get C20 timeseries, interpolated to current time domain
        c20=gravity.graceC20.interp(s.t);
%         figure
%         plot(c20.x_masked,c20.y_masked([],1),'x-','MarkerSize',10,'LineWidth',4), hold on
%         plot(c20.x,spline(c20.x_masked,c20.y_masked([],1),c20.x),'o-','MarkerSize',10,'LineWidth',2)
%         plot(c20.x,pchip(c20.x_masked,c20.y_masked([],1),c20.x),'+-','MarkerSize',10,'LineWidth',2)
%         legend('data','spline','pchip')
        %extrapolate in case there are gaps
        if c20.nr_gaps>0
          c20=c20.assign([...
            spline(c20.x_masked,c20.y_masked([],1),c20.x),...
            spline(c20.x_masked,c20.y_masked([],2),c20.x)...
          ],'t',c20.t,'mask',true(size(c20.x)));
        end
        for i=1:s.length
          s=s.setC(2,0,c20.y(i,1),s.t(i));
          e=e.setC(2,0,c20.y(i,2),s.t(i));
        end
%           figure
%           plot(c20.t,c20.y(:,1),'o-'), hold on
%           plot(m.t,m.C(2,0),'x-')
%           legend('GRACE',m.descriptor)
      end
      %remove C00 bias
      if delete_C00
        disp([product.str,': delete_C00'])
        for i=1:s.length
          s=s.setC(0,0,0);
          e=e.setC(0,0,0);
        end
      end
      %remove C20 bias
      if delete_C20
        disp([product.str,': delete_C20'])
        for i=1:s.length
          s=s.setC(2,0,0);
          e=e.setC(2,0,0);
        end
      end
      %remove static field (if requested)
      if ~strcmpi(static_field,'none')
        %load model (only if not already done)
        if isempty(dir([static_field,'.mat']))
          static=datastorage().init(static_field,'start',s.start,'stop',s.stop);
          save([static_field,'.mat'],'static')
        else
          load([static_field,'.mat'])
        end
        %subtract it
        s=s-static.data_get_scalar(datanames(static_field).set_field_path('signal'));
      end
      %propagate relevant data
      for i=1:numel(model_types)
        switch lower(model_types{i})
        case {'signal','sig','s'}
          obj=obj.data_set(product.dataname.set_field_path(model_types{i}),s);
        case {'error','err','e'}
          obj=obj.data_set(product.dataname.set_field_path(model_types{i}),e);
        otherwise
          error([mfilename,': unknown model type ''',model_types{i},'''.'])
        end
      end
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
      derived_quantity=product.mdget('derived_quantity');
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
    function obj=parametric_decomp(obj,product,varargin)
      %sanity
      assert(product.nr_sources==1,...
        [mfilename,': can only operate on a single source, not ',num2str(product.nr_sources),'.']...
      )
      %retrieve relevant parameters
      model_type =product.dataname.field_path{end};
      t_mod_f=product.mdget('t_mod_f','default',1);
      %sanity
      assert(numel(model_type)==1 && all(~strcmp(model_type,{'error','err','e'})),...
        [mfilename,': can not operate on models of type ''error''.']...
      )
      polynomial =ones(1,product.mdget('polyorder')+1);
      sinusoidal =time.num2duration(cell2mat(product.mdget('sin_period')),product.mdget('sin_period_unit'));
      %decompose
      s=obj.data_get_scalar(...
        [product.sources(1).name,'.signal']...
      ).parametric_decomposition(...
        't_mod_f',t_mod_f,...
        'polynomial',polynomial,...
        'sinusoidal',sinusoidal...
      );
      %propagate relevant data
      obj=obj.data_set(product,s);
    end
    function obj=plot_rms_ts(obj,product,varargin)
      obj.log('@','in','product',product)
      plot_functional =product.mdget('plot_functional');
      derived_quantity=product.mdget('derived_quantity');
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=product.plot_args(p,varargin{:});
      %retrieve plot elements (repetitive processing of parameters)
      e=obj.plot_elements(p,product);
      %loop over all data
      for t=1:numel(e.startlist)
        obj.log('@','iter','start',e.startlist(t),'stop',e.stoplist(t))
        %plot filename arguments
        filename_args=[product.file_args('plot'),{...
          'start',e.startlist(t),...
          'stop', e.stoplist(t),...
          'timestamp',true,...
          'remove_part','',...
          'prefix',p.Results.plot_file_prefix,...
          'suffix',e.suffix...
        }];
        %plot filename
        filename=product.dataname.file(filename_args{:});
        if isempty(dir(filename))
          %get the data for the current segment
          obj_curr=obj.trim('start',e.startlist(t),'stop',e.stoplist(t),'dn_list',e.sources);
          %make sure there is data 
          if any(cell2mat(obj_curr.vector_method_tr('all','nr_valid'))>1)
%             %remove C00 and C20
%             obj_curr=obj_curr.data_set('all',obj_curr.vector_method_tr('all','setC',[0 2],[0 0],[0 0]));
            %need all data to be in the same time domain
            assert(obj_curr.isteq('all'),...
              [mfilename,': the time domain of the input data is not in agreement.'])
            %build data array
            bar_y=zeros(max(obj_curr.length('all')),numel(e.sources));
            bar_t=zeros(max(obj_curr.length('all')),numel(e.sources));
            for di=1:numel(e.sources)
              tmp=obj_curr.data_get_scalar(e.sources{di}).scale(plot_functional,'functional').(derived_quantity);
              bar_y(:,di)=tmp(:,end);
              bar_t(:,di)=datenum(obj_curr.data_get_scalar(e.sources{1}).t);
            end
            %filter out nans
            good_idx=~all(isnan(bar_y),2);
            %plot it
            figure('visible',p.Results.plot_visible);
            h=bar(bar_t(good_idx,:),bar_y(good_idx,:),'EdgeColor','none');
            %enforce plot preferences
            product.enforce_plot
            %build plot annotation structure
            bh=cell(size(e.sources));
            for di=1:numel(e.sources)
              bh{di}=struct(...
                'mask',obj_curr.data_get_scalar(e.sources{di}).mask,...
                'y_mean',0,...
                'handle',h(di),...
                'title',e.sources{di}.str,...
                'xlabel','time',...
                'ylabel',plot_functional,...
                'y_units',gravity.functional_units(plot_functional),...
                'legend',{{e.sources{di}.name}}...
              );
            end
            %annotate plot
            obj.plot_annotate(bh,product,e.sources,varargin{:})
            %remove outline
            set(h,'edgecolor','none')
            %save this plot
            saveas(gcf,filename)
            obj.log('@','iter','plot saved',filename)
            % user feedback
            if strcmp(p.Results.plot_visible,'off')
              disp(['gswarm.plot_rms_ts: plotted ',product.name,' to file ',filename])
            end
          else
            disp(['gswarm.plot_rms_ts: not enough data to plot ',product.name,' to file ',filename,' (skipped)'])
          end
        else
          obj.log('@','iter','skipped',product,'file already exists',filename)
        end
      end
      obj.log('@','out')
    end
    function filelist=icgem(obj,product,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',@(i) isa(i,'dataproduct'));
      % parse it
      p.parse(product,varargin{:});
      %retrieve relevant data
      dat=obj.data_get_scalar(product);
      %call exporting routine
      filelist=dat.signal.icgem(...
        'prefix',product.name,...
        'path',  product.mdget('export_dir'),...
        'modelname',product.name...
      );
%         'error_obj',dat.error,...
    end
  end
end