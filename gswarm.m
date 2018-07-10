classdef gswarm
  methods(Static)
    function obj=load_models(obj,product,varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
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
        'static_model',   'none',@(i) ischar(i)...
      },product.metadata},varargin{:});
      %load all available data 
      [s,e]=gravity.load_dir(v.import_dir,v.model_format,str2func(v.date_parser),...
        'wilcarded_filename',v.wilcarded_filename,...
        'start',obj.start,...
        'stop',obj.stop,...
        'descriptor',product.name...
      );
      %apply model processing options
      [s,e]=gswarm.load_models_op('all',v,product,s,e);
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
          mod={'consistent_GM','consistent_R','max_degree','use_GRACE_C20','delete_C00','delete_C20','model_span','static_model'};
        case 'all'
          if exist('err','var')
            [mod,err]=gswarm.load_models_op(gswarm.load_models_op('mode_list'),v,product,mod,err);
          else
            mod=gswarm.load_models_op(gswarm.load_models_op('mode_list'),v,product,mod);
          end
        case 'consistent_GM'
          mod=mod.setGM(gravity.parameters('GM','value'));
          if exist('err','var')
            err=err.setGM(gravity.parameters('GM','value'));
          end
          show_msg=true;
        case 'consistent_R'
          mod=mod.setR( gravity.parameters('R' ,'value'));
          if exist('err','var')
            err=err.setR( gravity.parameters('R' ,'value'));
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
          if v.isparameter('use_GRACE_C20') && v.use_GRACE_C20
            %some sanity
            if strcmpi(v.date_parser,'static')
              error([mfilename,': there''s no point in replacing GRACE C20 coefficients in a static model.'])
            end
            %get C20 timeseries, interpolated to current time domain
            c20=gravity.graceC20.interp(mod.t);
  %         figure
  %         plot(c20.x_masked,c20.y_masked([],1),'x-','MarkerSize',10,'LineWidth',4), hold on
  %         plot(c20.x,spline(c20.x_masked,c20.y_masked([],1),c20.x),'o-','MarkerSize',10,'LineWidth',2)
  %         plot(c20.x,pchip(c20.x_masked,c20.y_masked([],1),c20.x),'+-','MarkerSize',10,'LineWidth',2)
  %         legend('data','spline','pchip')
            %extrapolate in case there are gaps
            if c20.nr_gaps>0
              c20=c20.assign([...
                interp1(c20.x_masked,c20.y_masked([],1),c20.x),...
                interp1(c20.x_masked,c20.y_masked([],2),c20.x)...
              ],'t',c20.t,'mask',true(size(c20.x)));
            end
            for i=1:mod.length
              mod=mod.setC(2,0,c20.y(i,1),mod.t(i));
              err=err.setC(2,0,c20.y(i,2),mod.t(i));
            end
  %           figure
  %           plot(c20.t,c20.y(:,1),'o-'), hold on
  %           plot(m.t,m.C(2,0),'x-')
  %           legend('GRACE',m.descriptor)
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
          msg={};
          %append extremeties, if requested
          if v.isparameter('model_span') && v.model_span{1}~=time.zero_date && v.model_span{1}<mod.start
            mod=mod.append(gravity.nan(mod.lmax,'t',v.model_span{1},'R',mod.R,'GM',mod.GM));
            if exist('err','var')
              err=mod.append(gravity.nan(err.lmax,'t',v.model_span{1},'R',err.R,'GM',err.GM));
            end
            msg{end+1}=['from ',datestr(v.model_span{1})];
            show_msg=true;
          end
          if v.isparameter('model_span') && v.model_span{2}~=time.inf_date && v.model_span{2}>mod.stop
            mod=mod.append(gravity.nan(mod.lmax,'t',v.model_span{2},'R',mod.R,'GM',mod.GM));
            if exist('err','var')
              err=mod.append(gravity.nan(err.lmax,'t',v.model_span{2},'R',err.R,'GM',err.GM));
            end
            msg{end+1}=['to ',datestr(v.model_span{2})];
            show_msg=true;
          end
          msg=strjoin(msg,' ');
        case 'static_model'
          %remove static field (if requested)
          if v.isparameter('static_model') && ~strcmpi(v.static_model,'none')
            %load model (only if not already done)
            if isempty(dir([v.static_model,'.mat']))
              static=datastorage().init(v.static_model);
              save([v.static_model,'.mat'],'static')
            else
              load([v.static_model,'.mat'])
            end
            %reduce to gravity class
            static=static.data_get_scalar(datanames(v.static_model).set_field_path('signal'));
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
    %% plots
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
            %need all data to be in the same time domain
            obj_curr=obj_curr.interp('all');
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
                'xlabel','',...
                'ylabel',plot_functional,...
                'y_units',gravity.functional_units(plot_functional),...
                'legend',{{e.sources{di}.name}}...
              );
            end
            %annotate plot
            obj.plot_annotate(bh,product,datanames.transmute(e.sources),varargin{:})
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
    end %better use plot_spatial_stats, more general
    function out=plot_ops(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
          'plot_min_degree',        2 ,@(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',      inf ,@(i) isnumeric(i) && isscalar(i);...
          'plot_smoothing_degree', [] ,@(i) isnumeric(i) && isscalar(i);...
          'plot_smoothing_method', '' ,@(i) ischar(i);...
        },plotting.default,...
        product.plot_args...
      },varargin{:});
      %collect the models 
      out.sources=cell(product.nr_sources,1);out.source_datanames=cell(size(out.sources));mc=0;
      for i=1:product.nr_sources
        dn_sources=obj.data_list(product.sources(i));
        for j=1:numel(dn_sources)
          %only plot the relevant model types
          if ~cells.isincluded(product.metadata.model_types,dn_sources{j}.field_path{end}); continue;end
          mc=mc+1;
          out.source_datanames(mc)=dn_sources(j);
          out.sources{mc}=obj.data_get_scalar(out.source_datanames{mc});
          %handle default value of plot_max_degree
          v.plot_max_degree=min([out.sources{mc}.lmax,v.plot_max_degree]);
        end
      end
      %enforce maximum degree
      out.sources=cellfun(@(i) i.set_lmax(v.plot_max_degree),out.sources,'UniformOutput',false);
      %enforce smoothing
      if ~isempty(v.plot_smoothing_degree) && ~isempty(v.plot_smoothing_method)
        out.sources=cellfun(@(i) i.scale(...
          v.plot_smoothing_degree,...
          v.plot_smoothing_method...
        ),out.sources,'UniformOutput',false);
        out.file_smooth=['.smooth-',v.plot_smoothing_method,'-',gravity.gauss_smoothing_name(v.plot_smoothing_degree)];
        out.title_smooth=[gravity.gauss_smoothing_name(v.plot_smoothing_degree),' ',...
                                gravity.smoothing_name(v.plot_smoothing_method),' smoothing'];
      else
        out.file_smooth='';
        out.title_smooth='';
      end
      %find reference product
      if product.ismdfield('stats_relative_to')
        for i=1:numel(out.source_datanames)
          if strfind(out.source_datanames{i}.str,product.metadata.stats_relative_to)
            out.mod_ref=out.sources{i};
            out.mod_ref_idx=i;
            out.mod_ref_name=out.source_datanames{out.mod_ref_idx};
            break
          end
        end
        %get the indexes of the sources to plot (i.e. not the reference)
        mod_idx=[1:out.mod_ref_idx-1,out.mod_ref_idx+1:product.nr_sources];
        %get product name difference between products to derive statistics from
        out.mod_names=cellfun(@(i) strjoin(i,' '),datanames.unique(out.source_datanames(mod_idx)),'UniformOutput',false);
        out.source_names(mod_idx)=cellfun(@(i) strrep(i,'_',' '),out.mod_names,'UniformOutput',false);
        out.source_names{out.mod_ref_idx}=upper(strtrim(strrep(str.clean(out.mod_ref_name.name,v.plot_title_suppress),'.',' ')));
        %reduce the models to plot
        out.mods=out.sources(mod_idx);
        %compute residual
        out.res=cellfun(@(i) out.mod_ref.interp(i.t)-i,out.mods,'UniformOutput',false);
        %title
        out.title_wrt=['wrt ',out.source_names{out.mod_ref_idx},' '];
        %get time domain common to all residuals
        out.t=simpletimeseries.t_mergev(out.res);
      else
        out.source_names=cellfun(...
          @(i) strtrim(strrep(str.clean(i,v.plot_title_suppress),'.',' ')),...
          cellfun(@(i) i.name,out.source_datanames,'UniformOutput',false),...
          'UniformOutput',false...
        );
        out.mod_names=out.source_names;
        %alias sources to mods and res
        out.mods=out.sources;
        out.res=out.sources;
        %patch remaining details
        out.title_wrt='';
        %get time domain common to all sources
        out.t=simpletimeseries.t_mergev(out.sources);
      end
      %filename particles
      out.file_root=product.file('plot','start',obj.start,'stop',obj.stop,'start_timestamp_only',false);
      out.file_deg=['deg',num2str(v.plot_min_degree),'-',num2str(v.plot_max_degree)];
      %legend
      out.source_legend_str=cellfun(@(i) strjoin(i,' '),datanames.unique(out.source_datanames),'UniformOutput',false);
      %title
      out.title_startstop=['(',datestr(out.sources{1}.t(1),'yyyy-mm'),' to ',datestr(out.sources{1}.t(end),'yyyy-mm'),')'];
      %start/stop
      [~,out.startlist,out.stoplist]=product.file('data',v.varargin{:},'start',obj.start,'stop',obj.stop);
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_low_degrees(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'degrees',    [2   3  ],@(i) isnumeric(i);...
        'orders',     [inf inf],@(i) isnumeric(i);...
        },...
        plotting.default,...
        rmfield(product.metadata,'sources')...
      },varargin{:});
      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      %convert degree-wise orders
      degrees_out=cell(size(v.degrees));
       orders_out=cell(size(v.orders));
      %need to convert keywords in orders ('inf'), if there
      v.orders=cell2mat(cells.num(v.orders));
      for i=1:numel(v.degrees)
        if ~isfinite(v.orders(i))
           orders_out{i}=-v.degrees(i):v.degrees(i);
          degrees_out{i}=ones(size(orders_out{i}))*v.degrees(i);
        else
           orders_out{i}=v.orders(i);
          degrees_out{i}=v.degrees(i);
        end
      end
      orders=cell2mat(orders_out);
      degrees=cell2mat(degrees_out);
      %loop over all requested degrees and orders
      for i=1:numel(degrees)
        d=degrees(i);
        o=orders(i);
        filename=strrep(out.file_root{1},'.png',['.C',num2str(d),',',num2str(o),'.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          %build legend string
          legend_str=cell(size(out.source_legend_str));
          legend_str_idx=true(size(legend_str));
          plotting.figure(v.varargin{:});
          %loop over all models
          for k=1:numel(out.sources)
            %plot the time series for this degree/order and model
            ts_now=out.sources{k}.ts_C(d,o);
            if isempty(ts_now) || ts_now.iszero
              legend_str_idx(k)=false;
              continue
            else
              ts_now.plot;
            end
            %add statistics to the legend (unless this is the produce from which the stats are derived)
            if k~=out.mod_ref_idx
              mod_ref_now=out.mod_ref.ts_C(d,o).interp(ts_now.t);
              stats=ts_now.stats2(mod_ref_now,'mode','struct','struct_fields',{'corrcoef','rmsdiff'},'period',seconds(inf));
              legend_str{k}=[out.source_legend_str{k},...
                ' corr=',num2str(stats.corrcoef,'%.2f'),...
                ', RMS{\Delta}=',num2str(stats.rmsdiff,'%.2g')];
            else
              legend_str{k}=[legend_str{k},' (reference)'];
            end
          end
          product.enforce_plot(varargin{:},...
            'plot_ylabel','[ ]',...
            'plot_legend',legend_str(legend_str_idx),...
            'plot_fontsize_legend',14,...
            'plot_title',['C',num2str(d),',',num2str(o)]...
          );
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_spatial_stats(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
          'plot_min_degree',         2, @(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',       inf, @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',   'geoid', @(i) gravity.isfunctional(i);...
          'derived_quantity','cumdrms', @(i) ismethod(gravity.unit(1),i);...
        },plotting.default,...
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
      %build filename
      filename=strrep(out.file_root{1},'.png',['.cumdrms',out.file_smooth,'.',strjoin(str.rep(out.source_names,' ','_'),'-'),'.png']);
      if isempty(dir(filename))
        %build data array
        bar_y=zeros(numel(out.t),numel(out.res));
        for di=1:numel(out.res)
          tmp=out.res{di}.interp(out.t).scale(v.plot_functional,'functional').(v.derived_quantity);
          bar_y(:,di)=tmp.y;
        end
        %plot it
        plotting.figure(v.varargin{:});
        h=bar(datenum(out.t),bar_y,'EdgeColor','none');
        %remove outline
        set(h,'edgecolor','none')
        %enforce it
        product.enforce_plot(v.varargin{:},...
          'plot_legend',out.mod_names,...
          'plot_ylabel',gravity.functional_label(v.plot_functional),...
          'plot_xdate',true,...
          'plot_xlimits',datenum([out.t(1),out.t(end)+days(1)]),...
          'plot_title',['Residual ',out.title_wrt,out.title_smooth]...
        );
        colormap jet
        saveas(gcf,filename)
        str.say('Plotted',filename)
      else
        disp(['NOTICE: plot already available: ',filename])
      end
      for i=1:numel(out.t)
        %gather models valid now
        dat    =cellfun(@(j) j.interp(out.t(i)),out.sources,'UniformOutput',false);
        dat_idx=cellfun(@(j) j.nr_valid>0,dat);
        %check if there's any data to plot
        if all(~dat_idx); continue;end
        %reduce data
        dat=dat(dat_idx);
        legend_str=out.source_names(dat_idx);
        %build filename
        filename=cells.scalar(product.file('plot',...
          'start',out.t(i),'stop',out.t(i),...
          'suffix',['drms',out.file_smooth,'.',strjoin(str.rep(out.source_names,' ','_'),'-')]...
        ),'get');
        if isempty(dir(filename))
          plotting.figure(v.varargin{:});
          for j=1:numel(dat)
            dat{j}.plot('mode','drms','functional','geoid');
          end
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',legend_str,...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_title',[datestr(out.t(i),'yyyy-mm'),' degree-RMS ',out.title_smooth]...
          );
          colormap jet
          saveas(gcf,filename)
          str.say('Plotted',filename)
        else
          disp(['NOTICE: plot already available: ',filename])
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_temporal_stats(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
          'plot_min_degree',     2 ,@(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',   inf ,@(i) isnumeric(i) && isscalar(i);...
          'plot_rmsdiff_caxis',[-inf inf] ,@(i) isnumeric(i) && isscalar(i);...
          'plot_temp_stat_list', {            'corrcoeff',           'rmsdiff',           'stddiff'},@(i) iscellstr(i); ...
          'plot_temp_stat_func', {               'nondim',             'geoid',             'geoid'},@(i) iscellstr(i); ...
          'plot_temp_stat_title',{'temporal corr. coeff.','temporal RMS\Delta','temporal STD\Delta'},@(i) iscellstr(i); ...
        },plotting.default,...
        product.plot_args...
      },varargin{:});
      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      %assemble title suffix
      title_suffix=[out.title_wrt,' ',out.title_startstop];
      if ~isempty(out.title_smooth);title_suffix=[title_suffix,newline,out.title_smooth];end
      %loop over all statistics
      for s=1:numel(v.plot_temp_stat_list)
        %maybe nothing is requested to be plotted
        if strcmp(v.plot_temp_stat_list{s},'none'); continue; end
        %loop over all sources
        for i=1:numel(out.mods);
          %build filename
          filename=strrep(out.file_root{1},'.png',['.',v.plot_temp_stat_list{s},'-triang',out.file_smooth,'.',out.source_names{i},'.png']);
          %plot only if not done yet
          if exist(filename,'file')==0
            plotting.figure(v.varargin{:});
            %compute the correlation coefficient between this model and mod_ref
            d=out.mod_ref.scale(v.plot_temp_stat_func{s},'functional').interp(out.mods{i}.t).stats2(...
              out.mods{i}.scale(v.plot_temp_stat_func{s},'functional'),...
              'mode','obj',...
              'struct_fields',v.plot_temp_stat_list(s),...
              'period',seconds(inf)...
            );
            %set y-units
            d.y_units(:)={gravity.functional_units(v.plot_temp_stat_func{s})};
            %plot it
            d.plot('method','triang');
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_caxis',v.(['plot_',v.plot_temp_stat_list{s},'_caxis']),...
              'plot_title',[out.mod_names{i},' ',v.plot_temp_stat_title{s},' ',title_suffix]...
            );
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
        %plot degree-mean corrcoeff
        filename=strrep(out.file_root{1},'.png',['.',v.plot_temp_stat_list{s},'-dmean',out.file_smooth,'.',strjoin(out.source_names,'-'),'.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          %init plot counter and data container
          d=zeros(numel(out.mods),out.mod_ref.lmax+1);
          %loop over all sources
          for i=1:numel(out.mods);
            %compute it
            d(i,:)=out.mod_ref.scale(v.plot_temp_stat_func{s},'functional').interp(out.mods{i}.t).stats2(...
                   out.mods{i}.scale(v.plot_temp_stat_func{s},'functional'),...
                   'mode','obj',...
                   'struct_fields',v.plot_temp_stat_list(s),...
                   'period',seconds(inf)...
            ).dmean;
          end
          %enforce minimum degree
          d=d(:,v.plot_min_degree+1:end);
          %filter out nans
          good_idx=~all(isnan(d),1);
          %plot it
          plotting.figure(v.varargin{:});
          bar(v.plot_min_degree:out.mod_ref.lmax,d(:,good_idx)','EdgeColor','none');
          colormap jet
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',out.mod_names,...
            'plot_ylabel',gravity.functional_label(v.plot_temp_stat_func{s}),...
            'plot_xlabel','SH degree',...
            'plot_xlimits',[v.plot_min_degree-1,v.plot_max_degree+1],...
            'plot_title',['degree-mean ',v.plot_temp_stat_title{s},' ',title_suffix]...
          );
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
%       %loop over all sources
%       for i=1:numel(out.mods);
%         %build filename
%         filename=strrep(out.file_root{1},'.png',['.rmsdiff-triang',out.file_smooth,'.',out.source_names{i},'.png']);
%         %plot only if not done yet
%         if exist(filename,'file')==0
%           plotting.figure(v.varargin{:});
%           %compute the correlation coefficient between this model and mod_ref
%           rmsdiff=out.mod_ref.scale('geoid','functional').interp(out.mods{i}.t).stats2(...
%              out.mods{i}.scale('geoid','functional'),'mode','obj','struct_fields',{'rmsdiff'},'period',seconds(inf)...
%           );
%           %set y-units
%           rmsdiff.y_units(:)={'m'};
%           %plot it
%           rmsdiff.plot('method','triang');
%           %enforce it
%           product.enforce_plot(v.varargin{:},...
%             'plot_caxis',v.plot_rmsdiff_caxis,...
%             'plot_title',[out.mod_names{i},' temporal RMS\Delta ',title_suffix]...
%           );
%           saveas(gcf,filename)
%           str.say('Plotted',filename)
%         end
%       end
%       %plot degree-mean corrcoeff
%       filename=strrep(out.file_root{1},'.png',['.rmsdiff-dmean',out.file_smooth,'.',strjoin(out.source_names,'-'),'.png']);
%       %plot only if not done yet
%       if exist(filename,'file')==0
%         %init plot counter and data container
%         d=zeros(numel(out.mods),out.mod_ref.lmax+1);
%         %loop over all sources
%         for i=1:numel(out.mods);
%           %compute it
%           d(i,:)=out.mod_ref.scale('geoid','functional').interp(out.mods{i}.t).stats2(...
%             out.mods{i}.scale('geoid','functional'),'mode','obj','struct_fields',{'rmsdiff'},'period',seconds(inf)...
%           ).dmean;
%         end
%         %enforce minimum degree
%         d=d(:,v.plot_min_degree+1:end);
%         %filter out nans
%         good_idx=~all(isnan(d),1);
%         %plot it
%         plotting.figure(v.varargin{:});
%         bar(v.plot_min_degree:out.mod_ref.lmax,d(:,good_idx)','EdgeColor','none');
%         colormap jet
%         %enforce it
%         product.enforce_plot(v.varargin{:},...
%           'plot_legend',out.mod_names,...
%           'plot_ylabel','geoid height [m]',...
%           'plot_xlabel','SH degree',...
%           'plot_xlimits',[v.plot_min_degree-1,v.plot_max_degree+1],...
%           'plot_title',['degree-mean temporal RMS\Delta ',title_suffix]...
%         );
%         saveas(gcf,filename)
%         str.say('Plotted',filename)
%       end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_stats(obj,product,varargin)
      obj=gswarm.plot_spatial_stats( obj,product,varargin{:});
      obj=gswarm.plot_temporal_stats(obj,product,varargin{:});
    end
    function obj=plot_parametric_decomposition(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
          't_mod_f',               1 ,@(i) isnumeric(i) && isscalar(i);...
          'plot_functional',  'eqwh', @(i) gravity.isfunctional(i);...
          'polyorder',                                     2 ,@(i) isnumeric(i) && isscalar(i);...
          'polynames',      {'constant','linear','quadratic'},@(i) iscellstr(i);...
          'plot_poly_range',[      inf,     inf,        inf],@(i) isnumeric(i);...
          'sin_period',              [      12,            6],@(i) isnumeric(i);...
          'sin_names',               {'yearly','semi-yearly'},@(i) iscellstr(i);...
          'plot_amplitude_range'     [inf,               inf],@(i) isnumeric(i);...
          'sin_period_unit',                        'months' ,@(i) ischar(i);...
        },plotting.default,...
        rmfield(product.metadata,'sources')...
      },varargin{:});
    
%tmp code
data_file='tmp_data.mat';
if exist(data_file,'file')    
  load(data_file)
else
%tmp code

      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      %loop over all models
      out.pd=cell(size(out.sources));
      for i=1:numel(out.sources)
        %compute parametric decompositions
        out.pd{i}=out.sources{i}.parametric_decomposition(...
          't_mod_f',v.t_mod_f,...
          'polynomial',0:v.polyorder,...
          'sinusoidal',time.num2duration(v.sin_period,v.sin_period_unit),...
          't0',seconds(obj.start-out.sources{i}.epoch)...
        );
        %assemble title suffix
        title_suffix=out.title_startstop;
        if ~isempty(out.title_smooth);title_suffix=[title_suffix,newline,out.title_smooth];end %#ok<AGROW>
        for j=0:v.polyorder
          par_name=['p',num2str(j)];
          %build filename
          filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.',par_name,out.file_smooth,'.',strrep(out.source_names{i},' ','-'),'.png']);
          %plot only if not done yet
          if exist(filename,'file')==0
            a=out.pd{i}.(par_name).scale(...
              v.plot_functional,'functional'...
            ).grid;
            %build cb titles
            if j==0
              cb_title=gravity.functional_label(v.plot_functional);
            else
              cb_title=[gravity.functional_names(v.plot_functional),' [',gravity.functional_units(v.plot_functional),...
                '/',v.sin_period_unit,'^',num2str(j),']'];
              %need to scale the time domain units
              time_scale_function=str2func(v.sin_period_unit);
              a=a.scale(simpletimeseries.timescale(time_scale_function(1))^j);
            end
            %plot it
            plotting.figure(v.varargin{:});
            a.imagesc(...
              'cb_title',cb_title...
            );
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_title',[v.polynames{j+1},' term for ',out.source_names{i},' ',title_suffix]...
            );
            if isfinite(v.plot_poly_range(j+1))
              caxis([-v.plot_poly_range(j+1),v.plot_poly_range(j+1)])
            end
            colormap jet
            cb.zero_white;
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
        for j=1:numel(v.sin_period)
          %build filename for amplitude
          filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.a',num2str(j),out.file_smooth,'.',strrep(out.source_names{i},' ','-'),'.png']);
          %plot only if not done yet
          if exist(filename,'file')==0
            %get the sine and cosine terms in the form of grids
            s=out.pd{i}.(['s',num2str(j)]).scale(v.plot_functional,'functional').grid;
            c=out.pd{i}.(['c',num2str(j)]).scale(v.plot_functional,'functional').grid;
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
              'plot_title',[v.sin_names{j},' amplitude for ',out.source_names{i},' ',title_suffix]...
            );
            if isfinite(v.plot_amplitude_range(j))
              caxis([0,v.plot_amplitude_range(j)])
            end
            colormap jet
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
          %build filename for phase
          filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.f',num2str(j),out.file_smooth,'.',strrep(out.source_names{i},' ','-'),'.png']);
          %plot only if not done yet
          if exist(filename,'file')==0
            %get the sine and cosine terms in the form of grids
            s=out.pd{i}.(['s',num2str(j)]).grid;
            c=out.pd{i}.(['c',num2str(j)]).grid;
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
              'plot_title',[v.sin_names{j},' phase for ',out.mod_names{i},' ',title_suffix]...
            );
            caxis([-v.sin_period(j),v.sin_period(j)]/2)
            colormap jet
            cb.zero_white;
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
      end
      
%tmp code
save(data_file,'out')
end
%tmp code
  
      for j=1:size(simplegrid.catchment_list,1)
        %build filename
        filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.',...
          strrep(simplegrid.catchment_list{j,1},' ','-'),out.file_smooth,...
          '.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          plotting.figure(v.varargin{:});
          legend_str=cell(1,numel(out.sources));
          for i=1:numel(out.sources)
            out.catch{i,j}=out.res{i}.scale(v.plot_functional,'functional').grid('Nlat',100,'Nlon',200).catchment(...
              simplegrid.catchment_list{j,1},...
              'parametric_decomposition',false,...
              v.varargin{:}...
            );
            legend_str{i}=out.source_names{i};
          end
          plotting.enforce(...
            v.varargin{:},...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_line_style',cells.deal({'o-'},size(out.sources)),...
            'plot_legend',legend_str...
          );
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      
      for i=1:numel(out.sources)
        %build filename
        filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.',...
          'std',out.file_smooth,'.',...
          strrep(out.source_names{i},' ','-'),'.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          plotting.figure(v.varargin{:});
          out.res{i}.scale(v.plot_functional,'functional').grid('Nlat',100,'Nlon',200).stats('mode','std').imagesc(...
            'boxes',simplegrid.catchment_list...
          );
          plotting.enforce(...
            v.varargin{:},...
            'plot_line_color','',...
            'plot_title',['temporal STD of ',out.mod_names{i}]...
          );
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      
    end
    %% utils
    function production
      %definitions
%       datafilename=file.unresolve('~/data/gswarm/analyses/2018-06-24/d.mat');
%       p=dataproduct('gswarm.swarm.all.res.plots','plot_dir',fileparts(datafilename));

      datafilename=file.unresolve('~/data/gswarm/analyses/2018-07-04/d.mat');
      p=dataproduct('gswarm.swarm.all.res-unsmoothed.plots','plot_dir',fileparts(datafilename));
      %load data if already available
      if exist(datafilename,'file')~=0
        load(datafilename)
      else
        d=datastorage('debug',true,'start',datetime('2013-12-01'),'stop',datetime('2016-12-31'));
        for i=1:p.nr_sources
          d=d.init(p.sources(i));
        end
        file.ensuredir(datafilename,true);
        save(datafilename,'d')
      end
      %plot it
      d.init(p);
    end
  end
end