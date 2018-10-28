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
        'static_model',   'none',@(i) ischar(i);...
        'overwrite_common_t',false,@(i) islogical(i);...
      },product.metadata},varargin{:});
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
            mod=mod.setC(2,0,c20.y(:,1),mod.t);
            err=err.setC(2,0,c20.y(:,2),mod.t);
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
            %load model (only if not already done)
            if isempty(dir([v.static_model,'.mat']))
              static=datastorage().init(v.static_model);
              %reduce to gravity class
              static=static.data_get_scalar(datanames(v.static_model).set_field_path('signal'));
              %save it
              save([v.static_model,'.mat'],'static')
            else
              load([v.static_model,'.mat'])
            end
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
      v=varargs.wrap('sources',{product.metadata},varargin{:});
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
          'plot_min_degree',        2 ,@(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',      inf ,@(i) isnumeric(i) && isscalar(i);...
          'plot_smoothing_degree', [] ,@(i) isnumeric(i) && isscalar(i);...
          'plot_smoothing_method', '' ,@(i) ischar(i);...
        },...
        product.plot_args...
      },varargin{:});
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
            if ~cells.isincluded(product.metadata.model_types,dn_sources{j}.field_path{end}); continue;end
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
        if isscalar(v.plot_smoothing_degree)
          v.plot_smoothing_degree=v.plot_smoothing_degree*ones(size(out.source.signal));
        else
          assert(numel(v.plot_smoothing_degree)==out.source.n,...
            'If plot_smoothing_degree is a vector, it needs to have the same number of elemets as sources')
        end
        for i=1:out.source.n
          out.source.signal{i}=out.source.signal{i}.scale(...
            v.plot_smoothing_degree(i),...
            v.plot_smoothing_method...
          );
        end
        out.file_smooth=['.smooth-',v.plot_smoothing_method,'-',gravity.gauss_smoothing_name(v.plot_smoothing_degree(end))];
        out.title_smooth=[gravity.gauss_smoothing_name(v.plot_smoothing_degree(end)),' ',...
                                gravity.smoothing_name(v.plot_smoothing_method),' smoothing'];
      else
        out.file_smooth='';
        out.title_smooth='';
      end
      %find reference product
      if product.ismdfield('stats_relative_to')
        for i=1:out.source.n
          str.say(out.source.datanames{i}.str,'?=',[product.metadata.stats_relative_to,'/'])
          if strfind(out.source.datanames{i}.str,[product.metadata.stats_relative_to,'/'])
            out.mod.ref=out.source.signal{i};
            out.source.ref_idx=i;
            out.mod.ref_name=out.source.datanames{out.source.ref_idx};
            break
          end
        end
        assert(isfield(out.source,'ref_idx'),['None of the sources of product ',product.str,...
          ' match ''stats_relative_to'' equal to ''',product.metadata.stats_relative_to,'''.'])
        %get the indexes of the sources to plot (i.e. not the reference)
        out.source.mod_idx=[1:out.source.ref_idx-1,out.source.ref_idx+1:out.source.n];
        %get product name difference between products to derive statistics from
        out.mod.names=cellfun(@(i) strjoin(i,' '),datanames.unique(out.source.datanames(out.source.mod_idx)),'UniformOutput',false);
        out.source.names(out.source.mod_idx)=cellfun(@(i) strrep(i,'_',' '),out.mod.names,'UniformOutput',false);
        out.source.names{out.source.ref_idx}=upper(strtrim(strrep(str.clean(out.mod.ref_name.name,v.plot_title_suppress),'.',' ')));
        %reduce the models to plot
        out.mod.dat=out.source.signal(out.source.mod_idx);
        %compute residual
        out.mod.res=cellfun(@(i) out.mod.ref.interp(i.t)-i,out.mod.dat,'UniformOutput',false);
        %title
        out.title_wrt=['wrt ',out.source.names{out.source.ref_idx},' '];
        %get time domain common to all residuals
        out.t=simpletimeseries.t_mergev(out.mod.res);
      else
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
      end
      %filename particles
      out.file_root=product.file('plot','start',obj.start,'stop',obj.stop,'start_timestamp_only',false);
      out.file_deg=['deg',num2str(v.plot_min_degree),'-',num2str(v.plot_max_degree)];
      %legend
      out.source.legend_str=cellfun(@(i) strjoin(i,' '),datanames.unique(out.source.datanames),'UniformOutput',false);
      %patch empty legend entries (this expects there to be only one empty legend entry) 
      if any(cells.isempty(out.source.legend_str))
        out.source.legend_str(cells.isempty(out.source.legend_str))={strjoin(datanames.common(out.source.datanames),' ')};
      end
      %title
      out.title_startstop=['(',datestr(out.source.signal{1}.t(1),'yyyy-mm'),' to ',datestr(out.source.signal{1}.t(end),'yyyy-mm'),')'];
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
        },...
        rmfield(product.metadata,'sources')...
      },varargin{:});
      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      %get collection of degrees/orders
      [degrees,orders]=gravity.resolve_degrees_orders(v.varargin{:});
      %loop over all requested degrees and orders
      for i=1:numel(degrees)
        d=degrees(i);
        o=orders(i);
        filename=strrep(out.file_root{1},'.png',['.C',num2str(d),',',num2str(o),'.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          %build legend string
          legend_str=cell(1,out.source.n);
          legend_str_idx=true(size(legend_str));
          plotting.figure(v.varargin{:});
          %loop over all models
          for k=1:out.source.n
            %plot the time series for this degree/order and model
            ts_now=out.source.signal{k}.ts_C(d,o);
            if isempty(ts_now) || ts_now.iszero
              legend_str_idx(k)=false;
              continue
            else
              ts_now.plot;
            end
            %add statistics to the legend (unless this is the produce from which the stats are derived)
            if k~=out.source.ref_idx
              mod_ref_now=out.mod.ref.ts_C(d,o).interp(ts_now.t);
              if all(isnan(mod_ref_now.y(:)))
                stats.corrcoef='NaN';stats.rmsdiff='NaN';
              else
                stats=ts_now.stats2(mod_ref_now,'mode','struct','struct_fields',{'corrcoef','rmsdiff'},'period',seconds(inf));
              end
                legend_str{k}=[out.source.legend_str{k},...
                  ' corr=',num2str(stats.corrcoef,'%.2f'),...
                  ', RMS{\Delta}=',num2str(stats.rmsdiff,'%.2g')];
            else
              legend_str{k}=[out.source.legend_str{k},' (reference)'];
            end
          end
          product.enforce_plot(varargin{:},...
            'plot_ylabel','[ ]',...
            'plot_legend',legend_str(legend_str_idx),...
            'plot_fontsize_legend',14,...
            'plot_title',['C',num2str(d),',',num2str(o)]...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
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
          'plot_derived_quantity',          'cumdrms', @(i) ismethod(gravity.unit(1),i);...
          'plot_spatial_stat_list',{'diff','monthly'}, @iscellstr; ... TODO: corr
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
        %build filename
        filename=strrep(out.file_root{1},'.png',['.cumdrms',out.file_smooth,'.',strjoin(str.rep(out.source.names,' ','_'),'-'),'.png']);
        if isempty(dir(filename))
          %build data array
          y=zeros(numel(out.t),numel(out.mod.res));
          for di=1:numel(out.mod.res)
            tmp=out.mod.res{di}.interp(out.t).scale(v.plot_functional,'functional').(v.plot_derived_quantity);
            y(:,di)=tmp.y;
          end
          %plot it
          plotting.figure(v.varargin{:});
          switch v.plot_type
          case 'bar'
            bar(out.t,y,'EdgeColor','none');
          case 'line'
            plot(out.t,y,'Marker','o');
          end
          %deal with legend stats
          if v.plot_show_legend_stats
            for i=1:numel(out.mod.res)
              out.mod.names{i}=[out.mod.names{i},...
                   ' ',num2str(mean(y(~isnan(y(:,i)),i)),2),...
               ' +/- ',num2str( std(y(~isnan(y(:,i)),i)),2),...
            ' \Sigma=',num2str( sum(y(~isnan(y(:,i)),i)),2)];
            end
          end
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',out.mod.names,...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_xdate',true,...
            'plot_xlimits',[out.t(1),out.t(end)+days(1)],...
            'plot_title',['Residual ',out.title_wrt,out.title_smooth]...
          );
          colormap jet
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        else
          disp(['NOTICE: plot already available: ',filename])
        end
      end
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
            'suffix',['drms',out.file_smooth,'.',strjoin(str.rep(out.source.names,' ','_'),'-')]...
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
              'plot_colormap','jet',...
              'plot_title',[datestr(out.t(i),'yyyy-mm'),' degree-RMS ',out.title_smooth]...
            );
            if v.plot_monthly_error
              %get previous lines
              lines_before=findobj(gca,'Type','line');
              %plot errors
              for j=1:numel(dat)
                dat_error{j}.plot('mode','drms','functional','geoid','line','--');
              end
              %get all lines
              lines_all=findobj(gca,'Type','line');
              %set consistent line clours
              for j=1:numel(dat)
                set(lines_all(j),'Color',get(lines_before(j),'Color'),'LineWidth',v.plot_line_width)
              end
              axis auto
              title([datestr(out.t(i),'yyyy-mm'),' degree-RMS ',out.title_smooth])
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
          'plot_min_degree',     2 ,@(i) isnumeric(i) && isscalar(i);...
          'plot_max_degree',   inf ,@(i) isnumeric(i) && isscalar(i);...
          'plot_rmsdiff_caxis',[-inf inf] ,@(i) isnumeric(i) && isscalar(i);...
          'plot_temp_stat_list', {            'corrcoeff',           'rmsdiff',           'stddiff'},@(i) iscellstr(i); ...
          'plot_temp_stat_func', {               'nondim',             'geoid',             'geoid'},@(i) iscellstr(i); ...
          'plot_temp_stat_title',{'temporal corr. coeff.','temporal RMS\Delta','temporal STD\Delta'},@(i) iscellstr(i); ...
          'plot_type',         'line', @(i) cells.isincluded(i,{'line','bar'});...
          'plot_pause_on_save', false, @islogical;...
        },...
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
        %NOTICE: out.mod.* have one less element than out.source.*, so keep that in mind and don't mix them!
        for i=1:numel(out.mod.dat);
          %build filename
          filename=strrep(out.file_root{1},'.png',...
            ['.',v.plot_temp_stat_list{s},'-triang',out.file_smooth,'.',str.rep(out.mod.names{i},' ','_'),'.png']...
          );
          %plot only if not done yet
          if exist(filename,'file')==0
            plotting.figure(v.varargin{:});
            %compute the correlation coefficient between this model and mod_ref
            d=out.mod.ref.scale(v.plot_temp_stat_func{s},'functional').interp(out.mod.dat{i}.t).stats2(...
              out.mod.dat{i}.scale(v.plot_temp_stat_func{s},'functional'),...
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
              'plot_title',[out.mod.names{i},' ',v.plot_temp_stat_title{s},' ',title_suffix]...
            );
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
        %plot degree-mean corrcoeff
        filename=strrep(out.file_root{1},'.png',...
          ['.',v.plot_temp_stat_list{s},'-dmean',out.file_smooth,'.',strjoin(str.rep(out.source.names,' ','_'),'-'),'.png']...
        );
        %plot only if not done yet
        if exist(filename,'file')==0
          %init plot counter and data container
          d=zeros(numel(out.mod.dat),out.mod.ref.lmax+1);
          %loop over all sources
          for i=1:numel(out.mod.dat);
            %compute it
            d(i,:)=out.mod.ref.scale(v.plot_temp_stat_func{s},'functional').interp(out.mod.dat{i}.t).stats2(...
                   out.mod.dat{i}.scale(v.plot_temp_stat_func{s},'functional'),...
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
          switch v.plot_type
          case 'bar'
            bar(v.plot_min_degree:out.mod.ref.lmax,d(:,good_idx)','EdgeColor','none');
            colormap jet
          case 'line'
            plot(v.plot_min_degree:out.mod.ref.lmax,d(:,good_idx)','Marker','o');
          end
          %need to adapt y label for correlation coefficients
          switch v.plot_temp_stat_list{s}
          case 'corrcoeff'
            y_label_str='Corr. Coeff. []';
          otherwise
            y_label_str=gravity.functional_label(v.plot_temp_stat_func{s});
          end
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_legend',out.mod.names,...
            'plot_ylabel',y_label_str,...
            'plot_xlabel','SH degree',...
            'plot_xlimits',[v.plot_min_degree-1,v.plot_max_degree+1],...
            'plot_title',['degree-mean ',v.plot_temp_stat_title{s},' ',title_suffix]...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=plot_stats(obj,product,varargin)
      obj=gswarm.plot_spatial_stats( obj,product,varargin{:});
      obj=gswarm.plot_temporal_stats(obj,product,varargin{:});
%       obj=gswarm.plot_low_degrees(obj,product,varargin{:});
    end
    function obj=plot_parametric_decomposition(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          't_mod_f',                                        1, @(i) isnumeric(i) && isscalar(i);...
          'polyorder',                                      2, @(i) isnumeric(i) && isscalar(i);...
          'plot_functional',                           'eqwh', @gravity.isfunctional;...
          'polynames',      {'constant','linear','quadratic'}, @iscellstr;...
          'plot_poly_range', [      inf,     inf,        inf], @isnumeric;...
          'sin_period',              [      12,            6], @isnumeric;...
          'sin_names',               {'yearly','semi-yearly'}, @iscellstr;...
          'plot_amplitude_range'     [inf,               inf], @isnumeric;...
          'sin_period_unit',                        'months' , @ischar;...
          'plot_pause_on_save',                         false, @islogical;...
        },...
        rmfield(product.metadata,'sources')...
      },varargin{:});
    
% %tmp code
% data_file='tmp_data.mat';
% if exist(data_file,'file')    
%   load(data_file)
% else
% %tmp code

      %collect the models 
      out=gswarm.plot_ops(obj,product,v.varargin{:});
      %loop over all models
      out.pd=cell(size(out.source.signal));
      for i=1:out.source.n
        %compute parametric decompositions
        out.pd{i}=out.source.signal{i}.parametric_decomposition(...
          't_mod_f',v.t_mod_f,...
          'polynomial',0:v.polyorder,...
          'sinusoidal',time.num2duration(v.sin_period,v.sin_period_unit),...
          't0',seconds(obj.start-out.source.signal{i}.epoch)...
        );
        %assemble title suffix
        title_suffix=out.title_startstop;
        if ~isempty(out.title_smooth);title_suffix=[title_suffix,newline,out.title_smooth];end %#ok<AGROW>
        for j=0:v.polyorder
          par_name=['p',num2str(j)];
          %build filename
          filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.',par_name,out.file_smooth,'.',strrep(out.source.names{i},' ','-'),'.png']);
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
              'plot_title',[v.polynames{j+1},' term for ',out.source.names{i},' ',title_suffix]...
            );
            if isfinite(v.plot_poly_range(j+1))
              caxis([-v.plot_poly_range(j+1),v.plot_poly_range(j+1)])
            end
            colormap jet
            cb.zero_white;
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
        for j=1:numel(v.sin_period)
          %build filename for amplitude
          filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.a',num2str(j),out.file_smooth,'.',strrep(out.source.names{i},' ','-'),'.png']);
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
              'plot_title',[v.sin_names{j},' amplitude for ',out.source.names{i},' ',title_suffix]...
            );
            if isfinite(v.plot_amplitude_range(j))
              caxis([0,v.plot_amplitude_range(j)])
            end
            colormap jet
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
          %build filename for phase
          filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.f',num2str(j),out.file_smooth,'.',strrep(out.source.names{i},' ','-'),'.png']);
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
              'plot_title',[v.sin_names{j},' phase for ',out.mod.names{i},' ',title_suffix]...
            );
            caxis([-v.sin_period(j),v.sin_period(j)]/2)
            colormap jet
            cb.zero_white;
            if v.plot_pause_on_save; keyboard; end
            saveas(gcf,filename)
            str.say('Plotted',filename)
          end
        end
      end
%       
% %tmp code
% save(data_file,'out')
% end
% %tmp code
  
      for j=1:size(simplegrid.catchment_list,1)
        %build filename
        filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.',...
          strrep(simplegrid.catchment_list{j,1},' ','-'),out.file_smooth,...
          '.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          plotting.figure(v.varargin{:});
          legend_str=cell(1,out.source.n);
          for i=1:out.source.n
            out.catch{i,j}=out.mod.res{i}.scale(v.plot_functional,'functional').grid('Nlat',100,'Nlon',200).catchment(...
              simplegrid.catchment_list{j,1},...
              'parametric_decomposition',false,...
              v.varargin{:}...
            );
            legend_str{i}=out.source.names{i};
          end
          plotting.enforce(...
            v.varargin{:},...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_line_style',cells.deal({'o-'},size(out.source.signal)),...
            'plot_legend',legend_str...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      
      for i=1:out.source.n
        %build filename
        filename=strrep(out.file_root{1},'.png',['.',v.plot_functional,'.',...
          'std',out.file_smooth,'.',...
          strrep(out.source.names{i},' ','-'),'.png']);
        %plot only if not done yet
        if exist(filename,'file')==0
          plotting.figure(v.varargin{:});
          out.mod.res{i}.scale(v.plot_functional,'functional').grid('Nlat',100,'Nlon',200).stats('mode','std').imagesc(...
            'boxes',simplegrid.catchment_list...
          );
          plotting.enforce(...
            v.varargin{:},...
            'plot_line_color','',...
            'plot_title',['temporal STD of ',out.mod.names{i}]...
          );
          if v.plot_pause_on_save; keyboard; end
          saveas(gcf,filename)
          str.say('Plotted',filename)
        end
      end
      
    end
    %% utils
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

      datafilename=file.unresolve('~/data/gswarm/analyses/2018-10-04/d.mat');
      p=cellfun(@(i) dataproduct(i,'plot_dir',fileparts(datafilename)),{...   
        'gswarm.swarm.all.TN-03_2.catchments',...
        'gswarm.swarm.all.TN-03_2.unsmoothed',...
        'gswarm.swarm.all.TN-03_2.smoothed'...
      },'UniformOutput',false);

    %NOTICE: There's problem with the handling of model errors!!!

      %save version numbers into latex table
      file.strsave(str.rep(datafilename,'d.mat','versions.tex'),gswarm.sourcelatextable(p{1}));
      
      %load data if already available
      if exist(datafilename,'file')~=0 && ~recompute
        load(datafilename)
      else
        d=datastorage('debug',true);
        for i=1:p{1}.nr_sources
          d=d.init(p{1}.sources(i),'recompute',recompute);
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
    %% legacy
    %     %better use plot_spatial_stats, more general
%     function obj=plot_rms_ts(obj,product,varargin)
%       obj.log('@','in','product',product)
%       plot_functional =product.mdget('plot_functional');
%       derived_quantity=product.mdget('plot_derived_quantity');
%       p=inputParser;
%       p.KeepUnmatched=true;
%       %parse optional parameters as defined in the metadata
%       p=product.plot_args(p,varargin{:});
%       %retrieve plot elements (repetitive processing of parameters)
%       e=obj.plot_elements(p,product);
%       %loop over all data
%       for t=1:numel(e.startlist)
%         obj.log('@','iter','start',e.startlist(t),'stop',e.stoplist(t))
%         %plot filename arguments
%         filename_args=[product.file_args('plot'),{...
%           'start',e.startlist(t),...
%           'stop', e.stoplist(t),...
%           'timestamp',true,...
%           'remove_part','',...
%           'prefix',p.Results.plot_file_prefix,...
%           'suffix',e.suffix...
%         }];
%         %plot filename
%         filename=product.dataname.file(filename_args{:});
%         if isempty(dir(filename))
%           %get the data for the current segment
%           obj_curr=obj.trim('start',e.startlist(t),'stop',e.stoplist(t),'dn_list',e.sources);
%           %make sure there is data 
%           if any(cell2mat(obj_curr.vector_method_tr('all','nr_valid'))>1)
%             %need all data to be in the same time domain
%             obj_curr=obj_curr.interp('all');
%             %build data array
%             bar_y=zeros(max(obj_curr.length('all')),numel(e.sources));
%             bar_t=zeros(max(obj_curr.length('all')),numel(e.sources));
%             for di=1:numel(e.sources)
%               tmp=obj_curr.data_get_scalar(e.sources{di}).scale(plot_functional,'functional').(derived_quantity);
%               bar_y(:,di)=tmp(:,end);
%               bar_t(:,di)=datenum(obj_curr.data_get_scalar(e.sources{1}).t);
%             end
%             %filter out nans
%             good_idx=~all(isnan(bar_y),2);
%             %plot it
%             figure('visible',p.Results.plot_visible);
%             h=bar(bar_t(good_idx,:),bar_y(good_idx,:),'EdgeColor','none');
%             %enforce plot preferences
%             product.enforce_plot
%             %build plot annotation structure
%             bh=cell(size(e.sources));
%             for di=1:numel(e.sources)
%               bh{di}=struct(...
%                 'mask',obj_curr.data_get_scalar(e.sources{di}).mask,...
%                 'y_mean',0,...
%                 'handle',h(di),...
%                 'title',e.sources{di}.str,...
%                 'xlabel','',...
%                 'ylabel',plot_functional,...
%                 'y_units',gravity.functional_units(plot_functional),...
%                 'legend',{{e.sources{di}.name}}...
%               );
%             end
%             %annotate plot
%             obj.plot_annotate(bh,product,datanames.transmute(e.sources),varargin{:})
%             %remove outline
%             set(h,'edgecolor','none')
%             %save this plot
%             if v.plot_pause_on_save; keyboard; end
%             saveas(gcf,filename)
%             obj.log('@','iter','plot saved',filename)
%             % user feedback
%             if strcmp(p.Results.plot_visible,'off')
%               disp(['gswarm.plot_rms_ts: plotted ',product.name,' to file ',filename])
%             end
%           else
%             disp(['gswarm.plot_rms_ts: not enough data to plot ',product.name,' to file ',filename,' (skipped)'])
%           end
%         else
%           obj.log('@','iter','skipped',product,'file already exists',filename)
%         end
%       end
%       obj.log('@','out')
%     end 
  end
end