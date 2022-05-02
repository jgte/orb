classdef gswarm
  properties(Constant)
    l1bdir_options={...
      '~/data/gswarm';...
      './data/gswarm';...
    };
    start_date=datetime('2013-12-01');
    % for i=1:12
    %   disp([num2str(i),' : ',datestr(...
    %     dateshift(dateshift(datetime(['2022-',num2str(i,'%02i'),'-15']),'start','months')-calmonths(2),'start','quarter')-days(1)...
    %   )]);
    % end
    %NOTICE: assigns the end date of the processing according to the current month as follows:
    % 1  : 30-Sep-2021
    % 2  : 30-Sep-2021
    % 3  : 31-Dec-2021
    % 4  : 31-Dec-2021
    % 5  : 31-Dec-2021
    % 6  : 31-Mar-2022
    % 7  : 31-Mar-2022
    % 8  : 31-Mar-2022
    % 9  : 30-Jun-2022
    % 10 : 30-Jun-2022
    % 11 : 30-Jun-2022
    % 12 : 30-Sep-2022
     stop_date=dateshift(dateshift(datetime('now'),'start','months')-calmonths(2),'start','quarter')-days(1);
  end
  methods(Static)
    function out=dir(type)
      switch type
        case 'data'
          for i=1:numel(gswarm.l1bdir_options)
            if file.exist(gswarm.l1bdir_options{i})
              out=gswarm.l1bdir_options{i};
              return
            end
          end
        %add more directories here
      end
    end
    function obj=load_models(obj,product,varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{...
        {...
          'import_dir',        '.',@(i) ischar(i) && exist(i, 'dir')~=0;...
          'model_format','unknonw',@ischar;...
          'date_parser', 'unknonw',@ischar;...
          'consistent_GM',   false,@islogical;...
          'consistent_R',    false,@islogical;...
          'max_degree'           0,@num.isscalar;...
          'use_GRACE_C20',  'none',@ischar;...
          'delete_C00',      false,@islogical;...
          'delete_C20',      false,@islogical;...
          'model_span', {time.zero_date,time.inf_date},@(i) isnumeric(i) && numel(i)==2;...
          'static_model',   'none',@ischar;...
          'overwrite_common_t',false,@islogical;...
        },...
        product.args...
      },varargin{:});
      %load all available data
        [s,e]=gravity.load_dir(...
        'datadir',v.import_dir,...
        'format',v.model_format,...
        'descriptor',product.name,...
        v.varargin{:}...
      );
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
      derived_quantity=product.mdget('plot_spatial_diff_quantity');
      functional      =product.mdget('functional');
      stats           =product.mdget('stats');
      stats_outlier_iter=product.mdget('stats_outlier_iter');
      stats_detrend   =product.mdget('stats_detrend');
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
        'outlier_iter',stats_outlier_iter,...
        'detrend',stats_detrend,...
        'period', stats_period,...
        'overlap',stats_overlap...
      );
      %add field with degree
      s.degree=0:g.lmax;
      %save data
      obj=obj.data_set(product,s);
    end
    %% plots
    function out=file_root(obj,product)
      %NOTICE: 'storage_period' 'direct' is needed so that only one name is returned, otherwise
      %         a cell array with numerous filenames may be returned, depending on the storage period
      out=product.mdset('storage_period','direct').file('plot',...
        'start',obj.start,...
        'stop',obj.stop,...
        'start_timestamp_only',false,...
        'no_extension',true,...
        'sub_dirs','single',...
        'use_storage_period',true...
      );
    end
    function pod=plot_ops(obj,product,varargin)
      %TODO: make inventory of the the operations done in this method, maybe wrap them in a switch loop, as is the case with load_models_op
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_min_degree',        2 ,@num.isscalar;... %NOTICE: this needs to be handled externally!
          'plot_max_degree',      inf ,@num.isscalar;...
          'plot_smoothing_degree', [] ,@num.isscalar;...
          'plot_smoothing_method', '' ,@ischar;...
          'plot_spatial_mask', 'none' ,@iscellstr;...
          'stats_relative_to', 'none' ,@ischar;...
          'model_types',        {'*'} ,@iscellstr;...
          'plot_lines_over_gaps_narrower_than',days(120),@isduration;...
          'plot_time_domain_source',0 ,@num.isscalar;...
          'inclusive',          false ,@islogical;... %generally we want plots to be exclusive, i.e. not show non-common data (e.g. C20 from 2002-04)
        },...
        product.args({'stats_relative_to','model_types'}),...
        product.plot_args...
      },varargin{:});
      % update start/stop considering the requested inclusive
      obj=obj.startstop_update('start',gswarm.production_date('start'),'stop',gswarm.production_date('stop'),'inclusive',v.inclusive);
      obj.log('@','startstop_update','product',product,'start',obj.start,'stop',obj.stop)
      %retrieve load/saving of plotdata: handles plot_force and plot_save_data
      [loaddata,savedata]=plotting.forcing(v.varargin{:});
      %first build data filename;
      %TODO: this datafilename is messed up, the smoothing only appears as 'gauss'
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
      if v.plot_time_domain_source==0
        %set time domain common to all sources
        pod.t=simpletimeseries.t_mergev(pod.source.dat);
      else
        %set time domain equal to the requested source
        pod.t=pod.source.dat{v.plot_time_domain_source}.t;
      end
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
        for i=1:numel(pod.source.dat)
          str.say('Applying',v.plot_spatial_mask,'mask to product',pod.source.datanames{i}.name)
          pod.source.dat{i}=pod.source.dat{i}.spatial_mask(v.plot_spatial_mask);
        end
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
      %filename particles
      pod.file_root=gswarm.file_root(obj,product);
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
        out=pod;
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
          'plot_max_lines',inf,@num.isscalar;...
          'show_legend_stats',  'yes' ,@str.islogical;...
          'plot_lines_over_gaps_narrower_than',days(120),@isduration;...
          'plot_legend_sorting','ascend', @ischar;...
          'plot_force',false,@islogical;...
          'plot_legend_include_smoothing',false,@islogical;...
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
        filename=file.build(v.pod.file_root,['C',num2str(d),',',num2str(o)],v.pod.file_smooth,'png');
        %plot only if not done yet
        if ~file.exist(filename,v.plot_force)
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
            %add statistics to the legend (unless this is the product from which the stats are derived)
            if k~=v.pod.source.ref_idx
              mod_ref_now=v.pod.mod.ref.ts_C(d,o).interp(...
                ts_now{k}.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than...
              );
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
          stats=stats(trivial_idx);
          %sort it (if requested)
          switch v.plot_legend_sorting
          case {'ascend','descend'};   [~,idx]=sort(cellfun(@(i) i.rms,stats),'ascend');
          case {'none','no'};             idx=1:numel(stats);
          otherwise
            error(['Cannot handle ''plot_legend_sorting'' with value ''',v.plot_legend_sorting,'''.'])
          end
          %truncate it
          idx=idx(1:min(numel(ts_now),v.plot_max_lines));
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
          %enforce it
          product.enforce_plot(v.varargin{:},...
            'plot_ylabel','[ ]',...
            'plot_legend',legend_str,...
            'plot_line_color_order',idx,...
            'plot_title',v.plot_title,...
            'plot_title_default',['C',num2str(d),',',num2str(o),' ',title_suffix]...
          );
          plotting.save(filename,v.varargin{:})
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
          'plot_min_degree',           2, @num.isscalar;...
          'plot_max_degree',         inf, @num.isscalar;...
          'plot_functional',     'geoid', @gravity.isfunctional;...
          'plot_type',            'line', @(i) cells.isincluded(i,{'line','bar'});...
          'plot_show_legend_stats',false, @islogical;...
          'plot_max_nr_lines',        20, @num.isscalar;...
          'plot_spatial_stat_list'       ,{'diff','monthly'}   , @iscellstr; ... TODO: corr
          'plot_spatial_diff_quantity'   ,{'cumdas','gridmean'}, @(i) all(cellfun(@(j) ismethod(gravity.unit(1),j),i));... %relevant if plot_spatial_stat_list has 'diff'
          'plot_spatial_monthly_quantity',{'das','triang'}     , @iscellstr;...   %relevant if plot_spatial_stat_list has 'monthly': one plot per month
          'plot_spatial_monthly_error'   ,false                , @islogical;...   %relevant if plot_spatial_stat_list has 'monthly': include format error in das monthly plots
          'plot_spatial_monthly_last'    ,inf                  , @num.isscalar;...%relevant if plot_spatial_stat_list has 'monthly': only plots these last months
          'plot_spatial_monthly_triang_caxis',[-inf inf]       , @(i) isnumeric(i) && numel(i)==2;...%relevant if plot_spatial_stat_list has 'monthly': define triang plots' caxis
          'plot_legend_include_smoothing',          false, @islogical;...
          'plot_lines_over_gaps_narrower_than', days(120), @isduration;...
          'plot_signal',           false, @islogical;...
          'plot_summary',           true, @islogical;...
          'plot_logy',             false, @islogical;...
          'plot_legend_sorting','ascend', @ischar;...
          'plot_force',            false, @islogical;...
          'plot_detrended',           '', @ischar;...   %see simpledata.detrend for modes, empty means no detrending
          'plot_outlier_iter',     false, @isfinite;... %number of outlier removal iters, see simpledata.stats

        },...
        product.plot_args...
      },varargin{:});
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
        for qi=1:numel(v.plot_spatial_diff_quantity)

          %do not detrend for gridmean
          switch v.plot_spatial_diff_quantity{qi}
          case 'gridmean'
            plot_detrended='none';
          otherwise
            plot_detrended=v.plot_detrended;
          end

          %% prepare data for (summary) diff rms plots

          %build filenames
          filenames={...
            file.build(v.pod.file_root, v.plot_spatial_diff_quantity{qi},            v.pod.file_smooth,'png');...
            file.build(v.pod.file_root,[v.plot_spatial_diff_quantity{qi},'_summary'],v.pod.file_smooth,'png');...
          };
          %prepare plot data only if needed
          if any(~file.exist(filenames,v.plot_force))
            %build data array (number of rows is short by one if v.plot_signal && stats_relative_to)
            y=cell(1,numel(v.pod.mod.res));t=cell(size(y));yc=zeros(size(y));
            for di=1:numel(v.pod.mod.res)+1
              %branch to retrieve data that the stats are derived relative to
              switch di
              case numel(v.pod.mod.res)+1
                if v.plot_signal && ~str.none(product.mdget('stats_relative_to','default','none'))
                  tmp=v.pod.source.dat{v.pod.source.ref_idx};
                else
                  continue
                end
              otherwise
                tmp=v.pod.mod.res{di};
              end
              %operate
              tmp=tmp.interp(...
                  v.pod.t,...                            %interp to common time domain
                  'interp_over_gaps_narrower_than',...   %and remove lines connecting over large gaps
                  v.plot_lines_over_gaps_narrower_than...
                ).addgaps(...                     %add explicit gaps: catches the case when v.pod.t is implicitly gapped
                  v.plot_lines_over_gaps_narrower_than...
                ).scale(...
                  v.plot_functional,'functional'...      %scale to functional
                ).outlier(...                            %maybe detrend and remove outliers, depends on value of :
                  'outlier_iter',v.plot_outlier_iter,... % - handled inside the obj.outlier method
                       'detrend',plot_detrended...       % - handled inside the obj.detrend method (called by obj.outlier)
                );
              %save time
              t{di}=tmp.t;
              tmp=tmp.(...
                v.plot_spatial_diff_quantity{qi}...         %compute derived quantity
              );
              %propagate relevant data point
              y{di}=tmp(:,end);
              %average it
              yc(di)=sum(y{di},1,'omitnan')./sum(~isnan(y{di}));
            end
            %sort it (if requested)
            switch v.plot_legend_sorting
            case {'ascend','descend'};   [yc_sorted,idx]=sort(yc,'ascend');
            case {'none','no'};           yc_sorted=yc;idx=1:numel(yc);
            otherwise
              error(['Cannot handle ''plot_legend_sorting'' with value ''',v.plot_legend_sorting,'''.'])
            end
            y_sorted=y(idx);
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
                y_sorted=  y_sorted(1:v.plot_max_nr_lines);
               yc_sorted= yc_sorted(1:v.plot_max_nr_lines);
              legend_str=legend_str(1:v.plot_max_nr_lines);
                     idx=idx(       1:v.plot_max_nr_lines);
            end
            %trim nans for plotting
            t_sorted=cell(size(y_sorted));
            for di=1:numel(y_sorted)
              plot_idx=any(num.trim_NaN(y_sorted{di}),2);
              y_sorted{di}=y_sorted{di}(plot_idx);
              t_sorted{di}=t{di}(plot_idx);
            end
            %compute abs value in case of log scale
            if str.logical(v.plot_logy)
              y_sorted=cellfun(@abs,y_sorted,'UniformOutput',false);
            end
            %plot abs value for gridmean
            switch v.plot_spatial_diff_quantity{qi}
            case 'gridmean'
              y_sorted=cellfun(@abs,y_sorted,'UniformOutput',false);
            otherwise
              %do nothing
            end
          end


          %% plot diff cumdas/gridmean

          fn_idx=1;
          if ~file.exist(filenames{fn_idx},v.plot_force)
            %plot it
            plotting.figure(v.varargin{:}); hold on
            for di=1:numel(y_sorted)
              switch v.plot_type
              case 'bar';   bar(t_sorted{di},y_sorted{di},'EdgeColor','none');
              case 'line'; plot(t_sorted{di},y_sorted{di},'Marker','o');
              end
            end
            %deal with legend stats
            if v.plot_show_legend_stats
              for di=1:numel(legend_str)
                legend_str{di}=[legend_str{di},...
                     ' ',num2str(mean(y_sorted{di}(~isnan(y_sorted{di}))),2),...
                 ' +/- ',num2str( std(y_sorted{di}(~isnan(y_sorted{di}))),2),...
              ' \Sigma=',num2str( sum(y_sorted{di}(~isnan(y_sorted{di}))),2)];
              end
            end
            %deal with title prefix
            switch v.plot_spatial_diff_quantity{qi}
            case 'gridmean'; title_prefix='Mean diff.';
            case 'cumdas';   title_prefix='RMS diff.';
            otherwise
              error(['cannot handle plot_spatial_diff_quantity with value ''',...
                v.plot_spatial_diff_quantity{qi},'''.'])
            end
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_legend',legend_str,...
              'plot_ylabel',gravity.functional_label(v.plot_functional),...
              'plot_xdate',true,...
              'plot_xlimits',[min(cellfun(@(i) i(1),t_sorted)),max(cellfun(@(i) i(end),t_sorted))+days(1)],...
              'plot_line_color_order',idx,...
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({title_prefix,v.pod.title_wrt,title_suffix})...
            );
            plotting.save(filenames{fn_idx},v.varargin{:})
          else
            disp(['NOTICE: plot already available: ',filenames{fn_idx}])
          end

          %% plot cumulative diff rms (summary)

          fn_idx=2;
          if ~file.exist(filenames{fn_idx},v.plot_force) && numel(yc_sorted)>1 && v.plot_summary
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
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({'Cum. residual',v.pod.title_wrt,v.pod.title_startstop,title_suffix})...
            );
            set(gca,...
              'YTick',1:numel(legend_str),...
              'yticklabels',plotting.legend_replace_clean(v.varargin{:},'plot_legend',flipud(legend_str(:)))...
            );
            plotting.save(filenames{fn_idx},v.varargin{:})
          else
            disp(['NOTICE: plot already available: ',filenames{fn_idx}])
          end

        end
      end

      %% epoch-wise plots

      %check if this plot is requested
      if cells.isincluded(v.plot_spatial_stat_list,'monthly')
        title_suffix=strjoin({v.pod.title_masking,v.pod.title_smooth},' ');
        if isfinite(v.plot_spatial_monthly_last)
          assert(v.plot_spatial_monthly_last>0,'the parameter ''plot_spatial_monthly_last'' must be positive')
          %NOTICE: there's no need to decrement start_idx because the last model is a empty by definition
          start_idx=numel(v.pod.t)-v.plot_spatial_monthly_last;
        else
          start_idx=1;
        end
        for i=start_idx:numel(v.pod.t)
          %gather models valid now
          dat    =cellfun(@(j) j.interp(v.pod.t(i)),v.pod.source.dat,'UniformOutput',false);
          dat_idx=cellfun(@(j) j.nr_valid>0,dat);
          %gather error if requested
          if v.plot_spatial_monthly_error
            try
              dat_error=cellfun(@(j) j.interp(v.pod.t(i)),v.pod.source.error,'UniformOutput',false);
            catch
              v.plot_spatial_monthly_error=false;
              disp('WARNING: cannot plot monthly error because it is not available. Implementation needed.')
            end
          end
          %check if there's any data to plot
          if all(~dat_idx); continue;end
          %reduce data
          dat=dat(dat_idx);
          if v.plot_spatial_monthly_error
            dat_error=dat_error(dat_idx);
          end
          legend_str=upper(v.pod.source.names(dat_idx));
          %loop over all requested derived quantities
          for qi=1:numel(v.plot_spatial_monthly_quantity)
            %build filename
            filename=file.build(...
              v.pod.file_root,v.plot_spatial_monthly_quantity{qi},...
               datestr(v.pod.t(i),'YYYYmmdd'),v.pod.file_smooth,'png');...
            if ~file.exist(filename,v.plot_force)
              %branch on the type of plot
              switch v.plot_spatial_monthly_quantity{qi}
              case 'das'
                plotting.figure(v.varargin{:});
                for j=1:numel(dat)
                  dat{j}.plot('mode',v.plot_spatial_monthly_quantity{qi},'functional',v.plot_functional);
                end
                %enforce it
                product.enforce_plot(v.varargin{:},...
                  'plot_legend',legend_str,...
                  'plot_ylabel',gravity.functional_label(v.plot_functional),...
                  'plot_colormap','jet',...
                  'plot_title',v.plot_title,...
                  'plot_title_default',str.show({datestr(v.pod.t(i),'yyyy-mm'),'degree-RMS',title_suffix})...
                );
                if v.plot_spatial_monthly_error
                  %get previous lines
                  lines_before=findobj(gca,'Type','line');
                  %plot errors
                  for j=1:numel(dat)
                    dat_error{j}.plot('method',v.plot_spatial_monthly_quantity{qi},'functional',v.plot_functional,'line','--');
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
              case {'triang','trianglog10'}
                [~,l,w]=plotting.subplot(0,numel(dat));
                plotting.figure(v.varargin{:},...
                  'plot_size',200+[0,0,21,9]*30.*[1 1 w l]);
                for j=1:numel(dat)
                  plotting.subplot(j,numel(dat),true);
                  dat{j}.plot('method',v.plot_spatial_monthly_quantity{qi},...
                    'colormap',str.clean(v.plot_colormap,{'opt','zero'}),... %TODO: fix this once plotting.colormap is implemented and used in gravity.plot
                    'functional',v.plot_functional...
                  );
                  %enforce it
                  product.enforce_plot(v.varargin{:},...
                    'plot_colormap','none',...
                    'plot_caxis',v.plot_spatial_monthly_triang_caxis,...
                    'plot_title',v.plot_title,...
                    'plot_title_default',str.show({legend_str{j},datestr(v.pod.t(i),'yyyy-mm'),title_suffix})...
                  );
                end
              otherwise
                error(['Cannot handle plot_spatial_monthly_quantity as ''',v.plot_spatial_monthly_quantity{qi},'''.'])
              end
              plotting.save(filename,v.varargin{:})
            else
              disp(['NOTICE: plot already available: ',filename])
            end
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
          'plot_min_degree',     2,          @num.isscalar;...
          'plot_max_degree',     inf,        @num.isscalar;...
          'plot_rms_caxis',      [-inf inf], @isnumeric;...
          'plot_functional',     'geoid',    @gravity.isfunctional;...
          'plot_type',           'line',     @(i) cells.isincluded(i,{'line','bar'});...
          'plot_max_nr_lines',   20,         @num.isscalar;...
          'plot_temp_stat_list', {...
                                    'corrcoeff',...
                                    'rms',...
                                    'std'...
                                  },         @iscellstr; ...
          'plot_temp_stat_title',{...
                                    'temporal corr. coeff.',...
                                    'temporal RMS\Delta',...
                                    'temporal STD\Delta'....
                                 },          @iscellstr; ...
          'plot_legend_include_smoothing', false, @islogical;...
          'plot_legend_sorting',        'ascend', @ischar;...
          'plot_corrcoeff_caxis',         [-1,1], @(i) isnumeric(i) && numel(i)==2 && all(abs(i)<=1)
          'plot_summary',                   true, @islogical;...
          'plot_force',                    false, @islogical;...
          'plot_detrended',       '', @ischar;...   %see simpledata.detrend for modes, empty means no detrending (defined in simpledata.stats)
          'plot_outlier_iter',     0, @isfinite;... %number of outlier removal iters, see simpledata.stats
          'plot_lines_over_gaps_narrower_than',     days(120), @isduration;...
        },...
        product.plot_args...
      },varargin{:});
      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end

      %loop over all statistics
      for s=1:numel(v.plot_temp_stat_list)
        %maybe nothing is requested to be plotted
        if strcmp(v.plot_temp_stat_list{s},'none'); continue; end
        %do not detrend for correlation coefficients
        switch v.plot_temp_stat_list{s}
        case {'corrcoeff','rms'}
          plot_detrended='none';
        otherwise
          plot_detrended=v.plot_detrended;
        end
        %collect arguments for call to stats/stats2
        stats_args={...
          'mode','obj',...
          'struct_fields',v.plot_temp_stat_list(s),...
          'outlier_iter',v.plot_outlier_iter,...
          'detrend',plot_detrended,...
          'period',seconds(inf)...
        };

        %% triangular plots

        %loop over all sources
        %NOTICE: out.mod.* have one less element than out.source.*, so keep that in mind and don't mix them!
        for i=1:numel(v.pod.mod.dat)
          %build filename
          filename=file.build(v.pod.file_root,...
            {v.plot_temp_stat_list{s},'triang'},v.pod.source.file_smooth{i},strsplit(v.pod.mod.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename,v.plot_force)
            plotting.figure(v.varargin{:});
            tmp=v.pod.mod.dat{i}.scale(v.plot_functional,'functional');
            if isfield(v.pod.mod,'ref')
              %compute the requested stat between this model and mod_ref
              d=tmp.stats2(...
                v.pod.mod.ref.scale(v.plot_functional,'functional').interp(...
                tmp.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than...
              ),...
                stats_args{:}...
              );
            else
              d=tmp.stats(...
                stats_args{:}...
              );
            end
            %set y-units
            d.units(:)={gravity.functional_units(v.plot_functional)};
            %plot it
            d.plot('method','triang');
            %enforce it
            product.enforce_plot(v.varargin{:},...
              'plot_caxis',cells.c2m(v.(['plot_',v.plot_temp_stat_list{s},'_caxis'])),...
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({...
                upper(v.pod.mod.names{i}),v.plot_temp_stat_title{s},v.pod.title_wrt,...
                v.pod.title_masking,v.pod.source.title_startstop{v.pod.source.mod_idx(i)},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{v.pod.source.mod_idx(i)}}),newline)...
              })...
            );cb.nan;
            %need to adapt caxis label for correlation coefficients
            switch v.plot_temp_stat_list{s}
            case 'corrcoeff'
              cb.label('Corr. Coeff. []');
            end
            plotting.save(filename,v.varargin{:})
          end
        end

        %% prepare data for (cumulative) degree-mean/mean corr coeff plots

        %plot degree-mean stat
        filenames={...
          file.build(v.pod.file_root,{v.plot_temp_stat_list{s},'dmean'        },v.pod.file_smooth,'png');...
          file.build(v.pod.file_root,{v.plot_temp_stat_list{s},'dmean_summary'},v.pod.file_smooth,'png');...
        };
        %prepare plot data only if needed
        if any(~file.exist(filenames,v.plot_force))
          %init plot counter and data container
          d=zeros(numel(v.pod.mod.dat),v.pod.mod.dat{1}.lmax+1);
          %loop over all sources
          for i=1:numel(v.pod.mod.dat)
            if isfield(v.pod.mod,'ref')
              %compute it
              d(i,:)=v.pod.mod.ref.scale(v.plot_functional,'functional').interp(...
                v.pod.mod.dat{i}.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than...
              ).stats2(...
                v.pod.mod.dat{i}.scale(v.plot_functional,'functional'),...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'detrend',plot_detrended,...
                'outlier_iter',v.plot_outlier_iter,...
                'period',seconds(inf)...
              ).dmean;
            else
              d(i,:)=v.pod.mod.dat{i}.scale(v.plot_functional,'functional').stats(...
                'mode','obj',...
                'struct_fields',v.plot_temp_stat_list(s),...
                'detrend',plot_detrended,...
                'outlier_iter',v.plot_outlier_iter,...
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
        if ~file.exist(filenames{fn_idx},v.plot_force)
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
            'plot_logy',plot_logy,...
            'plot_title',v.plot_title,...
            'plot_title_default',str.show({'degree mean',v.plot_temp_stat_title{s},newline,title_suffix})...
          );
          plotting.save(filenames{fn_idx},v.varargin{:})
        else
          disp(['NOTICE: plot already available: ',filenames{fn_idx}])
        end

        %% plot cumulative degree-mean/mean corr coeff stat  (only if not done yet)

        fn_idx=2;
        if ~file.exist(filenames{fn_idx},v.plot_force) && numel(dc_sorted)>1 && v.plot_summary
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
            'plot_title',v.plot_title,...
            'plot_title_default',str.show({...
              'Degrees',v.plot_min_degree,'-',v.plot_max_degree,'cum. degree mean',v.plot_temp_stat_title{s},newline,...
              title_suffix...
            })...
          );
          set(gca,...
            'YTick',1:numel(legend_str),...
            'yticklabels',plotting.legend_replace_clean(v.varargin{:},'plot_legend',flipud(legend_str(:)))...
          );
          plotting.save(filenames{fn_idx},v.varargin{:})
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
          'plot_spatial_step',                               1, @num.isscalar;...
          'catchment_list',     simplegrid.catchment_list(:,1), @iscellstr; ....
          'plot_std_range'                                 inf, @num.isscalar;...
          'plot_std_colormap'                     'parulazero', @ischar;...
          'plot_force',                                  false, @islogical;...
          'plot_detrended',       '', @ischar;...   %see simpledata.detrend for modes, empty means no detrending (defined in simpledata.stats)
          'plot_outlier_iter', false, @isfinite;... %number of outlier removal iters, see simpledata.stats
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
        if ~file.exist(filename,v.plot_force)
          plotting.figure(v.varargin{:});
          v.pod.source.dat{i...
            }.scale(v.plot_functional,'functional'...
            ).grid('spatial_step',v.plot_spatial_step...
            ).stats(...
              'detrend',v.plot_detrended,...
              'outlier_iter',v.plot_outlier_iter,...
              'mode','std'...
            ).imagesc(...
              'boxes',v.catchment_list,...
              'cb_title',gravity.functional_label(v.plot_functional)...
          );
          plotting.enforce(v.varargin{:},...
            'plot_legend_location','none',...
            'plot_caxis',[0,str.num(v.plot_std_range)],...
            'plot_colormap',v.plot_std_colormap,... %the "parula" colormap goes well with the red boxes
            'plot_xlabel','none','plot_ylabel','none',...
            'plot_title',v.plot_title,...
            'plot_title_default',str.show({...
              'temporal STD of ',upper(v.pod.source.names{i}),...
              v.pod.title_masking,v.pod.source.title_startstop{i},...
              strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
            }));
          plotting.no_ticks
          plotting.save(filename,v.varargin{:})
        end
      end
    end
    function obj=plot_stats(obj,product,varargin)
      %collect the models
      pod=gswarm.plot_ops(obj,product,varargin{:});
      %call the lower tier plotting routines, re-use pod
      obj=gswarm.plot_spatial_stats( obj,product,varargin{:},'pod',pod);
      obj=gswarm.plot_temporal_stats(obj,product,varargin{:},'pod',pod);
      obj=gswarm.plot_low_degrees(   obj,product,varargin{:},'pod',pod);
    end
    function obj=plot_parametric_decomposition(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          't_mod_f',                                        1, @num.isscalar;... %TODO: implement this
          'polyorder',                                      2, @num.isscalar;...
          'plot_functional',                           'eqwh', @gravity.isfunctional;...
          'polynames',      {'constant','linear','quadratic'}, @iscellstr;...
          'sin_period',              [      12,            6], @isnumeric;...
          'sin_names',               {'yearly','semi-yearly'}, @iscellstr;...
          'sin_period_unit',                        'months' , @ischar;...
          'timescale',                               'years' , @ischar;...
          't0',                        datetime('2000-01-01'), @isdatetime;...
          'plot_spatial_step',                              1, @num.isscalar;...
          'plot_catchment_list',simplegrid.catchment_list(:,1),@iscellstr; ....
          'plot_catchment_boxes',                       false, @islogical; ....
          'plot_std_colormap',                       'parula', @ischar;...
          'plot_std_range'                                inf, @num.isscalar;...
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
      %operate
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
      %retrieve load/saving of plotdata: handles plot_force and plot_save_data
      loaddata=plotting.forcing(v.varargin{:});
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
          if ~file.exist(filename,v.plot_force)
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
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({...
                v.polynames{j+1},'term for',upper(v.pod.source.names{i}),...
                v.pod.title_masking,v.pod.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
            }));
            plotting.no_ticks
            plotting.save(filename,v.varargin{:})
          end
        end
        for j=1:numel(v.sin_period)
          %build filename for amplitude
          filename=file.build(v.pod.file_root,...
            v.plot_functional,{'a',j},v.pod.source.file_smooth{i},strsplit(v.pod.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename,v.plot_force)
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
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({...
                v.sin_names{j},'amplitude for',upper(v.pod.source.names{i}),...
                v.pod.title_masking,v.pod.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
              }));
            plotting.no_ticks
            plotting.save(filename,v.varargin{:})
          end
          %build filename for phase
          filename=file.build(v.pod.file_root,...
            v.plot_functional,{'f',j},v.pod.source.file_smooth{i},strsplit(v.pod.source.names{i},' '),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename,v.plot_force)
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
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({...
                v.sin_names{j},'phase for',upper(v.pod.source.names{i}),...
                v.pod.title_masking,v.pod.source.title_startstop{i},...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
              }));
            plotting.no_ticks
            plotting.save(filename,v.varargin{:})
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
          'plot_spatial_step',                              1, @num.isscalar;...
          'catchment_list',    simplegrid.catchment_list(:,1), @iscellstr; ....
          'parametric_decomposition',                    true, @islogical;... %do pardecomp on the catchment timeseries
          'polyorder',                                      2, @num.isscalar;...
          'polynames',          {'bias','linear','quadratic'}, @iscellstr;...
          'sin_period',              [      12,            6], @isnumeric;...
          'sin_period_unit',                        'months' , @ischar;...
          'sin_names',               {'yearly','semi-yearly'}, @iscellstr;...
          'timescale',                               'years' , @ischar;...
          't0',                        datetime('2000-01-01'), @isdatetime;...
          'plot_parametric_components',           {'p0','p1'}, @iscellstr;... %plot these components
          'plot_parametric_timestep_value',                 7, @num.isscalar;...
          'plot_parametric_timestep_units',            'days', @ischar;...
          'plot_legend_include_smoothing',              false, @islogical;...
          'plot_line_colormap',                      'spiral', @ischar;...
          'plot_lines_over_gaps_narrower_than',     days(120), @isduration;...
          'plot_parametric_t_masked',                   false, @isduration;...
        },...
        product.args...
      },varargin{:});

      %collect the models, unless given externally
      v=varargs.wrap('sources',{v,...
        {...
          'pod',[],@isstruct;...
        }...
      },varargin{:});
      %operate
      if isempty(v.pod); v.pod=gswarm.plot_ops(obj,product,v.varargin{:}); end
      %retrieve load/saving of plotdata: handles plot_force and plot_save_data
      loaddata=plotting.forcing(v.varargin{:});
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
            'suffix',strjoin(cells.rm_empty(...
              {v.pod.source.datanames{i}.name,v.pod.source.file_smooth{i},v.catchment_list{j},'catchment-data'}...
            ),'.')...
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
        if ~file.exist(filename,v.plot_force)
          plotting.figure(v.varargin{:});
          legend_str=cell(1,v.pod.source.n);
          %the default color scheme for lines only has a limited number of entries, so to
          %be safe against that, just build a colormap from scratch
          plot_line_colormap=plotting.line_color_map(v.plot_line_colormap,v.pod.source.n);
          %retrieve time domain to plot
          t_step=time.num2duration(v.plot_parametric_timestep_value,v.plot_parametric_timestep_units);
          if v.plot_parametric_t_masked
            t_domain=time.domain(...
              v.pod.catch{i,j}.pws.t_masked(1),...
              v.pod.catch{i,j}.pws.t_masked(end),...
              t_step...
            );
          else
            t_domain=time.domain(...
              v.pod.catch{i,j}.pws.time(1),...
              v.pod.catch{i,j}.pws.time(end),...
              t_step...
            );
          end
          %loop over all sources
          for i=1:v.pod.source.n
            v.pod.catch{i,j}=simplegrid.catchment_plot(v.pod.catch{i,j},...
              'time',t_domain,...
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
          %enforce it
          plotting.enforce(...
            v.varargin{:},...
            'plot_line_style','none',... %these 2 nones are needed so that the formating defined in simplegrid.catchment_plot is not over written
            'plot_line_color','none',...
            'plot_ylabel',gravity.functional_label(v.plot_functional),...
            'plot_xlabel','none',...
            'plot_legend',legend_str,...
            'plot_title',v.plot_title,...
            'plot_title_default',v.catchment_list{j}...
          );
          plotting.save(filename,v.varargin{:})
        end
      end

      %show stats header (if relevant)
      msg=cell(v.pod.source.n*(numel(v.catchment_list)+1),1);mc=0;
      tab_len=20;
      mc=mc+1;msg{mc}=str.tablify(tab_len,'catch','solution','bias [cm]','bias diff [cm]','trend [cm/y]','trend diff [cm/y]','corr coeff');
      norm_stats=zeros(v.pod.source.n,5,numel(v.catchment_list));
      for j=1:numel(v.catchment_list)
        %build filename
        %NOTICE: the name of the tex file cannot be the same as the name of the plots (latex doesn't like it for some reason)
        %        that's why there's the 'stats' particle below
        filename=file.build(v.pod.file_root,...
          v.plot_functional,'catch',v.pod.file_smooth,strsplit(v.catchment_list{j},' '),'stats','tex'...
        );
        %plot only if not done yet
        if ~file.exist(filename,v.plot_force)
          %show some stats, if parametric decomposition was made
          if isfield(v.pod.catch{1,j},'pws')
            latex_data=cell(v.pod.source.n,6);
            latex_name=cell(v.pod.source.n,1);
            latex_scale=[100 100 100 100 1]; %converts m to cm
            %compute bias of the reference (the p0 value is only relevant at t0 and it doesn't really give a good idea of the bias)
            v.pod.catch{ref_idx,j}.mean=mean(v.pod.catch{ref_idx,j}.ws.y_masked);
            %loop over all sources
            for i=1:v.pod.source.n
              %compute bias of this catchment
              v.pod.catch{i,j}.mean=mean(v.pod.catch{i,j}.ws.y_masked);
              %get the names of this solutions and enforce string replacement and cleaning
              latex_name{i}=upper(v.pod.source.names{i});
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
                latex_scale(5)* v.pod.catch{i,j}.ws.interp(...
                  v.pod.catch{ref_idx,j}.ws.t,'interp_over_gaps_narrower_than',v.plot_lines_over_gaps_narrower_than...
                ).stats2(...
                  v.pod.catch{ref_idx,j}.ws,'mode','corrcoeff','detrend','none'...
                )...                                                                           %corrcoeff
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
      if ~all(norm_stats(i,:)==0)
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
    function obj=plot_maps(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % add input arguments and plot metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{....
        plotting.default,...
        {...
          'plot_functional',    'eqwh', @gravity.isfunctional;...
          'plot_force',          false, @islogical;...
          'plot_time',              [], @(i) isdatetime(i) || isnumeric(i);...
          'plot_time_fmt','yyyy-mm-dd', @isdatetime;...
          'plot_spatial_step',       1, @num.isscalar;...
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

      %resolve implied meanings of plot_idx
      if isempty(v.plot_time)
        v.plot_time=v.pod.t;
      elseif isnumeric(v.plot_time)
        %translate negative indexes to point to the end extremity
        negative_idx=v.plot_time<0;
        v.plot_time(negative_idx)=numel(v.pod.t)+v.plot_time(negative_idx)+1;
        v.plot_time=sort(v.pod.t(v.plot_time));
      end

      for i=1:v.pod.source.n
        clear grid_now
        for it=1:numel(v.plot_time)
          %build filename
          filename=file.build(v.pod.file_root,...
            v.plot_functional,'map',v.pod.source.file_smooth{i},...
            strsplit(v.pod.source.names{i},' '),datestr(v.plot_time(it),v.plot_time_fmt),'png'...
          );
          %plot only if not done yet
          if ~file.exist(filename,v.plot_force)
            if ~exist('grid_now','var')
              grid_now=v.pod.source.dat{i...
                }.scale(v.plot_functional,'functional'...
                ).interp(v.plot_time...
                ).grid('spatial_step',v.plot_spatial_step...
              );
            end
            plotting.figure(v.varargin{:});
            grid_now.at(v.plot_time(it)).imagesc(...
              'cb_title',gravity.functional_label(v.plot_functional)...
            );
            plotting.enforce(v.varargin{:},...
              'plot_legend_location','none',...
              'plot_xlabel','none','plot_ylabel','none',...
              'plot_title',v.plot_title,...
              'plot_title_default',str.show({...
                upper(v.pod.source.names{i}),datestr(v.plot_time(it),v.plot_time_fmt),...
                v.pod.title_masking,...
                strjoin(cells.rm_empty({' ',v.pod.source.title_smooth{i}}),newline)...
              }));
            plotting.no_ticks
            plotting.save(filename,v.varargin{:})
          end
        end
      end
    end
    %% utils
    function d=grace_model(stem,varargin)
      if ~exist('stem','var') || isempty(stem)
%         stem='grace';
        stem='gracefo';
%         stem='gswarm.gracefo';
      end
      %NOTICE: this function is useful to rebuilt the grace (+ grace-fo) climatological model
      %end dates change
      switch stem
        case 'grace';  stop=datetime('2017-06-30');
        case {'gswarm.gracefo','gracefo'}
          stop=dateshift(datetime('now'),'end','month');
      end
      v=varargs.wrap('sources',{....
        {...
          'debug',     true,                      @islogical;...
          'force',     false,                     @islogical;...
          'start',     datetime('2002-04-01'),    @isdatetime;...
          'stop',      stop,                      @isdatetime;...
          'functional','eqwh',                    @ischar;...
          'smoothing', 350e3,                     @numeric;...
          'mode',      'all',                     @(i)ischar(i)||iscellstr(i);...
          'products',  {...
            [stem,'.sh.rl06.csr'];...
            [stem,'.sh.rl06.csr.pd'];...
            [stem,'.sh.rl06.csr.pd.ts'];...
          },                                      @iscellstr;...
          'plot_dir',file.orbdir('plot'),         @ischar;...
          'plot_title','',                        @ischar;... %there's a default title
        },...
      },varargin{:});
      %this is the root of all files saved in this method
      fnroot=fullfile(v.plot_dir,[stem,'.sh.rl06.csr.pd.ts-']);
      modes_plot={'C20.png','mean.png','std.png','rms.png'};
      modes_data={'d.mat','a.mat','b.mat'};
      modes_all=[modes_data(:);modes_plot(:)];
      %unwrap special modes
      switch cells.scalar(v.mode,'get')
      case 'all';  v.mode=modes_all;
      case 'data'; v.mode=modes_data;
      case 'plot'; v.mode=modes_plot;
      case 'done'
        d=true;
        for i=1:numel(modes_all)
          d=d && file.exist([fnroot,modes_all{i}]);
        end
        return
      end
      %create vector with relevant data product
      p=cellfun(@(i) dataproduct(i),v.products,'UniformOutput',false);
      %one does not usually want to ignore an explicit force flag, so be sure that failsafes are off
      p=cellfun(@(i) i.mdset('never_force',false,'always_force',false),p,'UniformOutput',false);
      %compute climatological model
      if cells.isincluded(v.mode,{'d.mat'})
        filename=[fnroot,'d.mat'];
        if ~file.exist(filename)
          %may need to init grace climatological model if force is true
          %this needs to be done separately from the Swarm processing because
          %otherwise the Swarm period is what is used in the regression: no bueno)
          d=datastorage(...
            'inclusive',true,...
            'start',v.start,...
            'stop', v.stop,...
            'debug',v.debug...
          );
          for i=1:numel(p)
            d=d.init(p{i},'force',v.force);
          end
          save(filename,'d'); disp(['saved ',filename])
        else
          load(filename,'d'); disp(['loaded ',filename])
        end
      end
      %assemble derived data
      if cells.isincluded(v.mode,{'a.mat'})
        filename=[fnroot,'a.mat'];
        if ~file.exist(filename)
          a={
            d.data_get_scalar(p{1}.dataname.append_field_leaf('signal')),'original' ;...
       pardecomp.join(d.data.(p{2}.codename)),                           're-modelled' ;... %mode evaluated at original t
            d.data_get_scalar(p{3}),                                     'model' ;...
          };
          save(filename,'a'); disp(['saved ',filename])
        else
          load(filename,'a'); disp(['loaded ',filename])
        end
      end
      %detrend derived data
      if cells.isincluded(v.mode,{'b.mat'})
        filename=[fnroot,'b.mat'];
        if ~file.exist(filename)
          b=a;
          for i=1:size(b,1)
            b{i,4}=b{i,1}.scale(v.functional,'functional').detrend.grid;
          end
          save(filename,'b'); disp(['saved ',filename])
        else
          load(filename,'b'); disp(['loaded ',filename])
        end
      end
      %plot C20
      if cells.isincluded(v.mode,{'C20.png'})
        filename=[fnroot,'C20.png'];
        if ~exist(filename,'file')
          tn=slr.load(dataproduct(p{1}).mdget('use_GRACE_C20'),'static_model',dataproduct(p{1}).mdget('static_model')).trim(v.start,v.stop).addgaps(days(45));
          plotting.figure;
          for i=1:size(a,1)
            a{i,1}.ts_C(2,0).addgaps(days(45)).plot;
          end
          tn.ts_C(2,0).plot
          %enforce it
          plotting.enforce(...
            'plot_legend',[a(:,2);{tn.descriptor}],...
            'plot_ylabel',gravity.functional_label(v.functional),...
            'plot_title',v.plot_title,...
            'plot_title_default','C20'...
          );
          plotting.save(filename,v.varargin{:})
          disp(['plotted ',filename])
        else
          disp(['plot already available: ',filename])
        end
      end
%this shows an odd-looking parabula
%       if cells.isincluded(v.mode,{'freq-stats','all'})
%         for i=1:size(a,1)
%           a{i,3}=a{i,1}.scale(v.functional,'functional');
%         end
%         %plot cumulative degree stats
%         for s={'mean','std','rms'}
%           plotting.figure;
%           c=0;legend_str=cell(0);
%           i=1;
%           [~,ts]=a{i,3}.(['cumd',s{1}]); ts.addgaps(days(45)).plot;
%           c=c+1; legend_str{c}='with C20 (TN-11)';
%           [~,ts]=a{i,3}.scale(v.smoothing,'gauss').(['cumd',s{1}]); ts.addgaps(days(45)).plot;
%           c=c+1; legend_str{c}=str.show({'with C20',v.smoothing/1e3,'km gauss'});
%           [~,ts]=a{i,3}.setC(2,0,0).(['cumd',s{1}]); ts.addgaps(days(45)).plot;
%           c=c+1; legend_str{c}='without C20';
%           i=2;
%           [~,ts]=a{i,3}.(['cumd',s{1}]); ts.addgaps(days(45)).plot;
%           c=c+1; legend_str{c}='p12 with C20';
%           plotting.enforce('plot_legend',legend_str,'plot_title',['cumulative degree ',s{1}],...
%             'plot_ylabel',gravity.functional_label(v.functional));
%         end
%       end

      %plot spatial stats (TODO: the mean should be zero! )
      stats={'mean','std','rms'};
      for s=1:numel(stats)
        if cells.isincluded(v.mode,[stats{s},'.png'])
          filename=[fnroot,stats{s},'.png'];
          if ~exist(filename,'file')
            plotting.figure;
            for i=1:size(b,1)
              [~,~,ts]=b{i,4}.(stats{s})('tot'); ts.addgaps(days(45)).plot;
            end
            %enforce it
            plotting.enforce(...
              'plot_legend',b(:,2),...
              'plot_ylabel',gravity.functional_label(v.functional),...
              'plot_title',v.plot_title,...
              'plot_title_default',['spatial ',stats{s}]...
            );
            plotting.save(filename,v.varargin{:})
            disp(['plotted ',filename])
          else
            disp(['plot already available: ',filename])
          end
        end
      end
    end
    function plotdeepoceanmask(l,plotname)
      a=simplegrid.unit(l*90,l*45);
      a=a.spatial_mask('deep ocean');
      c=gray(l*16);
      c=flipud(c(l*4:l*12,:));
      plotting.figure;
      a.imagesc;
      colormap(c)
      colorbar off
      plotting.enforce(v.varargin{:},...
        'plot_ylabel','none',...
        'plot_xlabel','none',...
        'plot_legend_location','none'...
      );
      plotting.no_ticks;
      saveas(gcf,plotname)
    end
    function d=grace_animation(varargin)
      %need global project variable (forces the user to think about the context of this analysis)
      global PROJECT
      if ~isfield(PROJECT,'stop_date')
        PROJECT.stop_date=datetime('now');
      end
      PROJECT.name='2020-06-02.grace_animation';
      v=varargs.wrap('sources',{....
        {...
          'debug',     true,                        @islogical;...
          'inclusive', true,                        @islogical;...
          'stop',      datetime(PROJECT.stop_date), @isdatetime;...
          'get_input_data',    false,               @islogical;... %NOTICE: consider turning this on to update all input data
          'c20model',  true,                        @islogical;...
          'end_product','gswarm.gracefo.sh.rl06.csr.pd.ts',@ischar;...
          'frame_rate', 30                          @(i) isnumeric(i) && isscalar(i);...
          'decimate_rate', 5                        @(i) isnumeric(i) && isscalar(i);...
          'filename','GRACE_animation.gif'          @ischar;
        },... 
      },varargin{:});
      %get input data
      if v.get_input_data;gswarm.get_input_data('GRACE');end
      %ensure C20 model is available
      if v.c20model && ~gswarm.c20model('done',file.orbdir('plot'))
        gswarm.c20model('plot',file.orbdir('plot'));
        gswarm.c20model('latex',file.orbdir('plot'));
      end
      %need to be sure grace model is available up until the stop time and including all available grace-fo data
      if ~gswarm.grace_model('gswarm.gracefo','mode','done')
        gswarm.grace_model('gswarm.gracefo','stop',v.stop);
      end
      %load the data
      d=datastorage('debug',v.debug,'inclusive',v.inclusive,...
          'start',datetime('2002-04-01'),...
          'stop',v.stop...
        ).init(v.end_product...
        ).data_get_scalar(v.end_product);
      %load the metadata
      a=dataproduct(v.end_product);
      %build animation
      close all
      n_frames=floor(d.nr_valid/v.decimate_rate);
      fig=plotting.figure('plot_size',200+[0,0,21,9]*50,'plot_visible','off');
      s.msg=['Rendering ',v.filename]; s.n=n_frames;
      im=cell([1,n_frames]);
      for i=1:n_frames
        d.at(d.t_masked(i*v.decimate_rate)...
          ).scale(a.mdget('plot_smoothing_degree'),a.mdget('plot_smoothing_method')...
          ).scale(a.mdget('plot_functional'),'functional'...
          ).grid('Nfactor',2 ...
          ).imagesc('center_resample',false,'cb_title',a.mdget('plot_cb_title'));
        plotting.enforce(...
          'plot_title',datestr(d.t_masked(i*v.decimate_rate),'yyyy'),...
          'plot_colormap',a.mdget('plot_colormap'),...
          'plot_caxis',[-1,1],...
          'plot_legend_location','none',...
          'plot_ylabel','none',...
          'plot_xlabel','none'...
        );
        axis off
        im{i}=frame2im(getframe(fig));
        
        s=time.progress(s,i);
      end
      A=cell([1,n_frames]);map=cell([1,n_frames]);
      for i = 1:n_frames
        [A{i},map{i}] = rgb2ind(im{i},256);
      end
      for i = 1:n_frames
        if i == 1
          imwrite(A{i},map{i},v.filename,'gif','LoopCount',Inf,'DelayTime',0);
        else
          imwrite(A{i},map{i},v.filename,'gif','WriteMode','append','DelayTime',1/v.frame_rate);
        end
      end
    end
    %% realizations
    function out=paper(varargin)
%             'swarm.sh.gswarm.rl01.individual';... %figs 6 - 8
%             'swarm.sh.gswarm.rl01.smooth750';...  %fig 9
%             'swarm.sh.gswarm.rl01.smooth1500';... %fig 10
%             'swarm.sh.gswarm.rl01.smooth0';...
%             'swarm.sh.gswarm.rl01.smooth300';...
%             'swarm.sh.gswarm.rl01.smooth-ocean1000';...
%             'swarm.sh.gswarm.rl01.smooth-ocean1500';...
%             'swarm.sh.gswarm.rl01.smooth-ocean3000';...
%             'swarm.sh.gswarm.rl01.smooth-ocean5000';...
%             'swarm.sh.gswarm.rl01.lowdeg';...
%             'swarm.sh.gswarm.rl01.catchments';...
%             'swarm.sh.gswarm.rl01.quality';...
%             'swarm.sh.gswarm.rl01.quality-land';...
%             'swarm.sh.gswarm.rl01.quality-land-deg40';...


% 'swarm.sh.gswarm.rl01.tropical';... used in the rebuttal to reviewer 1

      %NOTICE: don't add plotting.default as a source to the varargs.wrap (it overrides the metadata)
      v=varargs.wrap('sources',{...
        {...
          'debug',                   true, @logical;...
          'start', datetime('2013-09-01'), @isdatetime;...
          'stop',  datetime('2019-09-30'), @isdatetime;...
          'type',                   'all', @ischar;...
          'figures_dir',fullfile(file.scriptdir,'..','article','figures'), @ischar;...
          'force',                  false, @islogical;... %NOTICE: this affects datastorage.init
          'plot_force',             false, @islogical;... %NOTICE: this is not the same as plot_save_data (which is irrelevant here)
          'products',                  {}, @iscellstr;
        },...
      },varargin{:});
      switch v.type
      case 'process'
        assert(~isempty(v.products),'If ''type'' is ''process'', then need non-empty ''products''.')
        out=datastorage(...
          'start',v.start,...
          'stop', v.stop,...
          'inclusive',true,...
          'debug',v.debug...
        );
        for i=1:numel(v.products)
          out=out.init(v.products{i},...
            'force',v.force,...
            v.subset('plot').varargin{:}...
          );
        end
      case 'checking'
        out=cellfun(@(i) gswarm.paper(varargin{:},'type',i),{...
          'C20';...
          'C20-sources';...
          'TN11';...
          'tropics';...
          'smooth750.complete';...%used to confirm that detrending is a good idea
          'formerr';...
        },'UniformOutput',false);
      case 'all'
%           'deepoceanmask';...     %fig:method:deepoceanmass
%           'gracelowdegrees';...   %fig:method:climatmod
%           'individual';...        %fig:res:ind:cumdrms, fig:res:ind:rms_dmean, fig:res:ind:corrcoeff_dmean
%           'gauss-global0750';...  %fig:res:comb:750
%           'gauss-global1500';...  %fig:res:comb:1500
%           'gauss-ocean0750';...   %fig:res:comb:ocean:cumdrms,fig:res:comb:ocean:rms_dmean
%           'gauss-land0750';...    %fig:res:comb:land:cumdrms,fig:res:comb:land:rms_dmean
%           'gauss-ocean3000';...   %fig:res:comb:ocean:cumdrms,fig:res:comb:ocean:rms_dmean
%           'lowdeg';...            %fig:res:signal:lowdeg:{corrswarm,corrgrace,rmsswarm,C*}
%           'partitioning';...      %fig:res:comb:partitioning
        out=cellfun(@(i) gswarm.paper(varargin{:},'type',i),{...
          'catchments';...        %fig:res:signal:*,fig:res:signal:var:{swarm,grace}
        },'UniformOutput',false);
      case 'gauss'
        out=cellfun(@(i) gswarm.paper(varargin{:},'type',i),...
          {'gauss-global','gauss-land','gauss-ocean'},...
        'UniformOutput',false);
      case {'gauss-global','gauss-land','gauss-ocean'}
        if strcmp(v.type,'gauss-ocean')
          pl={[v.type,'0750'];[v.type,'1500'];[v.type,'3000'];[v.type,'5000']};
        else
          pl={[v.type,'0300'];[v.type,'0750'];[v.type,'1500'];[v.type,'3000']};
        end
        out=cellfun(@(i) gswarm.paper(varargin{:},'type',i),pl,'UniformOutput',false);
      case 'smoothing-methods'    %no longer used (saved to ../issues/smoothing-methods/)
        pl={};
        for m={'gauss','spline','trunc'}
          for l={'global','land','ocean'}
            if strcmp(l{1},'ocean')
              r={'0750','1500','3000','5000'};
            else
              r={'0300','0750','1500','3000'};
            end
            pl{end+1}=[m{1},'-',l{1},r{1}]; %#ok<AGROW>
          end
        end
        out=cellfun(@(i) gswarm.paper(varargin{:},'type',i),pl,'UniformOutput',false);
      case 'C20-source'           %auxiliar
        out=dataproduct('model.processing.replaceC20.submetadata').metadata.use_GRACE_C20;
      case 'C20'                  %checking only
        plotfilename=fullfile(v.check_figs_dir,'checks','C20.png');
        if ~file.exist(plotfilename)
          plotting.figure;
          slr.load(gswarm.paper('type','C20-source')).plot;
          plotting.save(plotfilename,v.varargin{:});
        else
          out=[];
        end
      case 'C20-sources'          %checking only
        plotfilename=fullfile(v.figures_dir,'checks','C20-sources.png');
        if ~file.exist(plotfilename)
%         version_list={'GSFC-7day','GSFC','TN-14','TN-11','TN-07'};
          version_list={'GSFC-7day','TN-11','TN-11-C20-model'};
          plotting.figure;
          cellfun(@(i) slr.load(...
            i...
          ).plot(...
            'method','timeseries','degrees',2,'orders',0 ...
          ),version_list);
          plotting.enforce(...
            'plot_line_color','spiral',...
            'plot_legend',version_list...
          );
          plotting.save(plotfilename,v.varargin{:});
        else
          out=[];
        end
      case 'TN11'                 %checking only (not used anymore)
        plotfilename=fullfile(v.figures_dir,'checks','C20.TN11.png');
        if ~file.exist(plotfilename)
          slr.load('TN-11',...
            'start',time.zero_date,...
            'stop',datetime('2019-10-31')...
          ).plot(...
            'method','timeseries','degrees',2,'orders',0 ...
          )
          if str.none(v.plot_title); title(''); end
          plotting.save(plotfilename,v.varargin{:});
        else
          out=[];
        end
      case 'deepoceanmask'
        plotfilename=fullfile(v.figures_dir,'deepoceanmask.png');
        if ~file.exist(plotfilename,v.plot_force,true)
          gswarm.plotdeepoceanmask(4,plotfilename)
        end
        out=[];
      %NOTICE: the flag force will recompute the *timeseries* of the parametric decomposition even if never_force is true
      %NOTICE: use 'type','pd' if you want to recompute the parametric decomposition (consider 'force',true)
      case {'gracelowdegrees',...
            'gracefolowdegrees'}
        %resolve internal type
        type=strrep(v.type,'lowdegrees','');
        %define relevant products (the first product is assumed to have the 'signal' field)
        p=dataproduct.array({...
          'grace.sh.rl06.csr';...
          [type,'.sh.rl06.csr.pd.ts'];...
        });
        %define zonal degrees to plot
        degree_list=2:4;
        %one does not usually want to ignore an explicit force flag, so be sure that failsafes are off
        p=cellfun(@(i) i.mdset('never_force',false,'always_force',false),p,'UniformOutput',false);
        %build plot name list
        plotfilenames=arrayfun(@(i) fullfile(...
          v.figures_dir,[type,'.C',num2str(i),'0.png']...
        ),degree_list,'UniformOutput',false);
        %check if plots are already there
        if all(file.exist(plotfilenames,v.plot_force,true))
          out=[];
          return
        end
        %define datastorage and load data
        out=datastorage(...
          'start',datetime('2002-04-01'),...
          'stop', v.stop,...
          'debug',v.debug...
        );
        for i=1:numel(p); out=out.init(p{i},'force',v.force); end
        %make pretty plots
        grs=out.data_get_scalar(p{1}.dataname.set_field_path('signal'));
        grm=out.data_get_scalar(p{2});
        for i=1:numel(degree_list)
          plotting.figure(v.varargin{:});
          grsp=grs.ts_C(degree_list(i),0).addgaps(days(45)).plot('line',{'-o'});
          grmp=grm.ts_C(degree_list(i),0).plot;
          plotting.enforce(v.varargin{:},...
            'plot_legend',{...
              str.rep(grsp.legend{1},['C',num2str(degree_list(i)),',0'],'CSR RL06');...
              str.rep(grmp.legend{1},['C',num2str(degree_list(i)),',0'],'12-parameter model');...
            },...
            'plot_title','none',...
            'plot_ylabel',['C_{',num2str(degree_list(i)),',0} [ ]']...
          );
          saveas(gcf,plotfilenames{i})
        end
      %this is used to recompute the pd models
      %NOTICE: with force,'true' this will recompute parametric decomposition even if never_force is true
      %NOTICE: this is similar to above, except that the grace/gracefo type is in agreement
      case 'pd'                   %checking only (recompute climatological models)
        %define relevant products (the first product is assumed to have the 'signal' field)
        p=dataproduct.array({...
          'grace.sh.rl06.csr.pd';...
          'grace.sh.rl06.csr.pd.ts';...
          'gracefo.sh.rl06.csr.pd';...
          'gracefo.sh.rl06.csr.pd.ts';...
        });
        %one does not usually want to ignore an explicit force flag, so be sure that failsafes are off
        p=cellfun(@(i) i.mdset('never_force',false,'always_force',false),p,'UniformOutput',false);
        %init data
        out=datastorage(...
          'start',datetime('2002-04-01'),...
          'stop', v.stop,...
          'inclusive',true,...
          'debug',v.debug...
        );
        for i=1:numel(p); out=out.init(p{i},'force',v.force); end
      case 'pd-plot'              %checking only
        dn_list={...
          'grace.sh.rl06.csr/signal',...
          'grace.sh.rl06.csr.pd.ts',...
          'gracefo.sh.rl06.csr/signal',...
          'gracefo.sh.rl06.csr.pd.ts'...
        };
        %define degrees/orders to plot
        degree_list=2:4;
         order_list=zeros(size(degree_list));
        %build plot name
        plotfilenames=arrayfun(...
          @(i,j) fullfile(v.figures_dir,'checks',['grace.vs.gracefo.pd.C',num2str(i),num2str(j),'.png']),...
        degree_list,order_list,'UniformOutput',false);
        %check if plots are already there
        if ~all(file.exist(plotfilenames,v.plot_force,true))
          out=datastorage(...
            'start',v.start,...
            'stop', v.stop,...
            'inclusive',true,...
            'debug',v.debug...
          );
          g=cell(size(dn_list));
          for i=1:numel(dn_list)
            out=out.init(dn_list{i},'force',v.force);
            g{i}=out.data_get_scalar(dn_list{i});
          end
          for i=1:numel(degree_list)
            plotting.figure;
            for gi=1:numel(g)
              g{gi}.plot('method','timeseries','degrees',degree_list(i),'orders',order_list(i));
            end
            plotting.enforce(...
              'plot_legend',dn_list,...
              'plot_title',['C',num2str(degree_list(i)),num2str(order_list(i))]...
            );
            saveas(gcf,plotfilenames{i})
          end
        end
      case 'partitioning'
        functional='eqwh';
        gap_size=days(120);
        filenameroot=fullfile(v.figures_dir,'checks','partitioning');
        datafilename=[filenameroot,'.mat'];
        if ~file.exist(datafilename)
          out=gswarm.paper(...
            'type','process',...
            'products',{...
              'swarm.sh.gswarm.rl01.res.gracefo.sh.rl06.csr',...
              'gracefo.sh.rl06.csr',...
              'swarm.sh.gswarm.rl01'...
            }...
          );
          str.say('gathering')
          d0={...
            out.data.swarm_sh_gswarm_rl01_res_gracefo_sh_rl06_csr;...
            out.data.gracefo_sh_rl06_csr.signal;...
            out.data.swarm_sh_gswarm_rl01.signal;...
          };
          str.say('preprocessing')
          t=d0{2}.t_masked;
          di=cellfun(...
            @(i) i.set_lmax(20 ...
            ).scale(750e3,'gauss'...
            ).scale(functional,'functional'...
            ).interp(t,'interp_over_gaps_narrower_than',gap_size),...
          d0,'UniformOutput',false);
          str.say('Spatial masking')
          d.gm.o={...
            di{1};...
            di{1}.spatial_mask(    'tropical');...
            di{1}.spatial_mask('non-tropical');...
            di{1}.spatial_mask(        'land');...
            di{1}.spatial_mask(  'deep ocean');...
          };
          d.gr.o={...
            di{2};...
            di{2}.spatial_mask(    'tropical');...
            di{2}.spatial_mask('non-tropical');...
            di{2}.spatial_mask(        'land');...
            di{2}.spatial_mask(  'deep ocean');...
          };
          d.sw.o={...
            di{3};...
            di{3}.spatial_mask(    'tropical');...
            di{3}.spatial_mask('non-tropical');...
            di{3}.spatial_mask(        'land');...
            di{3}.spatial_mask(  'deep ocean');...
          };
          str.say('Detrending and computing DAS')
          d.gm.n='\Delta(Swarm,GRACE)';
          d.gr.n='GRACE';
          d.sw.n='Swarm';
          for j={'gm','gr','sw'}
            dnow=cellfun(@(i) i.detrend.addgaps(gap_size),d.(j{1}).o,'UniformOutput',false);
            d.(j{1}).d=cellfun(@(i) i.cumdas,dnow,'UniformOutput',false);
            d.(j{1}).t=cellfun(@(i) i.t     ,dnow,'UniformOutput',false);
          end
          str.say('Created data file',datafilename)
          save(datafilename,'d')
        else
          load(datafilename,'d')
        end
        plotfilename=[filenameroot,'.all.png'];
        if ~file.exist(plotfilename)
          str.say('Plotting',plotfilename)
          plotting.figure; hold on
          legend_str={};
          type_str={'global','tropical','non-tropical','land','deep ocean'};
          fn=fieldnames(d);
          mf={'o','+','x'};
          for f=1:numel(fn)
            for i=1:numel(d.(fn{f}).t)
              plot(d.(fn{f}).t{i},d.(fn{f}).d{i}(:,end),[mf{f},'-'])
              legend_str{end+1}=[d.(fn{f}).n,' ',type_str{i}]; %#ok<AGROW>
            end
          end
          plotting.enforce(...
            'plot_title','none',...
            'plot_ylabel',gravity.functional_label(functional),...
            'plot_xdateformat', 'yy/mm',...
            'plot_legend_sorting', 'none',...
            'plot_legend_include_smoothing', false,...
            'plot_fontsize_axis', 28,...
            'plot_fontsize_title', 32,...
            'plot_fontsize_label', 28,...
            'plot_fontsize_legend', 28,...
            'plot_legend_location','eastoutside',...
            'plot_legend',legend_str...
          );
          plotting.save(plotfilename,v.varargin{:});
        end
        plotfilename=fullfile(v.figures_dir,[v.type,'.res.png']);
        if ~file.exist(plotfilename)
          str.say('Plotting',plotfilename)
          plotting.figure; hold on
          legend_str={};
          type_str={'global','tropical','non-tropical','land','deep ocean'};
          fn={'gm'};
          for f=1:numel(fn)
            for i=1:numel(d.(fn{f}).t)
              plot(d.(fn{f}).t{i},d.(fn{f}).d{i}(:,end))
              legend_str{end+1}=type_str{i}; %#ok<AGROW>
            end
          end
          plotting.enforce(...
            'plot_title','none',...
            'plot_ylabel',gravity.functional_label(functional),...
            'plot_xdateformat', 'yy/mm',...
            'plot_legend_sorting', 'none',...
            'plot_legend_include_smoothing', false,...
            'plot_fontsize_axis', 28,...
            'plot_fontsize_title', 32,...
            'plot_fontsize_label', 28,...
            'plot_fontsize_legend', 28,...
            'plot_legend',legend_str...
          );
          plotting.save(plotfilename,v.varargin{:});
        end
        out=d;
      case 'res'                  %processing only (called implicitly where relevant)
        out=gswarm.paper(varargin{:},'type','res.gracefo.sh.rl06.csr');
      case 'formerr'              %checking only
        out=gswarm.paper(varargin{:},'type','individual.formerr');
      otherwise
        out=gswarm.paper(...
          'type','process',...
          'products',{['swarm.sh.gswarm.rl01.',v.type]}...
        );
      end
    end
    %% Swarm-ITT
    function d=quality(varargin)
      global PROJECT
      %NOTICE: this is used to produce the plots in ~/data/gswarm/dissemination/quality
      %workflow:
      % - you need to delete the last data file of (if not run from a dedicated dir):
      %   - swarm.sh.gswarm.rl01
      %   - swarm.sh.gswarm.rl01.err
      % - check the last entries in the time domain correspond to the most recent models

      %parse inputs
      v=varargs.wrap('sources',{....
        {...
          'mask',                 'deep ocean', @ischar;...
          'debug',                        true, @ischar;...
          'start',datetime(PROJECT.start_date), @isdatetime;...
          'stop', datetime(PROJECT.stop_date ), @isdatetime;...
          'smoothing',                   750e3, @isnumeric;...
          'functional',                'geoid', @ischar;...
          'quality_dir',fullfile(getenv('HOME'),'data','gswarm','dissemination','quality')...
                                              , @ischar;...
          'git_ci',                       true, @str.islogical;...
        },...
      },varargin{:});
      %changes will be saved in git at the end, make sure git is ok in the quality_dir
      assert(file.git('isuptodate',v.quality_dir),['There are changes in ',v.quality_dir,...
        ' that need to be committed before this method can continue.'])
      %define Swarm products
      p.err = dataproduct('swarm.sh.gswarm.rl01.err');
      p.swm = dataproduct(p.err.sources(1)); %swarm.sh.gswarm.rl01
      p.grm = dataproduct(p.err.sources(2)); %gswarm.gracefo.sh.rl06.csr.pd.ts
      p.grs = dataproduct(dataproduct(p.grm.sources(1)).sources(1)); %gswarm.gracefo.sh.rl06.csr
      %sanity on the C20 data
      assert(strcmp( p.grs.metadata.use_GRACE_C20, p.swm.metadata.use_GRACE_C20),...
        ['C20 model of products ',p.grs.name,' and ',p.swm.name,' are not the same.'])
      %retrienvee C20 data name
      C20_name=p.grs.metadata.use_GRACE_C20;
      %report input models
      [~,f]=fileparts(strrep(p.swm.metadata.wildcarded_filename,'*_',''));
      str.say(' --- Input combined models are: ',f)
      %filename stem
      f=fullfile(v.quality_dir,...
        strjoin({...
          f,datestr(now,'yyyy-mm-dd'),[datestr(v.start,'yyyymm'),'-',datestr(v.stop,'yyyymm')],...
          num2str(v.smoothing/1e3)...
        },'.')...
      );
      str.say(' --- File name stem is: ',f)

      %show the grace model
      fn={[f,'.grace-C20.png'],[f,'.grace-C30.png']};
      if any(~file.exist(fn))
        start=datetime('2001-04-01');
        str.say(' --- plotting GRACE model from',start,'to',v.stop)
        d=datastorage(...
          'start',start,...
          'stop', v.stop,...
          'debug',v.debug...
        ).init(p.grs).init(p.grm);
        %make pretty plots
        if ~file.exist(fn{1})
          plotting.figure(varargin{:});
          grs=d.data_get_scalar(p.grs.dataname).ts_C(2,0).plot('zeromean',true);
          grm=d.data_get_scalar(p.grm.dataname).ts_C(2,0).plot('zeromean',true);
          c20=slr.load(               C20_name).ts_C(2,0).plot('zeromean',true);
          plotting.enforce('plot_legend',{...
            str.rep(grs.legend{1},'C2,0','GRACE data');...
            str.rep(grm.legend{1},'C2,0','GRACE model');...
            str.rep( c20.legend{1},'C20',[C20_name,' (w/ static)']);...
          },'plot_title','C20','plot_ylabel','non-dim');
          saveas(gcf,fn{1})
          str.say('plotted',fn{1})
        end
        if ~file.exist(fn{2})
          plotting.figure(varargin{:});
          grs=d.data_get_scalar(p.grs.dataname).ts_C(3,0).plot('zeromean',true);
          grm=d.data_get_scalar(p.grm.dataname).ts_C(3,0).plot('zeromean',true);
          plotting.enforce('plot_legend',{...
            str.rep(grs.legend{1},'C3,0','GRACE data');...
            str.rep(grm.legend{1},'C3,0','GRACE model');...
          },'plot_title','C30','plot_ylabel','non-dim');
          saveas(gcf,fn{2})
          str.say('plotted',fn{2})
        end
      end

      %show the GRACE and Swarm data
      fn={[f,'.swarm-C20.png'],[f,'.swarm-C30.png']};
      if any(~file.exist(fn))
        start=datetime('2013-10-01');
        str.say(' --- plotting GRACE and Swarm data from',start,'to',v.stop)
        d=datastorage(...
          'start',start,...
          'stop', v.stop,...
          'debug',v.debug...
        ).init(p.grs).init(p.grm).init(p.swm);
        %make pretty plots
        if ~file.exist(fn{1})
          plotting.figure(varargin{:});
          swp=d.data_get_scalar(p.swm.dataname).ts_C(2,0).addgaps(days(45)).plot('zeromean',false);
          gsp=d.data_get_scalar(p.grs.dataname).ts_C(2,0).addgaps(days(45)).plot('zeromean',false);
          gmp=d.data_get_scalar(p.grm.dataname).ts_C(2,0).addgaps(days(45)).plot('zeromean',false);
          c20=slr.load(               C20_name).ts_C(2,0).addgaps(days(45)).plot('zeromean',false);
          %NOTICE: if C20_name does not contain '-model', there are more legend entries than plotted lines
          plotting.enforce('plot_legend',{...
            str.rep(swp.legend{1},'C2,0','Swarm (repl.)');...
            str.rep(gsp.legend{1},'C2,0','GRACE data');...
            str.rep(gmp.legend{1},'C2,0','GRACE model (p12)');...
            strrep(C20_name,'-model','');...
            C20_name;...
          },'plot_title','C20','plot_ylabel','non-dim');
          saveas(gcf,fn{1})
          str.say('plotted',fn{1})
        end
        if ~file.exist(fn{2})
          plotting.figure(varargin{:});
          sw=d.data_get_scalar(p.swm.dataname).ts_C(3,0).addgaps(days(45)).plot;
          gs=d.data_get_scalar(p.grs.dataname).ts_C(3,0).addgaps(days(45)).plot;
          gm=d.data_get_scalar(p.grm.dataname).ts_C(3,0).addgaps(days(45)).plot;
          plotting.enforce('plot_legend',{...
            str.rep(sw.legend{1},'C3,0','Swarm (repl.)');...
            str.rep(gs.legend{1},'C3,0','GRACE data');...
            str.rep(gm.legend{1},'C3,0','GRACE model');...
          },'plot_title','C30','plot_ylabel','non-dim');
          saveas(gcf,fn{2})
          str.say('plotted',fn{2})
        end
      end

      %show the GRACE and Swarm data
      fn=strjoin({f,v.functional,'quality','png'},'.');
      if ~file.exist(fn)
        start=datetime('2013-10-01');
        str.say(' --- plotting quality from',start,'to',v.stop)
        d=datastorage(...
          'inclusive',true,... %this needs to be true, otherwise the last month may be gone
          'start',start,...
          'stop', v.stop,...
          'debug',v.debug...
        ).init(p.swm).init(p.err);
        str.say(' --- Getting gravity:')
        g =d.data_get_scalar(p.err.dataname.append_field_leaf('signal'));
        gs=d.data_get_scalar(p.swm.dataname.append_field_leaf('signal'));
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
        saveas(gcf,fn)
        str.say(' --- Exporting  the time series:')
        cts.export(strrep(fn,'.png','.dat'),'ascii')
      end

      %save new data in git
      if v.git_ci
        str.say(' --- Committing new data:')
        disp(file.git('add',[v.quality_dir,'/*']))
        disp(file.git('ci' ,[v.quality_dir,'/*'],['updated quality data up until ',datestr(v.stop)]))
      end

    end
    function get_input_data(mode)
      switch mode
      case('GRACE')
        %NOTICE: the download-l2.sh script iterates over specific years (currently 2021)
        disp('Downloading GRACE data')
        file.system('./download-l2.sh CSR 06','disp',true,'cd',grace.dir('l1b'));
      case ('Swarm')
        disp('Downloading Swarm data')
        file.system('~/bin/rsyncf.sh aristarchos --no-l2r','disp',true,'cd',gswarm.dir('data'));
      case('all')
        gswarm.get_input_data('GRACE')
        gswarm.get_input_data('Swarm')
      otherwise
        error(['Cannot handle ''mode'' with value ''',mode,'''.'])
      end
    end
    function out=c20model(mode,plot_dir,version)
      if ~exist('version','var') || isempty(version)
        version=dataproduct('model.processing.replaceC20.submetadata').metadata.use_GRACE_C20;
      end
      %check if this is a C20 model
      model=contains(version,'-model');
      plotfile=fullfile(plot_dir,['C20_',version,'.png']);
      %document the C20 model
      switch mode
      case 'done'
        out=file.exist(plotfile);
        if model
          out=out && file.exist(strrep(plotfile,'.png','.tex'));
        end
      case 'clear'
        out=delete(gswarm.c20model('filename',plot_dir));
      case 'filename'
        out=plotfile;
      case 'plot'
        if ~gswarm.c20model('done',plot_dir)
          %this updates the coefficient time series, to make sure it has recent-enough data
          slr.graceC20('mode','set','version',strrep(version,'-model',''));
          if model
            %this updates the model itself
            slr.graceC20('mode','set','version',strrep(version,'-model',''));
            out=slr.graceC20('mode','model-plot','version',version);
          else
            out=slr.graceC20('mode','plot-all','version',unique({'GSFC-7DAY','GSFC','TN-14','TN-11',upper(version)}));
          end
          plotting.save(gswarm.c20model('filename',plot_dir))
        else
          out=[];
        end
      case 'latex'
        latexfile=strrep(plotfile,'.png','.tex');
        if model
          out=slr.graceC20('mode','model-list-tex','version',version);
          if ~file.exist(latexfile)
            file.strsave(latexfile,out);
            disp(['Saved C20 coefficients to ',latexfile])
          end
        else
          file.strsave(latexfile,'');
          disp(['Saved empty string (',version,' is not a C20 model) to ',latexfile])
        end
      otherwise
        error(['Unknown mode ''',mode,'''.'])
      end
    end
    function d=precombval(varargin)
      %NOTICE: this method produces the plots in the dir defined in the project.yaml file, which are
      %        needed to produce the pre-combination report in the 'report' dir at the same location,
      %        usually: ~/data/gswarm/analyses/<date>-precombval
      %process
      d=gswarm.production(...
        'products',  {...
          'swarm.sh.gswarm.rl01.individual';...
          'swarm.sh.gswarm.rl01.individual.maps';...
        },...
        'get_input_data',false,...
        'plot_pause_on_save',false,...
        varargin{:}...
      );
    end
    function d=validation(varargin)
      %NOTICE: this method produces the plots in the dir defined in the project.yaml file,
      %        which are needed to produce the pre-combination report in the 'report' dir
      %        at the same location, usually: ~/data/gswarm/analyses/<date>-validation

      %WORKFLOW Workflow of the TYPE={precombval|validation} report:
      %WORKFLOW 1.  OPTIONAL: plug the thumb drive in (no need to mount-disk.sh)
      %WORKFLOW 2.  sure everything is in synch.sh:
      %WORKFLOW     2.1: in ~/data/gswarm/analyses/report-TYPE/
      %WORKFLOW     2.2: rsyncf.sh remotes-file=~/data/gswarm/rsyncf.list --no-l2r --dry-run aristarchos
      %WORKFLOW 3.  create a new TYPE dir: ~/data/gswarm/analyses/new-analysis.sh TYPE
      %WORKFLOW 4.  cd to the orb dir in the new TYPE dir (shown by new-analysis.sh
      %WORKFLOW     script) and check the stop_date in project.yaml is the last day of the
      %WORKFLOW     last available model (the new-analysis.sh script automatically updates
      %WORKFLOW     this but better check if it's OK)
      %WORKFLOW 5.  fire up matlab and look at the gswarm.TYPE method:
      %WORKFLOW     5.1: check if all products are being used, some may be commented
      %WORKFLOW     5.2: check if the 'get_input_data' option is true:
      %WORKFLOW         5.2.1: the swarm data is downloaded from aristarchos (need
      %WORKFLOW                rsyncf.sh and ~/data/gswarm/rsyncf.list)
      %WORKFLOW         5.2.2: the GRACE data is downloaded from PODACC (need
      %WORKFLOW                ~/data/grace/download-l2.sh, which
      %WORKFLOW                iterates over specific years, currently 2021)
      %WORKFLOW         5.2.2: NOTICE: when doing tests, it's quicker to set 'get_input_data' to false.
      %WORKFLOW     5.3: if TYPE=validation, check if the 'git_ci' option is true:
      %WORKFLOW         5.3.1: after the swarm data is processed, the quality is computed in
      %WORKFLOW                the gswarm.quality method, where it is added to the git repo
      %WORKFLOW         5.3.2: NOTICE: when doing tests, it's best to set 'git_ci' to false.
      %WORKFLOW     5.4: check that the C20 data is updated and the model evaluated at the
      %WORKFLOW          last 3 months
      %WORKFLOW         5.4.1: The easiest way to be sure is to run:
      %WORKFLOW                'gswarm.c20model('plot',file.orbdir('plot'))'
      %WORKFLOW         5.4.2: For TYPE=precombval, a TN-14 model is used, so this is a good 
      %WORKFLOW                opportunity to send the email to Bryant Loomis and ask for 
      %WORKFLOW                the updated weekly C20 data.
      %WORKFLOW         5.4.3: For TYPE=validation, make sure the data Bryant sent is saved as
      %WORKFLOW                ~/data/gswarm/analyses/<date>-validation/orb/auxiliary/GSFC_SLR_C20_7day.txt
      %WORKFLOW     5.5: run the gswarm.TYPE method and keep an eye the last epoch of the
      %WORKFLOW          data as it is being loaded, it has to be the same as the last
      %WORKFLOW          available month; otherwise the analysis is incomplete
      %WORKFLOW     5.6: things that may go wrong:
      %WORKFLOW         5.6.1: the C20 data is not up-to-date
      %WORKFLOW         5.6.2: IfG releases a new version of their models and the metadata
      %WORKFLOW                was not updated to that
      %WORKFLOW         5.6.3: AIUB names the models incorrectly or does not compress them
      %WORKFLOW         5.6.4: You changed a matlab class and the *.mat files in the GRACE
      %WORKFLOW                L2/AIUB/ASU/IfG/OSU data dirs are now outdated (you can tell
      %WORKFLOW                this is the case when the error happens only on the first new
      %WORKFLOW                models); just delete the offending *.mat files:
      %WORKFLOW                rm -fv ~/data/gswarm/*/gravity/*.mat ~/data/grace/L2/CSR/RL06/*.mat
      %WORKFLOW                and re-import everything (by simply re-running gswarm.TYPE).
      %WORKFLOW         5.6.5: Some analysis start in 2002-04 instead of 2016-01, particularly
      %WORKFLOW                the product gswarm.swarm.validation.unsmoothed. Not sure why
      %WORKFLOW                this is happening but running it another time fixes the problem.
      %WORKFLOW 6.  go through the report and update all %NEEDS UPDATING lines (some are
      %WORKFLOW     automatically updated by new-analysis.sh)
      %WORKFLOW     6.1: There are plots that refer to the most recent months and cannot be
      %WORKFLOW          named automatically, look for '%NEEDS UPDATING (MAPS)' and run
      %WORKFLOW          ./ls-missing-figures.sh to see which plots are being wrongly picked.
      %WORKFLOW 7.  compile it and compare this report with the previous one
      %WORKFLOW 8.  go to the report dir in the new TYPE dir and make sure everything
      %WORKFLOW     is in synch.sh (don't synch %NEEDS UPDATING lines)
      %WORKFLOW 9.  Submit changes to GitHub:
      %WORKFLOW     9.1: Run the script link-dup-metadatafiles.sh in directory
      %WORKFLOW          ~/data/gswarm/analyses/<date>-TYPE/orb/metadata/ with arguments:
      %WORKFLOW          source=<previous date>-TYPE sink=<date>-TYPE [echo]
      %WORKFLOW     9.2: check the following git repos to make sure things are up to date (should
      %WORKFLOW          be open in smerge):
      %WORKFLOW          - ~/data/gswarm/analyses/<date>-TYPE/orb/
      %WORKFLOW          - ~/data/gswarm/analyses/
      %WORKFLOW 10. put the data into aristarchos (remove --dry-run, as usual):
      %WORKFLOW     ~/data/gswarm/rsync.local2remote-subset.sh --delete --dry-run
      %WORKFLOW
      %WORKFLOW If this is a precombval, then mail it to colleagues and you're done.
      %WORKFLOW
      %WORKFLOW 11. run publish.sh [echo] (inside the dir of the latest report)
      %WORKFLOW 12. email the report to colleagues
      %WORKFLOW
      %WORKFLOW     (wait for email reponses)
      %WORKFLOW
      %WORKFLOW 13. make sure data is in aristarchos (remove --dry-run, as usual):
      %WORKFLOW     ~/data/gswarm/rsync.local2remote-subset.sh --delete --dry-run
      %WORKFLOW 14. login to aristarchos and:
      %WORKFLOW     cd /home/gswarm/data/dissemination
      %WORKFLOW     tail *list && ./op.sh get-L1B-GPS get-L1B-ATT force && tail *list
      %WORKFLOW     ./op.sh update-source update-raw
      %WORKFLOW     ./op.sh make -val=yyyy-mm (run this for as many months as needed)
      %WORKFLOW     [...]
      %WORKFLOW     follow the instructions, additional commands needed.
      %WORKFLOW
      %WORKFLOW     After sending the email to ESA and ICGEM (as reported in the instructions):
      %WORKFLOW 15  update the dissemination dir git repo:
      %WORKFLOW     git st && git add *; git ci -m 'added new models'
      %WORKFLOW 16. update local data (remove --dry-run, as usual):
      %WORKFLOW     ~/data/gswarm/rsync.remote2local-subset.sh --delete --dry-run
      %WORKFLOW 17. plug in the data disk, mount it and sync the thumbs drive:
      %WORKFLOW     mount-disk.sh data
      %WORKFLOW     ~/media/data/data/rsync.thumb.sh reverse --dry-run

      %TODO on next release(s):
      % - remove trends from gswarm.swarm.validation.maps (also plot trends)
      % - ensure EWH axis are consistent for global, land and ocean plots
      % - fix titles of maps contain 'AIUB AIUB' and 'IFG IFG'
      
      %TODO: gswarm.swarm.validation.unsmoothed is broken in the first run, the time series 
      %      starts in April 2002
      
      %produce plots for the report
      d=gswarm.production(...
        'products',  {...
          'gswarm.swarm.validation.land';...
          'gswarm.swarm.validation.deepocean';...
          'gswarm.swarm.validation.smoothed';...
          'gswarm.swarm.validation.pardecomp';...
          'gswarm.swarm.validation.catchments';...
          'gswarm.swarm.validation.maps';...
          'gswarm.swarm.validation.unsmoothed';...
        },...
        'get_input_data',false,... #this is usually more problems that what it's worth; just copy the data manually
        varargin{:}...
      );
      %export quality metrics
      gswarm.quality(...
        'git_ci',true,...
        varargin{:}...
      );
    end
    function out=production_date(mode)
      global PROJECT
      assert(strcmp(mode,'start') || strcmp(mode,'stop'),...
        'input ''mode'' must either ''start'' or ''stop''.')
      mode=[mode,'_date'];
      if isfield(PROJECT,mode)
        out=datetime(PROJECT.(mode));
      else
        out=gswarm.(mode);
      end
    end
    function d=production(varargin)
      %need global project variable (forces the user to think about the context of this analysis)

      %NOTICE: this method expects some input arguments, notably:
      % - products (cellstr)

      %parse input args
      %NOTICE: gracefo.sh.rl06.csr.ld.ts has metadata never_force set as true (usually!)
      %        so 'force' as true will only reload the Swarm individual models
      %NOTICE: 'inclusive' can be false, because the GRACE data is only used to derive
      %        gracefo.sh.rl06.csr.ld.ts, separately
      v=varargs.wrap('sources',{....
        {...
          'products',  {},                          @iscellstr;...
          'start',     gswarm.production_date('start'),@isdatetime;...
          'stop',      gswarm.production_date('stop') ,@isdatetime;...
          'debug',     true,                        @islogical;...
          'inclusive', true,                        @islogical;...
          'force',     false,                       @islogical;... %this affects datastorage.init
          'force_d',   false,                       @islogical;... %this affects load(datafilename,'d') if force is false
          'get_input_data',false,                   @islogical;... %NOTICE: consider turning this on to update all input data
          'c20model',  true,                        @islogical;...
         'grace_model',true,                        @islogical;...
        },...
      },varargin{:});
      %check if all needed arguments are available
      required_args={'products'};
      for i=1:numel(required_args)
        assert(~isempty(v.(required_args{i})),['need argument ''',required_args{i},'''.'])
      end
      %get input data
      if v.get_input_data;gswarm.get_input_data('all'); end
      %ensure C20 model is available
      if v.c20model && ~gswarm.c20model('done',file.orbdir('plot'))
        gswarm.c20model('plot',file.orbdir('plot'));
        gswarm.c20model('latex',file.orbdir('plot'));
      end
      %need to be sure grace model is available up until the stop time and including all available grace-fo data
      if v.grace_model && ~gswarm.grace_model('gswarm.gracefo','mode','done')
        gswarm.grace_model('gswarm.gracefo','stop',v.stop);
      end
      %create vector with relevant data product
      p=cellfun(@(i) dataproduct(i),v.products,'UniformOutput',false);
      %define data filename
      datafilename=fullfile(file.orbdir('plot'),'d.mat');
      %load data if already available
      if file.exist(datafilename) && ~v.force_d && ~v.force
        str.say('Loading analysis data from',datafilename)
        load(datafilename,'d')
      else
        d=datastorage('debug',v.debug,'inclusive',v.inclusive,...
          'start',v.start,...
          'stop',v.stop);
        for i=1:p{1}.nr_sources
          d=d.init(p{1}.sources(i),'force',v.force);
        end
        file.ensuredir(datafilename,true);
        save(datafilename,'d')
      end
      %plot it
      tstart=tic;
      for i=1:numel(p)
        d.init(p{i},v.varargin{:});
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