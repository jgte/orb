classdef grace
  %static
  properties(Constant)
    %default value of some internal parameters
    default_list=struct(...
      'calpar_csr_defs',struct(... %calpar name; plot scale; 
        'AC0X_'     ,1e-7,...
        'AC0Z_'     ,5e-7,...
        'AC0Y_____1',1e-5,...
        'AC0Y_____2',1e-5,...
        'AC0Y_____3',1e-5,...
        'AC0Y_____4',2e-5,...
        'AC0Y_____5',2e-5,...
        'AC0Y_____6',5e-5,...
        'AC0Y_____7',5e-5,...
        'AC0Y_____8',5e-5,...
        'AC0XD'     ,1e-7,...
        'AC0ZD'     ,1e-6,...
        'AC0YD____1',1e-5,...
        'AC0YD____2',5e-5,...
        'AC0YD____3',5e-5,...
        'AC0YD____4',1e-4,...
        'AC0YD____5',1e-4,...
        'AC0YD____6',1e-4,...
        'AC0YD____7',1e-4,...
        'AC0YD____8',2e-4,...
        'AC0XQ'     ,1e-7,...
        'AC0ZQ'     ,1e-6,...
        'AC0YQ____1',1e-4,...
        'AC0YQ____2',1e-4,...
        'AC0YQ____3',1e-4,...
        'AC0YQ____4',1e-4,...
        'AC0YQ____5',1e-4,...
        'AC0YQ____6',1e-4,...
        'AC0YQ____7',1e-4,...
        'AC0YQ____8',1e-4...
      ) ...
    );
    plot_pars=struct(...
      'size',200+[0,0,16,9]*65,...
      'units','points',...
      'visible','on',...
      'fontsize',struct(...
        'axis', 24,...
        'title',32,...
        'label',28),...
      'line',struct(...
        'width',2)...
    );
    parameter_list=struct(...
      'calpar_csr_defs',struct('default',grace.default_list.calpar_csr_defs,'validation',@(i) isstruct(i)),...
      'start',        struct('default',datetime([0 0 31]),                'validation',@(i) isdatetime(i)),...
      'stop',         struct('default',datetime([0 0 31]),                'validation',@(i) isdatetime(i))...
    );
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={};
  end
  %read only
  properties(SetAccess=private)
    calpar_csr_defs
    start
    stop
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    calpar_csr_reg
    acc_csr_reg
  end
  %calculated only when asked for
  properties(Dependent)
    calpar_csr
    acc_csr
  end
  methods(Static)
    function out=parameters
      out=fieldnames(grace.parameter_list);
    end
    function enforce_plot_pars(par,fig_handle,axis_handle)
      set(axis_handle,'FontSize',par.fontsize.axis)
      set(get(axis_handle,'Title'),'FontSize',par.fontsize.title);
      set(get(axis_handle,'XLabel'),'FontSize',par.fontsize.label);
      set(get(axis_handle,'YLabel'),'FontSize',par.fontsize.label);
      set(fig_handle,'Position',par.size,'PaperUnits',par.units,'PaperPosition',par.size);
    end
  end
  methods
    %% constructor
    function obj=grace(varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('calpar_csr', @(i) ischar(i));
      p.addParameter('acc_csr',    @(i) ischar(i));
      %declare parameters
      for j=1:numel(grace.parameters)
        %shorter names
        pn=grace.parameters{j};
        %declare parameters
        p.addParameter(pn,grace.parameter_list.(pn).default,grace.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(varargin{:});
      % save parameters with defaults first, they may be needed below
      for i=1:numel(grace.parameters)
        %shorter names
        pn=grace.parameters{i};
        if ~isscalar(p.Results.(pn)) && ~iscell(p.Results.(pn))
          %vectors are always lines (easier to handle strings)
          obj.(pn)=transpose(p.Results.(pn)(:));
        else
          obj.(pn)=p.Results.(pn);
        end
      end
      % gather input argument names
      in=fieldnames(p.Results);
      for i=1:numel(in)
        switch lower(in{i})
          case 'calpar_csr'
            obj.calpar_csr=p.Results.(in{i});
          case 'acc_csr'
            obj.acc_csr=p.Results.(in{i});
          case grace.parameters{:}
            % already done
          otherwise
            error([mfilename,': cannot handle parameter ''',in{i},'''.'])
        end
      end
    end
    %% calpar_csr
    function obj=set.calpar_csr(obj,storage)
      jobID_idx=2;
      datafile=fullfile(storage,'calpar_csr.mat');
      if isempty(dir(datafile))
        calpars=fieldnames(obj.calpar_csr_defs);
        obj.calpar_csr_reg.a=cell(numel(calpars),1);
        obj.calpar_csr_reg.b=cell(numel(calpars),1);
        for j=1:size(calpars)
          for sat={'A','B'}
            f=fullfile(storage,['GRC-',sat{1},'____',calpars{j},'.GraceAccCal']);
            tmp=simpletimeseries.import(f,true);
            obj.calpar_csr_reg.(lower(sat{1})){j}=tmp.resample(days(1));
          end
        end
        %ensure the time domain is the same in all calpars (in each sat)
        obj.calpar_csr_reg.a=simpledata.merge_multiple(obj.calpar_csr_reg.a,'GRC-A calpar_csr');
        obj.calpar_csr_reg.b=simpledata.merge_multiple(obj.calpar_csr_reg.b,'GRC-B calpar_csr');
        %ensure job IDs are consistent
        equal_idx=simpledata.isequal_multiple(obj.calpar_csr_reg.a,jobID_idx,'GRC-A calpar_csr');
        if any(~equal_idx)
          error([mfilename,': Job ID inconsistency between ',calpars{find(~equal_idx,1,'first')},...
            ' and ',calpars{find(~equal_idx,1,'first')+1},' of GRC-A calpar_csr (possibly more).'])
        end
        equal_idx=simpledata.isequal_multiple(obj.calpar_csr_reg.b,jobID_idx,'GRC-B calpar_csr');
        if any(~equal_idx)
          error([mfilename,': Job ID inconsistency between ',calpars{find(~equal_idx,1,'first')},...
            ' and ',calpars{find(~equal_idx,1,'first')+1},' of GRC-B calpar_csr (possibly more).'])
        end
        %save data
        save(datafile,'obj.calpar_csr_reg');
      else
        %load data
        obj.calpar_csr_reg=load(datafile,'obj.calpar_csr_reg');
      end
    end
    function out=get.calpar_csr(obj)
      out=obj.calpar_csr_reg;
    end
    %% acc_csr
    function obj=set.acc_csr(obj,storage)
      RL='61';
      datelist=simpletimeseries.list(obj.start,obj.stop,days(1));
      for sat={'A','B'}
        %build required file list
        filelist=cell(size(datelist));
        for i=1:numel(datelist)
          filelist{i}=fullfile(storage,['ACC1B_',datestr(datelist{i},'yyyy-mm-dd'),'_',sat{1},'_',RL,'2.asc']);
        end
        %load the data
        obj.acc_csr.(lower(sat{1}))=simpletimeseries.import(filelist);
      end
    end
    function out=get.acc_csr(obj)
      out=obj.acc_csr_reg;
    end
    %% operate dataype-wise
    function obj=op(obj,operation,datatype,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'operation',                   @(i) ischar(i));
      p.addRequired( 'datatype',                    @(i) ischar(i));
      % parse it
      p.parse(operation,datatype);
      %gather datatype names
      dt_names=fieldnames(obj.(datatype));
      %loop over all fields of this datatype
      for i=1:numel(dt_names)
        obj.(datatype).(dt_names{i}).a=obj.(datatype).(dt_names{i}).a.(p.Results.operation)(varargin{:});
        obj.(datatype).(dt_names{i}).b=obj.(datatype).(dt_names{i}).b.(p.Results.operation)(varargin{:});
      end
    end
    %% costumized operations
    function out=calpar_csr_op(obj,opname,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('plot_dir',    '.',@(i) ischar(i));
      p.addParameter('plot_prefix', '', @(i) ischar(i));
      p.addParameter('plot_suffix', '', @(i) ischar(i));
      p.addParameter('plot_xlimits',[], @(i) isempty(i) || isdatetime(i));
      p.addParameter('plot_pars',  grace.plot_pars, @(i) isstruct(i));
      % parse it
      p.parse(varargin{:});
      switch lower(opname)
      case 'stats-modes'
        out={'mean','std','rms','length'};
      case 'stats2-modes'
        out={'corrcoef','length'};
      case 'column-idx'
        out=1;
      case 'load-stats'
        % case-specific parameters
        stats_modes=obj.calpar_csr_op('stats-modes');
        % loop over CSR cal pars
        calpars=fieldnames(obj.calpar_csr_defs);
        stats_filename=fullfile(p.Results.plot_dir,'stats.mat');
        if isempty(dir(stats_filename))
          s=cell(size(calpars));
          for i=1:size(calpars)
            for sat={'a','b'};
              if i>1 && ~obj.calpar_csr.(sat{1}){1}.istequal(obj.calpar_csr.(sat{1}){i})
                error([mfilename,': time domain discrepancy in ',calpars{i}])
              end
              %compute stats per period
              s{i}.(sat{1})=obj.calpar_csr.(sat{1}){i}.stats(...
                'period',years(1)/12,...
                'overlap',seconds(0),...
                'struct_fields',stats_modes...
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
        % case-specific parameters
        stats_modes=obj.calpar_csr_op('stats2-modes');
        % loop over CSR cal pars
        calpars=fieldnames(obj.calpar_csr_defs);
        stats_filename=fullfile(p.Results.plot_dir,'stats2.mat');
        if isempty(dir(stats_filename))
          s=cell(size(calpars));
          for i=1:size(calpars)
            %compute stats per period
            s{i}=simpletimeseries.stats2(...
              obj.calpar_csr.a{i},...
              obj.calpar_csr.b{i},...
              'period',years(1)/12,...
              'overlap',seconds(0),...
              'struct_fields',stats_modes...
            );
          end
          save(stats_filename,'s')
        else
          load(stats_filename)
        end
        %outputs
        out=s;
      end
    end
    %% costumized plotting routine
    function plot(obj,datatype,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('plot_dir',    '.',@(i) ischar(i));
      p.addParameter('plot_prefix', '', @(i) ischar(i));
      p.addParameter('plot_suffix', '', @(i) ischar(i));
      p.addParameter('plot_xlimits',[], @(i) isempty(i) || isdatetime(i));
      p.addParameter('plot_pars',  grace.plot_pars, @(i) isstruct(i));
      % parse it
      p.parse(varargin{:});
      switch lower(datatype)
      case 'calpar_csr'
        % loop over CSR cal pars
        calpars=fieldnames(obj.calpar_csr_defs);
        for j=1:size(calpars)
          calpar_clean=str.clean(calpars{j},'_');
          filename=fullfile(p.Results.plot_dir,...
            strrep([p.Results.plot_prefix,calpar_clean,p.Results.plot_suffix,'.png'],' ','-')...
          );
          if isempty(dir(filename))
            figure('visible',p.Results.plot_pars.visible);
            %make sure data is compatible
            obj.calpar_csr.a{j}.compatible(obj.calpar_csr.b{j})
            %plot data
            obj.calpar_csr.a{j}.plot(varargin{:},'LineWidth',p.Results.plot_pars.line.width)
            obj.calpar_csr.b{j}.plot(varargin{:},'LineWidth',p.Results.plot_pars.line.width)
            %legend
            legend({...
              obj.calpar_csr.a{j}.labels{1},...
              obj.calpar_csr.b{j}.labels{1}...
            });
            %fix axis (if needed/requested)
            scale=obj.calpar_csr_defs.(calpars{j});
            v=axis;
            if ~isempty(p.Results.plot_xlimits)
              v(1:2)=datenum(p.Results.plot_xlimits);
            end
            if any(abs(v(3:4))>scale)
              v(3:4)=[-scale,scale];
            end
            axis(v)
            grid on
            ylabel(obj.calpar_csr.a{j}.y_units{1}) %could also be obj.calpar_csr.b{j}.y_units{1}
            title(['GRACE Acc Cal ',calpar_clean])
            grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
            saveas(gcf,filename)
          end
        end
      case 'calpar_csr-count'
        for sat={'a','b'}
          %some shortcuts
          obj1=obj.calpar_csr.(sat{1}){1};
          n=numel(obj.calpar_csr.(sat{1}));
          %collect masks from all calpars
          m=struct(...
            'x',true(obj1.length,3), 'x_idx',[1,   11,   21   ],...
            'y',true(obj1.length,24),'y_idx',[3:10,13:20,23:30],...
            'z',true(obj1.length,3), 'z_idx',[2,   12,   22   ]...
          );
          for i=1:n
            if i>1 && ~obj1.istequal(obj.calpar_csr.(sat{1}){i})
              error([mfilename,': time domain discrepancy'])
            end
            if any(i==m.x_idx)
              m.x(:,(m.x_idx==i))=obj.calpar_csr.(sat{1}){i}.mask;
            end
            if any(i==m.y_idx)
              m.y(:,(m.y_idx==i))=obj.calpar_csr.(sat{1}){i}.mask;
            end
            if any(i==m.z_idx)
              m.z(:,(m.z_idx==i))=obj.calpar_csr.(sat{1}){i}.mask;
            end
          end
          filename=fullfile(p.Results.plot_dir,strrep(...
            [p.Results.plot_prefix,lower(datatype),'.GRACE-',upper(sat{1}),p.Results.plot_suffix,'.png'],...
          ' ','-'));
          if isempty(dir(filename))
            figure('visible',p.Results.plot_pars.visible);
            plot(datenum(obj1.t),sum(m.x,2),'x'), hold on
            plot(datenum(obj1.t),sum(m.y,2),'o')
            plot(datenum(obj1.t),sum(m.z,2),'+')
            %fix axis (if requested)
            v=axis;
            if ~isempty(p.Results.plot_xlimits)
              v(1:2)=datenum(p.Results.plot_xlimits);
            end
            axis(v)
            legend({'x','y','z'})
            datetick('x','yyyy-mm')
            xlabel('time')
            grid on
            ylabel('nr cal pars')
            title(['GRACE ',upper(sat{1})])
            grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
            saveas(gcf,filename)
          end
        end
      case 'calpar_csr-stats'
        s=obj.calpar_csr_op('load-stats',varargin{:});
        calpars=fieldnames(obj.calpar_csr_defs);
        column_idx=obj.calpar_csr_op('column-idx');
%         %plot number of data
%         mode='length';
%         for i=1:numel(calpars)
%           calpar_clean=str.clean(calpars{i},'_');
%           filename=fullfile(p.Results.plot_dir,...
%             strrep([p.Results.plot_prefix,mode,'.',calpar_clean,p.Results.plot_suffix,'.png'],' ','-')...
%           );
%           if isempty(dir(filename))
%             figure('visible',p.Results.plot_pars.visible);
%             b=bar(mean([datenum(s{i}.a.t),datenum(s{i}.b.t)],2),[s{i}.a.length(:,1),s{i}.b.length(:,1)]);
%             b(1).EdgeColor = 'red';
%             b(1).FaceColor = 'red';
%             b(2).EdgeColor = 'black';
%             b(2).FaceColor = 'black';
%             %fix axis (if requested)
%             v=axis;
%             if ~isempty(p.Results.plot_xlimits)
%               v(1:2)=datenum(p.Results.plot_xlimits);
%             end
%             axis(v)
%             legend('GRACE-A','GRACE-B')
%             datetick('x','yyyy-mm')
%             xlabel('time')
%             grid on
%             title(['nr cal pars/month for ',calpar_clean])
%             grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
%             saveas(gcf,filename)
%           end
%         end
        %plot more stats
        for mode=obj.calpar_csr_op('stats-modes');
          for i=1:numel(calpars)
            calpar_clean=str.clean(calpars{i},'_');
            filename=fullfile(p.Results.plot_dir,...
              strrep([p.Results.plot_prefix,mode{1},'.',calpar_clean,p.Results.plot_suffix,'.png'],' ','-')...
            );
            if isempty(dir(filename))
              switch mode{1}
              case 'length'
                mask=true(s{i}.a.(mode{1}).length,1);
              otherwise
                mask=s{i}.a.length.y(:,1)>=10;
              end    
              figure('visible',p.Results.plot_pars.visible);
              h=s{i}.a.(mode{1}).mask_and(mask).mask_update.plot('columns',column_idx);
              set(h.handle{1},'LineWidth',p.Results.plot_pars.line.width);
              h=s{i}.b.(mode{1}).mask_and(mask).mask_update.plot('columns',column_idx);
              set(h.handle{1},'LineWidth',p.Results.plot_pars.line.width);
              %fix axis (if needed/requested)
              scale=obj.calpar_csr_defs.(calpars{i});
              v=axis;
              if ~isempty(p.Results.plot_xlimits)
                v(1:2)=datenum(p.Results.plot_xlimits);
              end
              if any(abs(v(3:4))>scale)
                switch mode{1}
                case 'mean'
                  v(3:4)=[-scale,scale];
                case {'std','rms'}
                  v(3:4)=[0,scale];
                end
              end
              axis(v)
              legend('GRACE-A','GRACE-B')
              ylabel(['[',s{i}.b.(mode{1}).y_units{column_idx},']'])
              grid on
              title([mode{1},' for ',calpar_clean])
              grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
              saveas(gcf,filename)
            end
          end
        end
      case 'calpar_csr-corr'
        s=obj.calpar_csr_op('load-stats2',varargin{:});
        column_idx=obj.calpar_csr_op('column-idx');
        % loop over all required modes
        for mode=obj.calpar_csr_op('stats2-modes',varargin{:});
          calpars=fieldnames(obj.calpar_csr_defs);
          % loop over CSR cal pars
          for i=1:size(calpars)
            calpar_clean=str.clean(calpars{i},'_');
            filename=fullfile(p.Results.plot_dir,...
              strrep([p.Results.plot_prefix,mode{1},'.',calpar_clean,p.Results.plot_suffix,'.png'],' ','-')...
            );
            if isempty(dir(filename))
              switch mode{1}
              case 'length'
                mask=true(s{i}.(mode{1}).length,1);
              otherwise
                mask=s{i}.length.y(:,1)>=10;
              end    
              figure('visible',p.Results.plot_pars.visible);
              h=s{i}.(mode{1}).mask_and(mask).mask_update.plot('columns',column_idx);
              set(h.handle{1},'LineWidth',p.Results.plot_pars.line.width);
              %fix axis (if needed/requested)
              v=axis;
              if ~isempty(p.Results.plot_xlimits)
                v(1:2)=datenum(p.Results.plot_xlimits);
              end
              axis(v)
              datetick('x','yyyy-mm')
              xlabel('time')
              grid on
              title([mode{1},' for ',calpar_clean])
              grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
              saveas(gcf,filename)
            end
          end
        end
      case 'calpar_csr-tables'
        calpars=fieldnames(obj.calpar_csr_defs);
        column_idx=obj.calpar_csr_op('column-idx');
        tab=12;
        %single-sat stats
        s=obj.calpar_csr_op('load-stats',varargin{:});
        stats_modes=obj.calpar_csr_op('stats-modes');
        for sat={'a','b'};
          disp(['GRACE-',upper(sat{1})])
          for i=0:numel(calpars)
            out=cell(1,numel(stats_modes)+1);
            if i==0
              calpar_clean='cal par';
            else
              calpar_clean=str.clean(calpars{i},'_');
            end
            out{1}=str.tabbed(calpar_clean,tab);
            for j=1:numel(stats_modes);
              if i==0
                out{j+1}=str.tabbed(stats_modes{j},tab,true);
              else
                d=s{i}.(sat{1}).(stats_modes{j}).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
                out{j+1}=str.tabbed(num2str(d(column_idx),'% .3g'),tab,true);
              end
            end
            disp(strjoin(out,' '))
          end
        end
        %all-sat stats
        s=obj.calpar_csr_op('load-stats2',varargin{:});
        stats_modes=obj.calpar_csr_op('stats2-modes');
        for i=0:numel(calpars)
          out=cell(1,numel(stats_modes)+1);
          if i==0
            calpar_clean='cal par';
          else
            calpar_clean=str.clean(calpars{i},'_');
          end
          out{1}=str.tabbed(calpar_clean,tab);
          for j=1:numel(stats_modes);
            if i==0
              out{j+1}=str.tabbed(stats_modes{j},tab,true);
            else
              d=s{i}.(stats_modes{j})(:,1).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
              out{j+1}=str.tabbed(num2str(mean(d(column_idx)),'% .3f'),tab,true);
            end
          end
          disp(strjoin(out,' '))
        end
      case 'calpar_csr-overview'
        calpars=fieldnames(obj.calpar_csr_defs);
        column_idx=obj.calpar_csr_op('column-idx');
        tab=12;
        %single-sat stats
        s=obj.calpar_csr_op('load-stats',varargin{:});
        stats_modes=obj.calpar_csr_op('stats-modes');
        for j=1:numel(stats_modes);
          filename=fullfile(p.Results.plot_dir,...
            strrep([p.Results.plot_prefix,'overview.',stats_modes{j},p.Results.plot_suffix,'.png'],' ','-')...
          );
          if isempty(dir(filename))
            figure('visible',p.Results.plot_pars.visible);
            %gathering data
            y=ones(numel(calpars),2);
            for i=1:numel(calpars)
              y(i,1)=s{i}.a.(stats_modes{j}).trim(...
                p.Results.plot_xlimits(1),...
                p.Results.plot_xlimits(2)...
              ).stats(...
                'period',days(inf),...
                'outlier',true,...
                'nsigma',3 ...
              ).mean(column_idx);
              y(i,2)=s{i}.b.(stats_modes{j}).trim(...
                p.Results.plot_xlimits(1),...
                p.Results.plot_xlimits(2)...
              ).stats(...
                'period',days(inf),...
                'outlier',true,...
                'nsigma',3 ...
              ).mean(column_idx);
            end
            pos_idx=(y>0);
            %plot negative parameters downwards
            y_now=abs(y); y_now(pos_idx)=nan;
            b=bar(log10(y_now)); hold on
            b(1).EdgeColor = 'red';
            b(1).FaceColor = 'red';
            b(2).EdgeColor = 'black';
            b(2).FaceColor = 'black';
            %plot positive parameters upwards
            y_now=abs(y); y_now(~pos_idx)=nan;
            b=bar(-log10(y_now)); hold on
            b(1).EdgeColor = 'red';
            b(1).FaceColor = 'red';
            b(2).EdgeColor = 'black';
            b(2).FaceColor = 'black';

            %getting parameter names
            x=cell(1,numel(calpars));
            for i=1:numel(calpars)
              x{i}=str.clean(calpars{i},'_');
            end
            legend('GRACE-A','GRACE-B')
            set(gca,'XTick',1:numel(x))
            set(gca,'XTickLabel',x)
            set(gca,'XTickLabelRotation',45)
            grid on
            title([...
              stats_modes{j},' for ',...
              datestr(p.Results.plot_xlimits(1)),' to ',...
              datestr(p.Results.plot_xlimits(2))...
            ]);
            axis([0 numel(calpars)+1 -10 10])
            grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
            saveas(gcf,filename)
          end
        end
        %all-sat stats
        s=obj.calpar_csr_op('load-stats2',varargin{:});
        stats_modes=obj.calpar_csr_op('stats2-modes');
        for j=1:numel(stats_modes);
          filename=fullfile(p.Results.plot_dir,...
            strrep([p.Results.plot_prefix,'overview.',stats_modes{j},p.Results.plot_suffix,'.png'],' ','-')...
          );
          if isempty(dir(filename))
            figure('visible',p.Results.plot_pars.visible);
            %gathering data
            y=ones(numel(calpars),1);
            for i=1:numel(calpars)
              y(i,1)=s{i}.(stats_modes{j}).trim(...
                p.Results.plot_xlimits(1),...
                p.Results.plot_xlimits(2)...
              ).stats(...
                'period',days(inf),...
                'outlier',true,...
                'nsigma',3 ...
              ).mean(column_idx);
            end
            bar(y)
            %getting parameter names
            x=cell(1,numel(calpars));
            for i=1:numel(calpars)
              x{i}=str.clean(calpars{i},'_');
            end
            set(gca,'XTick',1:numel(x))
            set(gca,'XTickLabel',x)
            set(gca,'XTickLabelRotation',45)
            grid on
            title([...
              stats_modes{j},' for ',...
              datestr(p.Results.plot_xlimits(1)),' to ',...
              datestr(p.Results.plot_xlimits(2))...
            ]);
            axis([0 numel(calpars)+1 -1 1])
            grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
            saveas(gcf,filename)
          end
        end        
      otherwise
        error([mfilename,': unknown data of type ''',datatype,'''.'])
      end
    end
  end
end
