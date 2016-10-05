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
      'calpar_csr_defs',   struct('default',grace.default_list.calpar_csr_defs,     'validation',@(i) isstruct(i))...
    );
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={};
  end
  %read only
  properties(SetAccess=private)
    calpar_csr_defs
    calpar_csr
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    
  end
  %calculated only when asked for
  properties(Dependent)
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
            storage=p.Results.(in{i});
            calpars=fieldnames(obj.calpar_csr_defs);
            jobID_idx=2;
            obj.calpar_csr.a=cell(numel(calpars),1);
            obj.calpar_csr.b=cell(numel(calpars),1);
            for j=1:size(calpars)
              for sat={'A','B'}
                f=fullfile(storage,['GRC-',sat{1},'____',calpars{j},'.GraceAccCal']);
                tmp=simpletimeseries.import(f);
                obj.calpar_csr.(lower(sat{1})){j}=tmp.resample(days(1));
              end
            end
            %ensure the time domain is the same in all calpars (in each sat)
            obj.calpar_csr.a=simpledata.merge_multiple(obj.calpar_csr.a,['GRC-A ',in{i}]);
            obj.calpar_csr.b=simpledata.merge_multiple(obj.calpar_csr.b,['GRC-B ',in{i}]);
            %ensure job IDs are consistent
            equal_idx=simpledata.isequal_multiple(obj.calpar_csr.a,jobID_idx,['GRC-A ',in{i}]);
            if any(~equal_idx)
              error([mfilename,': Job ID inconsistency between ',calpars{find(~equal_idx,1,'first')},...
                ' and ',calpars{find(~equal_idx,1,'first')+1},' of GRC-A ',in{i},' (possibly more).'])
            end
            equal_idx=simpledata.isequal_multiple(obj.calpar_csr.b,jobID_idx,['GRC-B ',in{i}]);
            if any(~equal_idx)
              error([mfilename,': Job ID inconsistency between ',calpars{find(~equal_idx,1,'first')},...
                ' and ',calpars{find(~equal_idx,1,'first')+1},' of GRC-B ',in{i},' (possibly more).'])
            end
          case grace.parameters{:}
            % already done
          otherwise
            error([mfilename,': cannot handle parameter ''',in{i},'''.'])
        end
      end
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
    %% costumized plotting routine
    function plot(obj,datatype,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('plot_dir',   '.',@(i) ischar(i));
      p.addParameter('plot_prefix','', @(i) ischar(i));
      p.addParameter('plot_suffix','', @(i) ischar(i));
      p.addParameter('plot_pars',  grace.plot_pars, @(i) isstruct(i));
      p.addParameter('sat',        'a',@(i) ischar(i));
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
            %fix scale (if needed)
            scale=obj.calpar_csr_defs.(calpars{j});
            v=axis;
            if any(abs(v(3:4))>scale)
              axis([v(1) v(2) -scale scale])
            end
            grid on
            ylabel(obj.calpar_csr.a{j}.y_units{1}) %could also be obj.calpar_csr.b{j}.y_units{1}
            title(['GRACE Acc Cal ',calpar_clean])
            grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
            saveas(gcf,filename)
          end
        end
      case 'calpar_csr-count'
        %some shortcuts
        obj1=obj.calpar_csr.(p.Results.sat){1};
        n=numel(obj.calpar_csr.(p.Results.sat));
        %collect masks from all calpars
        m=struct(...
          'x',true(obj1.length,3),...
          'x_idx',[1,11,21],...
          'y',true(obj1.length,24),...
          'y_idx',[3:10,13:20,23:30],...
          'z',true(obj1.length,3),...
          'z_idx',[2,12,22]...
        );
        for i=1:n
          if i>1 && ~obj1.istequal(obj.calpar_csr.(p.Results.sat){i})
            error([mfilename,': time domain discrepancy'])
          end
          if any(i==m.x_idx)
            m.x(:,(m.x_idx==i))=obj.calpar_csr.(p.Results.sat){i}.mask;
          end
          if any(i==m.y_idx)
            m.y(:,(m.y_idx==i))=obj.calpar_csr.(p.Results.sat){i}.mask;
          end
          if any(i==m.z_idx)
            m.z(:,(m.z_idx==i))=obj.calpar_csr.(p.Results.sat){i}.mask;
          end
        end
        filename=fullfile(p.Results.plot_dir,strrep(...
          [p.Results.plot_prefix,lower(datatype),'.GRACE-',upper(p.Results.sat),p.Results.plot_suffix,'.png'],...
        ' ','-'));
        if isempty(dir(filename))
          figure('visible',p.Results.plot_pars.visible);
          plot(datenum(obj1.t),sum(m.x,2),'x'), hold on
          plot(datenum(obj1.t),sum(m.y,2),'o')
          plot(datenum(obj1.t),sum(m.z,2),'+')
          legend({'x','y','z'})
          datetick('x','yyyy-mm')
          xlabel('time')
          grid on
          ylabel('nr cal pars')
          title(['GRACE ',upper(p.Results.sat)])
          grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
          saveas(gcf,filename)
        end
      case 'calpar_csr-stats'
        % loop over CSR cal pars
        calpars=fieldnames(obj.calpar_csr_defs);
        stats_filename=fullfile(p.Results.plot_dir,'stats.mat');
        if isempty(dir(stats_filename))
          s=cell(size(calpars));
          for i=1:size(calpars)
            if i>1 && ~obj.calpar_csr.(p.Results.sat){1}.istequal(obj.calpar_csr.(p.Results.sat){i})
              error([mfilename,': time domain discrepancy'])
            end
            %compute stats per period
            s{i}.a=obj.calpar_csr.a{i}.stats(years(1)/12);
            s{i}.b=obj.calpar_csr.b{i}.stats(years(1)/12);
          end
          save(stats_filename,'s')
        else
          load(stats_filename)
        end
        %plot number of data
        mode='length';
        for i=1:numel(calpars)
          calpar_clean=str.clean(calpars{i},'_');
          filename=fullfile(p.Results.plot_dir,...
            strrep([p.Results.plot_prefix,mode,'.',calpar_clean,p.Results.plot_suffix,'.png'],' ','-')...
          );
          if isempty(dir(filename))
            figure('visible',p.Results.plot_pars.visible);
            for sat={'a','b'};
              for j=1:numel(s{i}.(sat{1}))
                n=size(s{i}.(sat{1})(:));
                %build time and data domains
                t.(sat{1})=zeros(n);
                d.(sat{1})=zeros(n);
                for k=1:numel(t.(sat{1}))
                  t.(sat{1})(k)=datenum(s{i}.(sat{1}){k}.t);
                  if isempty(s{i}.(sat{1}){k}.(mode))
                    d.(sat{1})(k)=NaN;
                  else
                    d.(sat{1})(k)=s{i}.(sat{1}){k}.(mode);
                  end
                end
              end
            end
            b=bar(mean([t.a,t.b],2),[d.a,d.b]);
            b(1).EdgeColor = 'red';
            b(1).FaceColor = 'red';
            b(2).EdgeColor = 'black';
            b(2).FaceColor = 'black';
            v=axis;
            axis([datenum(datetime('2002-01-01')),datenum(datetime('2016-12-31')),v(3),v(4)])
            legend('GRACE-A','GRACE-B')
            datetick('x','yyyy-mm')
            xlabel('time')
            grid on
            title(['nr cal pars/month for ',calpar_clean])
            grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
            saveas(gcf,filename)
          end
        end
        %plot more stats
        for mode={'mean','std','rms'};
          for i=1:numel(calpars)
            calpar_clean=str.clean(calpars{i},'_');
            filename=fullfile(p.Results.plot_dir,...
              strrep([p.Results.plot_prefix,mode{1},'.',calpar_clean,p.Results.plot_suffix,'.png'],' ','-')...
            );
            if isempty(dir(filename))
              figure('visible',p.Results.plot_pars.visible);
              for sat={'a','b'};
                for j=1:numel(s{i}.(sat{1}))
                  n=size(s{i}.(sat{1})(:));
                  %build time and data domains
                  t.(sat{1})=zeros(n);
                  d.(sat{1})=zeros(n);
                  for k=1:numel(t.(sat{1}))
                    t.(sat{1})(k)=datenum(s{i}.(sat{1}){k}.t);
                    if isempty(s{i}.(sat{1}){k}.(mode{1})) || s{i}.(sat{1}){k}.length<5
                      d.(sat{1})(k)=NaN;
                    else
                      d.(sat{1})(k)=s{i}.(sat{1}){k}.(mode{1})(1);
                    end
                  end
                end
              end
              plot(t.a,d.a,t.b,d.b,,'LineWidth',p.Results.plot_pars.line.width), hold on
              %fix scale
              scale=obj.calpar_csr_defs.(calpars{i});
              v=axis;
              if any(abs(v(3:4))>scale)
                switch mode{1}
                case 'mean'
                  axis([datenum(datetime('2002-01-01')),datenum(datetime('2016-12-31')),-scale,scale])
                case {'std','rms'}
                  axis([datenum(datetime('2002-01-01')),datenum(datetime('2016-12-31')),0,scale])
                end
              else
                axis([datenum(datetime('2002-01-01')),datenum(datetime('2016-12-31')),v(3),v(4)])
              end
              legend('GRACE-A','GRACE-B')
              datetick('x','yyyy-mm')
              xlabel('time')
              grid on
              title([mode{1},' for ',calpar_clean])
              grace.enforce_plot_pars(p.Results.plot_pars,gcf,gca)
              saveas(gcf,filename)
            end
          end
        end
      otherwise
        error([mfilename,': unknown data of type ''',datatype,'''.'])
      end
    end
  end
end
