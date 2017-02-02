classdef csr
  methods(Static)
    function obj=import_calpar(obj,dataname,varargin)
      %retrieve product info
      product=obj.mdget(dataname);
      %check if data is already in matlab format
      if ~product.isfile('data')
        %get names of parameters and levels
        levels    =product.mdget('levels');
        fields    =product.mdget('fields');
        sats      =product.mdget('sats');
        bias_files=product.mdget('bias_files');
        %need to get long-term biases
        for s=1:numel(sats)
          ltb.(sats{s})=flipud(transpose(dlmread(bias_files{s})));
        end
        %load data
        for i=1:numel(levels)
          for j=1:numel(fields)
            tmp=struct('A',[],'B',[]);
            %read data
            for s=1:numel(sats)
              f=fullfile(product.mdget('import_dir'),['gr',sats{s},'.',fields{j},'.',levels{i},'.GraceAccCal']);
              tmp.(sats{s})=simpletimeseries.import(f,'cut24hrs',false);
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
                t=tmp.(sats{s}).mjd-ltb.(sats{s})(2,1);
                tmp.(sats{s})=tmp.(sats{s}).assign(...
                  [tmp.(sats{s}).y(:,1)+polyval(ltb.(sats{s})(:,lbt_idx),t),tmp.(sats{s}).y(:,2:end)]...
                );
              end
            end
            %ensure date is compatible between the satellites
            if ~tmp.A.isteq(tmp.B)
              [tmp.A,tmp.B]=tmp.A.merge(tmp.B);
            end
            %propagate data to object
            for s=1:numel(sats)
              obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},tmp.(sats{s}));
            end
            disp(['loaded data for ',str.just(levels{i},6),' and ',str.just(fields{j},6)])
          end
        end
        %loop over all sat and level to check for consistent time domain and Job IDs agreement
        for i=1:numel(levels)
          for s=1:numel(sats)
            %gather names for this sat and level
            names=obj.vector_names(product.dataname.type,levels{i},'',sats{s});
            %ensure the time domain is the same for all fields (in each sat and level)
            obj=obj.vector_op(@simpledata.merge_multiple,...
              names,...
              ['GRACE-',sats{s},' ',str.just(levels{i},6),' ',product.dataname.type]...
            );
            %ensure job IDs are consistent for all fields (in each sat and level)
            equal_idx=obj.vector_tr(@simpledata.isequal_multiple,...
              names,...
              product.mdget('jobid_col'),['GRACE-',sats{s},' ',str.just(levels{i},6),' ',product.dataname.type]...
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
        %loop over all sats, levels and fields to:
        % - in case of estim: ensure that there are no arcs with lenghts longer than consecutive time stamps
        % - in case of aak and accatt: ensure that the t0 value is the same as the start of the arc
        for s=1:numel(sats)
          %loop over all required levels
          for i=1:numel(levels)
            switch levels{i}
            case 'estim'
              %this check ensures that there are no arcs with lenghts longer than consecutive time stamps
              for j=1:numel(fields)
                disp(['Checking ',product.dataname.type,' ',str.just(levels{i},6),' ',str.just(fields{j},6),' ',sats{s}])
                %save time series into dedicated var
                ts_now=obj.sat_get(product.dataname.type,levels{i},fields{j},sats{s});
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
                  msg=cell(1,min([numel(bad_idx),10])+1);
                  msg{1}='idx: arc init time; arc length; succ time diff; delta arc len (should be zero)';
                  for k=1:numel(msg)-1
                    idx=idx1(bad_idx(k));
                    msg{k+1}=[...
                      num2str(idx1(bad_idx(k))),': ',...
                      datestr(ts_now.t(idx)),'; ',...
                      num2str(al(bad_idx(k)),'%.5d'),'; ',...
                      num2str(dt(bad_idx(k))),' ',...
                      num2str(al(bad_idx(k))-dt(bad_idx(k)))...
                    ];
                  end
                  disp([....
                    ': found ',num2str(numel(bad_idx)),' arc lengths (3rd column) longer than ',...
                    ' difference between consecutive time stamps (4th column):',10,...
                    strjoin(msg,'\n'),10,...
                    'These data have been discarded!'
                  ])
                  mask=ts_now.mask;
                  mask(idx1(bad_idx))=false;
                  obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},ts_now.mask_and(mask).mask_update);
                end
              end
            case {'aak','accatt'}
              %this check ensures that the t0 value is the same as the start of the arc
              for j=1:numel(fields)
                %some fields do not have t0
                if ~any(fields{j}(end)=='DQ')
                  disp(['Skipping ',product.dataname.name,' ',str.just(levels{i},6),' ',str.just(fields{j},6),' ',sats{s}])
                  continue
                end
                disp(['Checking ',product.dataname.name,' ',str.just(levels{i},6),' ',str.just(fields{j},6),' ',sats{s}])
                %save time series into dedicated var
                ts_now=obj.sat_get(product.dataname.type,levels{i},fields{j},sats{s});
                %forget about epochs that have been artificially inserted to represent forward steps
                idx1=find(diff(ts_now.t)>seconds(1));
                %get t0
                t0=simpletimeseries.utc2gps(datetime(ts_now.y(idx1,3),'convertfrom','modifiedjuliandate'));
                %find arcs that have (much) t0 different than their first epoch
                bad_idx=find(abs(ts_now.t(idx1)-t0)>seconds(1) & ts_now.mask(idx1));
                %report if any such epochs have been found
                if ~isempty(bad_idx)
                  msg=cell(1,min([numel(bad_idx),10])+1);
                  msg{1}='idx: arc init time - MJD = delta time (should be zero)';
                  for k=1:numel(msg)-1
                    idx=idx1(bad_idx(k));
                    msg{k+1}=[...
                      num2str(idx1(bad_idx(k))),': ',...
                      datestr(ts_now.t(idx1(bad_idx(k))),'yyyy-mm-dd HH:MM:SS'),' - ',...
                      datestr(t0(bad_idx(k)),'yyyy-mm-dd HH:MM:SS'),' = ',...
                      char(ts_now.t(idx1(bad_idx(k)))-t0(bad_idx(k)))...
                    ];
                  end
                  disp([...
                    'found ',num2str(numel(bad_idx)),' arc init time (2nd column) different than the',...
                    ' MJD reported in the data (3rd column):',10,...
                    strjoin(msg,'\n'),10,...
                    'These data have been discarded!'
                  ])
                  mask=ts_now.mask;
                  mask(idx1(bad_idx))=false;
                  obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},ts_now.mask_and(mask).mask_update);
                end
              end
            end
          end
        end
        %save data
        s=obj.datatype_get(product.dataname.type); %#ok<*NASGU>
        save(char(product.file('data')),'s');
        clear s
      else
        %load data
        load(char(product.file('data')),'s');
        levels=fieldnames(s); %#ok<NODEF>
        for i=1:numel(levels)
          obj=obj.level_set(product.dataname.type,levels{i},s.(levels{i}));
        end
      end
    end
    function obj=compute_calmod(obj,dataname,varargin)
      %retrieve products info
      product=obj.mdget(dataname);
      %paranoid sanity
      if product.nr_sources~=2
        error([mfilename,': number of sources in product ',dataname,...
          ' is expected to be 2, not ',num2str(product.nr_sources),'.'])
      end
      %get sources
      calparp=obj.mdget(product.sources(1));
      l1baccp=obj.mdget(product.sources(2));
      %retrieve relevant parameters
      levels    =calparp.mdget('levels');
      sats      =calparp.mdget('sats');
      param_col =calparp.mdget('param_col');
      coords     =calparp.mdget('coords');
      for s=1:numel(sats)
        %gather quantities
        acc=obj.sat_get(l1baccp.dataname.type,l1baccp.dataname.level,l1baccp.dataname.field,sats{s});
        if ~isa(acc,'simpletimeseries')
          %patch nan calibration model
          calmod=simpletimeseries(...
            [obj.start;obj.stop],...
            nan(2,numel(obj.par.acc.data_col_name))...
          );
        else
          %loop over all 
          for l=1:numel(levels)
            %init models container
            calmod=simpletimeseries(acc.t,zeros(acc.length,numel(coords))).copy_metadata(acc);
            calmod.descriptor=['calibration model ',levels{l},' GRACE-',upper(sats{s})];
            disp(['Computing the ',calmod.descriptor])
            for i=1:numel(coords)
              %skip Y coordinate for now
              if coords{i}=='Y',continue,end
              %build nice structure with the relevant calibration parameters
              cal=struct(...
                'ac0' ,obj.sat_get(calparp.dataname.type,levels{l},['AC0',coords{i}    ],sats{s}).interp(acc.t),...
                'ac0d',obj.sat_get(calparp.dataname.type,levels{l},['AC0',coords{i},'D'],sats{s}).interp(acc.t),...
                'ac0q',obj.sat_get(calparp.dataname.type,levels{l},['AC0',coords{i},'Q'],sats{s}).interp(acc.t)...
              );
              %sanity
              if any([isempty(acc),isempty(cal.ac0),isempty(cal.ac0d),isempty(cal.ac0q)])
                error([mfilename,': not all data is available to perform this operation.'])
              end
              %retrieve time domain (it is the same for all cal pars)
              fields=fieldnames(cal);
              for f=1:numel(fields)
                t.(fields{f})=days(acc.t-simpletimeseries.ToDateTime(cal.(fields{f}).y(:,end),'modifiedjuliandate'));
              end
              %paranoid sanity check
              good_idx=~isnan(t.ac0);
              if any(t.ac0(good_idx)~=t.ac0d(good_idx)) || any(t.ac0(good_idx)~=t.ac0q(good_idx))
                error([mfilename,': calibration time domain inconsistent between parameters, debug needed!'])
              end
              %build calibration model
              calmod=calmod.set_cols(i,...
                cal.ac0.cols( param_col)+...
                cal.ac0d.cols(param_col).times(t.ac0d)+...
                cal.ac0q.cols(param_col).times(t.ac0q.^2)...
              );
            end
          end
        end
        %propagate it
        obj=obj.sat_set(dataname.type,dataname.level,levels{l},sats{s},calmod);
      end
    end
    function obj=import_acc_l1b(obj,dataname,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('dataname',@(i) ischar(i) || isa(i,'datanames'));
      p.addParameter('start', [], @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  [], @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(dataname,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
      %retrieve product info
      product=obj.mdget(dataname);
      %retrieve relevant parameters
      sats =product.mdget('sats');
      indir=product.mdget('import_dir');
      version=product.mdget('version');
      %gather list of daily data files
      [~,timestamplist]=product.file('data',varargin{:});
      %loop over the satellites
      for s=1:numel(sats)
        infile=cell(size(timestamplist));
        %loop over all dates
        for i=1:numel(timestamplist)
          %build input data filename
          infile{i}=fullfile(indir,datestr(timestamplist(i),'yy'),'acc','asc',...
            ['ACC1B_',datestr(timestamplist(i),'yyyy-mm-dd'),'_',sats{s},'_',version,'.asc']...
          );
        end
        %load (and save the data in mat format, as handled by simpletimeseries.import)
        obj=obj.sat_set(dataname.type,dataname.level,dataname.field,sats{s},...
          simpletimeseries.import(infile,'cut24hrs',false)...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
    function obj=import_acc_mod(obj,dataname,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('dataname',@(i) ischar(i) || isa(i,'datanames'));
      p.addParameter('start', [], @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  [], @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(dataname,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
      %retrieve product info
      product=obj.mdget(dataname);
      %retrieve relevant parameters
      sats =product.mdget('sats');
      indir=product.mdget('import_dir');
      acc_version =product.mdget('acc_version' );
      gps_version =product.mdget('gps_version' );
      grav_version=product.mdget('grav_version');
      %gather list of daily data files
      [~,timestamplist]=product.file('data',varargin{:});
      %loop over the satellites
      for s=1:numel(sats)
        infile=cell(size(timestamplist));
        %loop over all dates
        for i=1:numel(timestamplist)
          %build input data filename
          infile{i}=fullfile(indir,datestr(timestamplist(i),'yy'),datestr(timestamplist(i),'mm'),'gps_orb_l',...
            ['grc',sats{s},'_gps_orb_',datestr(timestamplist(i),'yyyy-mm-dd'),...
            '_RL',acc_version,'_GPSRL',gps_version,'_RL',grav_version,'.*.acc']...
          );
        end
        %load (and save the data in mat format, as handled by simpletimeseries.import)
        obj=obj.sat_set(dataname.type,dataname.level,dataname.field,sats{s},...
          simpletimeseries.import(infile,'cut24hrs',true)...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
  end
end