classdef csr
  methods(Static)
    function obj=import_calpar(obj,dataname,varargin)
      % parse mandatory arguments
      p=inputParser;
      p.addRequired('dataname', @(i) isa(i,'datanames'));
      p.parse(dataname);
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
              
              %additional processing: add end of arcs
              switch levels{i}
              case {'aak','accatt'}
                %get arc stars
                arc_starts=tmp.(sats{s}).t;
                %build arc ends
                arc_ends=[arc_starts(2:end);dateshift(arc_starts(end),'end','day')]-seconds(1);
                %arc ends are at maximum 24 hours after arc starts
                fix_idx=arc_ends-arc_starts>days(1);
                arc_ends(fix_idx)=arc_starts(fix_idx)+days(1);
              case 'estim'
                %get arc stars
                arc_starts=tmp.(sats{s}).t;
                %build arc ends
                arc_ends=arc_starts+seconds(tmp.(sats{s}).y(:,3))-seconds(1);
                %get seconds-of-day of arc ends
                sod_arc_ends=seconds(arc_ends-dateshift(arc_ends,'start','day'));
                %find the 24hrs arcs (those that have ~0 seconds of days)
                idx=find(sod_arc_ends<tmp.(sats{s}).t_tol);
                %push those arcs to midnight and remove 1 second
                arc_ends(idx)=dateshift(arc_ends(idx),'start','day')-seconds(1);
              end
              
              %surpress over-lapping arcs
              ov_idx=find(arc_starts(2:end)-arc_ends(1:end-1)<0);
              if ~isempty(ov_idx)
                arc_ends(ov_idx)=arc_starts(ov_idx+1)-seconds(1);
              end
              
              %get data and set the arc length to zero
              arc_end_y=tmp.(sats{s}).y;
              arc_end_y(:,3)=0;
              %build timeseries with arc ends
              arc_end_ts=simpletimeseries(arc_ends,arc_end_y,...
                'format','datetime',...
                'labels',tmp.(sats{s}).labels,...
                'units',tmp.(sats{s}).y_units,...
                'timesystem',tmp.(sats{s}).timesystem,...
                'descriptor',['end of arcs for ',tmp.(sats{s}).descriptor]...
              );
              
              %additional processing: add gaps
              gap_idx=[...
                abs(seconds(arc_ends(1:end-1)+seconds(1)-arc_starts(2:end))) > tmp.(sats{s}).t_tol;...
              false];
              gap_t=arc_ends(gap_idx)+seconds(1);
              %build timeseries with arc ends
              gap_ts=simpletimeseries(gap_t,nan(numel(gap_t),tmp.(sats{s}).width),...
                'format','datetime',...
                'labels',tmp.(sats{s}).labels,...
                'units',tmp.(sats{s}).y_units,...
                'timesystem',tmp.(sats{s}).timesystem,...
                'descriptor',['gaps for ',tmp.(sats{s}).descriptor]...
              );
            
              if obj.debug
%                 t0=datetime('16-Aug-2002');t1=datetime('21-Aug-2002');
%                 t0=datetime('04-Apr-2002');t1=datetime('08-Apr-2002');
%                 t0=datetime('12-Jan-2003');t1=datetime('16-Jan-2003');
%                 t0=datetime('2003-11-29');t1=datetime('2003-12-02');
                t0=datetime('2012-06-29');t1=datetime('2012-07-03');
                o=tmp.(sats{s}).trim(t0,t1);
                disp(str.tablify(22,'orignal t','original y'))
                for di=1:o.length
                  disp(str.tablify(22,o.t(di),o.y(di,1)))
                end
                 as=tmp.(sats{s}).trim(t0,t1);
                 ae=arc_end_ts.trim(t0,t1);
                 disp(str.tablify(22,'arc start t','arc start y','arc end t','arc end y'))
                 for di=1:min([as.length,ae.length])
                   disp(str.tablify(22,as.t(di),as.y(di,1),ae.t(di),ae.y(di,1)))
                 end
                 g=gap_ts.trim(t0,t1);
                 if ~isempty(g)
                   disp(str.tablify(22,'gap t','gap y'))
                   for di=1:g.length
                     disp(str.tablify(22,g.t(di),g.y(di,1)))
                   end
                 end
              end

              %augment the original timeseries with the end-of-arcs and gaps (only new data)
              tmp.(sats{s})=tmp.(sats{s}).augment(arc_end_ts,true).augment(gap_ts,true);
              
              
              if obj.debug
                au=tmp.(sats{s}).trim(t0,t1);
                disp(str.tablify(22,'augmented t','augmented y'))
                for di=1:au.length
                  disp(str.tablify(22,au.t(di),au.y(di,1)))
                end
                keyboard
              end

            end
            
%             %ensure date is compatible between the satellites
%             if ~tmp.A.isteq(tmp.B)
%               [tmp.A,tmp.B]=tmp.A.merge(tmp.B);
%             end
            %propagate data to object
            for s=1:numel(sats)
              obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},tmp.(sats{s}));
            end
            disp(str.tablify([15,6,3,6],'loaded data for',levels{i},'and',fields{j}))
          end
        end
%         %for each sat, ensure consistent time domain 
%         for s=1:numel(sats)
%           %gather names for this sat and level
%           names=obj.vector_names(product.dataname.type,'','',sats{s});
%           %gather info for user feedback
%           length_start=cellfun(@(i)(obj.data_get(i).length),names);
%           %ensure the time domain is the same for all fields (in each sat and level)
%           obj=obj.vector_op(@simpledata.merge_multiple,...
%             names,...
%             str.tablify([7,6,10],['GRACE-',sats{s}],levels{i},product.dataname.type)...
%           );
%           %user feedback
%           length_stop=cellfun(@(i)(obj.data_get(i).length),names);
%           cellfun(...
%             @(i,j,k)(disp(str.tablify([32,6,24,4],i.name,j,'data entries, changed by',k))),...
%             names,num2cell(length_start),num2cell(length_stop-length_start))
%         end
        %loop over all sat and level to check Job IDs agreement
        for i=1:numel(levels)
          for s=1:numel(sats)
            %gather names for this sat and level
            names=obj.vector_names(product.dataname.type,levels{i},'',sats{s});
            %ensure job IDs are consistent for all fields (in each sat and level)
            equal_idx=obj.vector_tr(@simpledata.isequal_multiple,...
              names,...
              product.mdget('jobid_col'),...
              str.tablify([7,6,10],['GRACE-',sats{s}],levels{i},product.dataname.type)...
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
%         %loop over all sat, add epochs at day boundaries and build fstep time domain
%         for i=1:numel(levels)
%           for j=1:numel(fields)
%             for s=1:numel(sats)
%               %save time series into dedicated var
%               ts_now=obj.sat_get(product.dataname.type,levels{i},fields{j},sats{s});
%               %create epochs at day boundaries
%               t_days=dateshift(ts_now.t(1),'start','day'):days(1):dateshift(ts_now.t(end),'end','day');
%               %add day boundaries
%               ts_now=ts_now.t_merge(t_days);
%               %build fstep time domain
%               ts_now=ts_now.fstep(seconds(1));
%               %save time series back to object
%               obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},ts_now);
%             end
%           end
%         end

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
                disp(str.tablify([8,10,6,6,1],'Checking',product.dataname.type,levels{i},fields{j},sats{s}))
                %save time series into dedicated var
                ts_now=obj.sat_get(product.dataname.type,levels{i},fields{j},sats{s});
                %forget about epochs that have been artificially inserted to represent gaps and end of arcs
                idx1=find(diff(ts_now.t)>seconds(1));
                %get arc lenths
                al=ts_now.y(idx1,3);
                %get consecutive time difference
                dt=seconds(diff(ts_now.t(idx1)));
                %find arcs that span over time stamps
                bad_idx=find(al(1:end-1)-dt>1); %no abs here!
                %report if any such epochs have been found
                if ~isempty(bad_idx)
                  msg=cell(1,min([10,numel(bad_idx)])+1);
                  msg{1}=str.tablify(20,'idx','arc init t','arc length','succ time diff','delta arc len');
                  for k=1:numel(msg)-1
                    idx=idx1(bad_idx(k));
                    msg{k+1}=str.tablify(20,...
                      idx1(bad_idx(k)),...
                      ts_now.t(idx),...
                      al(bad_idx(k)),...
                      dt(bad_idx(k)),...
                      al(bad_idx(k))-dt(bad_idx(k))...
                    );
                  end
                  disp([....
                    ': found ',num2str(numel(bad_idx)),' arc lengths (3rd column) longer than ',...
                    ' difference between consecutive time stamps (4th column):',10,...
                    strjoin(msg,'\n')
                  ])
%                   mask=ts_now.mask;
%                   mask(idx1(bad_idx))=false;
%                   obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},ts_now.mask_and(mask).mask_update);
                end
              end
            case {'aak','accatt'}
              %this check ensures that the t0 value is the same as the start of the arc
              for j=1:numel(fields)
                %some fields do not have t0
                if ~any(fields{j}(end)=='DQ')
                  disp(str.tablify([8,10,6,6,1],'Skipping',product.dataname.type,levels{i},fields{j},sats{s}))
                  continue
                end
                disp(str.tablify([8,10,6,6,1],'Checking',product.dataname.type,levels{i},fields{j},sats{s}))
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
                  msg{1}=str.tablify(20,'idx','arc init time','MJD','delta time');
                  for k=1:numel(msg)-1
                    idx=idx1(bad_idx(k));
                    msg{k+1}=str.tablify(20,...
                      idx1(bad_idx(k)),...
                      ts_now.t(idx1(bad_idx(k))),...
                      t0(bad_idx(k)),...
                      ts_now.t(idx1(bad_idx(k)))-t0(bad_idx(k))...
                    );
                  end
                  disp([...
                    'found ',num2str(numel(bad_idx)),' arc init time (2nd column) different than the',...
                    ' MJD reported in the data (3rd column):',10,...
                    strjoin(msg,'\n')...
                  ])
%                   mask=ts_now.mask;
%                   mask(idx1(bad_idx))=false;
%                   obj=obj.sat_set(product.dataname.type,levels{i},fields{j},sats{s},ts_now.mask_and(mask).mask_update);
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
    function import_calpar_debug_plots(debug)
      if ~exist('debug','var') || isempty(debug)
        debug=false;
      end
      %get current git version
      [status,timetag]=system(['git log -1 --format=%cd --date=iso-local ',mfilename,'.m']);
      %get rid of timezone and leading trash
      timetag=timetag(9:27);
      %sanity
      assert(status==0,[mfilename,': could not determine git time tag'])
      %create dir for plots
      plot_dir=fullfile('plot','import_calpar_debug_plots',timetag);
      if isempty(dir(plot_dir)); mkdir(plot_dir); end

      %load calibration parameters
      a=datastorage('debug',debug).init('grace.calpar_csr','plot_dir',plot_dir);
      %retrieve product info
      product=a.mdget(datanames('grace.calpar_csr'));
      %define start/stop pairs and level
      i=0;ssl=struct([]);
      i=i+1; ssl(i).field='AC0X';
      ssl(i).start=datetime('2002-08-06 00:00:00');
      ssl(i).stop =datetime('2002-08-06 23:59:59');
      i=i+1; ssl(i).field='AC0X';
      ssl(i).start=datetime('2002-08-16 00:00:00');
      ssl(i).stop =datetime('2002-08-18 23:59:59');
      i=i+1; ssl(i).field='AC0X';
      ssl(i).start=datetime('2003-01-12 00:00:00');
      ssl(i).stop =datetime('2003-01-15 00:00:00');
      i=i+1; ssl(i).field='AC0X';
      ssl(i).start=datetime('2003-11-29 00:00:00');
      ssl(i).stop =datetime('2003-12-02 00:00:00');
      i=i+1; ssl(i).field='AC0XD';
      ssl(i).start=datetime('2006-06-03 00:00:00');
      ssl(i).stop =datetime('2006-06-07 00:00:00');
      i=i+1; ssl(i).field='AC0XD';
      ssl(i).start=datetime('2006-06-12 00:00:00');
      ssl(i).stop =datetime('2006-06-19 00:00:00');
      i=i+1; ssl(i).field='AC0X';
      ssl(i).start=datetime('2012-06-30 00:00:00');
      ssl(i).stop =datetime('2012-07-03 00:00:00');
      sats  =product.mdget('sats');
      %loop over the data
      for i=1:numel(ssl)
        for s=1:numel(sats)
          dataname_now={'grace','calpar_csr','',ssl(i).field,sats{s}};
          p=a.trim(ssl(i).start,ssl(i).stop).plot(...
            datanames(dataname_now),...
            'plot_file_full_path',false,...
            'plot_together','level'...
          );
        end
      end
    end
    function obj=compute_calmod(obj,dataname,varargin)
      % parse mandatory arguments
      p=inputParser;
      p.addRequired('dataname', @(i) isa(i,'datanames'));
      p.parse(dataname);
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
      coords    =calparp.mdget('coords');
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
            %propagate it
            obj=obj.sat_set(dataname.type,dataname.level,levels{l},sats{s},calmod);
          end
        end
      end
    end
    function obj=import_acc_l1b(obj,dataname,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('dataname',@(i) isa(i,'datanames'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
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
      [~,timestamplist]=product.file('data',varargin{:},...
        'start',p.Results.start,...
        'stop', p.Results.stop...
      );
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
      p.addRequired('dataname',@(i) isa(i,'datanames'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
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
      [~,timestamplist]=product.file('data',varargin{:},...
        'start',p.Results.start,...
        'stop', p.Results.stop...
      );
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