classdef csr
  methods(Static)
    function timetag=gitversion
      %get current git version
      [status,timetag]=system(['git -C ',fileparts(which(mfilename)),' log -1 --format=%cd --date=iso-local ',mfilename,'.m']);
      %get rid of timezone and leading trash
      timetag=timetag(9:27);
      %sanity
      assert(status==0,[mfilename,': could not determine git time tag'])
    end
    function log(msg)
      logname=fullfile(fileparts(mfilename),'import_calpar.log');
      if ~exist('msg','var')
        if ~isempty(dir(logname))
          system(['mv -v ',logname,' ',strrep(logname,'.log',''),'.',datestr(datetime('now'),30),'.log']);
        end
      else
        fid = fopen(logname,'a');  
        fprintf(fid,[strjoin(msg,'\n'),'\n']);
        fclose(fid);
      end
    end
    function report(debug,idx,context,id,labels,data,log_flag)
      if ~exist('log_flag','var') ||isempty(log_flag)
        log_flag=true;
      end
      if isempty(idx); return; end
      [~,ids]=fileparts(id);
      if isempty(ids); ids=id; end
      msg=cell(1,numel(idx)+2);
      msg{1}=str.tablify([36,30,1,5],[context,' for'],ids,':',num2str(numel(idx)));
      msg{2}=str.tablify(20,'data','idx',labels{:});
      for k=1:numel(idx)
        msg_data=cell(1,numel(data));
        for l=1:numel(data)
          if numel(data{l})==numel(idx)
            msg_data{l}=data{l}(k);
          else
            msg_data{l}=data{l}(idx(k));
          end
        end
        msg{k+2}=str.tablify(20,ids,idx(k),msg_data{:});
      end
      if debug;disp(strjoin(msg(1:min([20,numel(msg)])),'\n')); else disp(msg{1}); end; 
      if log_flag;csr.log(msg);end
    end
    function obj=import_calpar(obj,product,varargin)
      %open log file
      csr.log
      % parse mandatory arguments
      p=inputParser;
      p.addRequired('product',  @(i) isa(i,'dataproduct'));
      p.addParameter('debugdate', [], @(i) ischar(i) || isempty(i));
      p.parse(product);
      %get names of parameters and levels
      levels     =product.mdget('levels');
      fields     =product.mdget('fields');
      fields_out =product.mdget('fields_out');
      sats       =product.mdget('sats');
      bias_files =product.mdget('bias_files');
      param_col  =product.mdget('param_col');
      jobid_col  =product.mdget('jobid_col');
      arclen_col =product.mdget('arclen_col');
      t0_col     =product.mdget('t0_col');
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
              switch levels{i}
              case {'aak','accatt'}
                %need to ensure timestamps are in agreement with t0
                assert(tmp.(sats{s}).width < t0_col || all(tmp.(sats{s}).mjd==tmp.(sats{s}).y(:,t0_col)),[mfilename,':',...
                  'discrepancy between time domain and t0.'])
              end
              %add long-term biases
              t=tmp.(sats{s}).mjd-ltb.(sats{s})(2,1);
              tmp.(sats{s})=tmp.(sats{s}).assign(...
                [tmp.(sats{s}).y(:,param_col)+polyval(...
                  ltb.(sats{s})(:,lbt_idx),t),...
                  tmp.(sats{s}).y(:,[1:param_col-1,param_col+1:end])...
                ]...
              );
            end

            %additional processing: add end of arcs
            switch levels{i}
            case {'aak','accatt'}
              %get arc stars
              arc_starts=tmp.(sats{s}).t;
              %build arc ends
              arc_ends=[arc_starts(2:end);dateshift(arc_starts(end),'end','day')]-seconds(1);
%                 %arc ends are at maximum 24 hours after arc starts (only for those arcs starting at mid-night)
%                 fix_idx=arc_ends-arc_starts>days(1) & ...
%                   seconds(arc_starts-dateshift(arc_starts,'start','day'))<tmp.(sats{s}).t_tol;
%                 arc_ends(fix_idx)=arc_starts(fix_idx)+days(1)-seconds(1);
            case 'estim'
              %get arc stars
              arc_starts=tmp.(sats{s}).t;
              %build arc ends (arc duration given explicitly)
              arc_ends=arc_starts+seconds(tmp.(sats{s}).y(:,arclen_col))-seconds(1);

              %patch missing arc durations
              idx=find(isnat(arc_ends));
              %report edge cases
              csr.report(obj.debug,idx,'Arcs without arc length',f,...
                {'arc start','arc duraction'},...
                {arc_starts,  tmp.(sats{s}).y(:,arclen_col)}...
              )
              %fix it
              if ~isempty(idx);
                arc_ends(idx)=dateshift(arc_starts(idx),'end','day')-seconds(1);
              end
%                 
%                 %get seconds-of-day of arc ends
%                 sod_arc_ends=seconds(arc_ends-dateshift(arc_ends,'start','day'));
%                 %find the 24hrs arcs (those that have ~0 seconds of days)
%                 idx=find(sod_arc_ends<tmp.(sats{s}).t_tol);
%                 %push those arcs to midnight and remove 1 second TODOL this was 'start'
%                 arc_ends(idx)=dateshift(arc_ends(idx),'end','day')-seconds(1);
%                 %check for ilegal arc durations
%                 idx=find(sod_arc_ends>86400);
%                 csr.report(obj.debug,idx,'Arcs ending after midnight',f,...
%                   {'arc start','arc end','sod arc end'},...
%                   {arc_starts,arc_ends,sod_arc_ends}...
%                 )
%                 %fix it
%                 if ~isempty(idx)
%                   arc_ends(idx)=dateshift(arc_ends(idx),'start','day')-seconds(1);
%                 end
            end

            %bug trap
            assert(all(~isnat(arc_starts)),...
              [mfilename,': found NaT in the arc starts'])

            %compute arc day start and end
            day_starts=dateshift(arc_starts,'start','day');
            day_ends  =dateshift(arc_starts,'end',  'day');

            %arc ends cannot go over day boundaries
            idx=find(arc_ends>=day_ends);
            csr.report(obj.debug,idx,'Arc ends over day boundary',f,...
              {'curr arc start','curr arc end','day ends'},...
              {      arc_starts,      arc_ends, day_ends}...
            )
            %fix it
            if ~isempty(idx)
              arc_ends(idx)=day_ends(idx)-seconds(1);
            end

%               %fix NaTs in arc ends
%               idx=find(isnat(arc_ends));
%               %fix it
%               if ~isempty(idx)
%                 arc_ends(idx)=day_ends(idx)-seconds(1);
%               end
            %bug trap
            assert(all(~isnat(arc_ends)),...
              [mfilename,': found NaT in the arc starts/ends'])

            %surpress over-lapping arcs
            idx=find(arc_starts(2:end)-arc_ends(1:end-1)<0);
            csr.report(obj.debug,idx,'Over-lapping arcs',f,...
              {'curr arc start','curr arc end','next arc start'},...
              {     arc_starts,      arc_ends, [arc_starts(2:end);arc_starts(1)]}...
            )
            %fix it
            if ~isempty(idx)
              arc_ends(idx)=arc_starts(idx+1)-seconds(1);
            end

            %fancy stuff: handle parameters defined as arc segments
            if ~isempty(strfind(fields{j},'AC0Y'))
              %there are 8 segments per day
              periodicity=days(1)/8;
              %get day location for this parameter
              day_loc=str2double(fields{j}(end));
              %get sub-arc starts/ends
              sub_arc_starts=arc_starts+periodicity*(day_loc-1);
                sub_arc_ends=arc_starts+periodicity*(day_loc  )-seconds(1);
              %get sub arc boundaries
              sub_arc_bound_starts=max([arc_starts,day_starts],[],2);
              sub_arc_bound_ends  =min([arc_ends,  day_ends  ],[],2);
              %cap sub-arc start/ends to be within the current day
              idx={...
                find( sub_arc_starts>sub_arc_bound_ends   ),...
                find( sub_arc_starts<sub_arc_bound_starts ),...
                find( sub_arc_ends  >sub_arc_bound_ends   ),...
                find( sub_arc_ends  <sub_arc_bound_starts )...
              };
              msg={...
                'Sub-arc starts after day/arc ends',...
                'Sub-arc starts before day/arc starts',...
                'Sub-arc ends after day/arc ends',...
                'Sub-arc ends before day/arc starts'...
              };
              for k=1:numel(idx)
                csr.report(obj.debug,idx{k},msg{k},f,...
                  {'sub-arc start','sub-arc end','day start','day end'},...
                  { sub_arc_starts, sub_arc_ends, day_starts, day_ends}...
                )
                %fix it
                if ~isempty(idx{k})
                  switch k
                  case 1; sub_arc_starts(idx{k})=sub_arc_bound_ends(  idx{k});
                  case 2; sub_arc_starts(idx{k})=sub_arc_bound_starts(idx{k});
                  case 3; sub_arc_ends(  idx{k})=sub_arc_bound_ends(  idx{k});
                  case 4; sub_arc_ends(  idx{k})=sub_arc_bound_starts(idx{k});
                  end
                end
              end
              %propagate the arc extremeties
              arc_starts=sub_arc_starts;
                arc_ends=sub_arc_ends;
            end

%               %arc ends cannot be at day starts (that's the next arc start)
%               idx=find(arc_ends==day_ends);
%               if ~isempty(idx)
%                 arc_ends(idx)=arc_ends(idx)-seconds(1);
%               end

            %propagate data
            arc_start_y=tmp.(sats{s}).y;
              arc_end_y=tmp.(sats{s}).y;

            %remove arcs with zero length (only applicable to AC0Y*2-8)
            zero_len_idx=find(arc_ends-arc_starts<=0);
            csr.report(obj.debug,zero_len_idx,'Non-positive arc length',f,...
              {'arc start','arc end','arc_length'},...
              { arc_starts, arc_ends, arc_ends-arc_starts}...
            )
            if ~isempty(zero_len_idx)
              good_idx=(arc_ends-arc_starts>0);
              arc_starts =arc_starts( good_idx);
              arc_ends   =arc_ends(   good_idx);
              arc_start_y=arc_start_y(good_idx,:);
              arc_end_y  =arc_end_y(  good_idx,:);
            end

            %debug date report
            if ~isempty(p.Results.debugdate)
              rep_date=datetime(p.Results.debugdate);
              rep_delta=arc_starts-rep_date;
              rep_idx=find(abs(rep_delta)==min(abs(rep_delta)));
              rep_idx=(rep_idx(1)-8):(rep_idx(end)+8);
              csr.report(true,rep_idx,['DEBUG DATE: Arcs around ',datestr(rep_date)],f,...
                {'arc start','arc end','arc length','inter-arc gap'},...
                {arc_starts(rep_idx  ),arc_ends(  rep_idx),...
                 arc_ends(  rep_idx  )-arc_starts(rep_idx),...
                 arc_starts(rep_idx+1)-arc_ends(  rep_idx)...
                 },...
              false)
            end

%               % set the arc length to zero for arc ends
%               switch levels{i}
%               case 'estim'
%                 arc_end_y(:,arclen_col)=0;
%               end

            %build timeseries with arc starts
            arc_start_ts=simpletimeseries(arc_starts,arc_start_y,...
              'format','datetime',...
              'labels',tmp.(sats{s}).labels,...
              'units',tmp.(sats{s}).y_units,...
              'timesystem',tmp.(sats{s}).timesystem,...
              'descriptor',tmp.(sats{s}).descriptor...
            );
            %build timeseries with arc ends
            arc_end_ts=simpletimeseries(arc_ends,arc_end_y,...
              'format','datetime',...
              'labels',tmp.(sats{s}).labels,...
              'units',tmp.(sats{s}).y_units,...
              'timesystem',tmp.(sats{s}).timesystem,...
              'descriptor',['end of arcs for ',tmp.(sats{s}).descriptor]...
            );

            %augment arc starts with arc ends (only new data)
            tmp.(sats{s})=arc_start_ts.augment(arc_end_ts,'old',true);

          end

          %propagate data to object
          for s=1:numel(sats)
            obj=obj.data_set(product.dataname.set_field_path([levels(i),fields(j),sats(s)]),tmp.(sats{s}));     
          end
          %user feedback
          disp(str.tablify([15,6,3,6],'loaded data for',levels{i},'and',fields{j}))
        end
      end

      %merge cross-track accelerations together
      ac0y='AC0Y';
      field_part_list={'','D','Q'};
      %loop over all levels and sats
      for i=1:numel(levels)
        for s=1:numel(sats)
          for f=1:numel(field_part_list)
            %start with first field
            field=[ac0y,field_part_list{f},'1'];
            ts_now=obj.data_get_scalar(product.dataname.set_field_path([levels(i),field,sats(s)]));
            %loop over all other fields
            for fpl=2:8
              field=[ac0y,field_part_list{f},num2str(fpl)];
              ts_now=ts_now.augment(...
                obj.data_get_scalar(product.dataname.set_field_path([levels(i),field,sats(s)])),...
                'quiet',true,...
                'old',true,...
                'skip_gaps',true...
              );
              %debug date report
              if ~isempty(p.Results.debugdate)
                rep_date=datetime(p.Results.debugdate);
                disp(['DEBUG DATE: merge AC0Y*:',levels{i},':',sats{s},':',field,' @ ',datestr(rep_date)]);
                idx=ts_now.idx(rep_date);
                ts_now.peek((idx-10):(idx+10));
              end
            end
            %save the data
            obj=obj.data_set(product.dataname.set_field_path([levels(i),{[ac0y,field_part_list{f}]},sats(s)]),ts_now);
            %user feedback
            disp(str.tablify([29,5,3,6,3,7],'merged cross-track parameter',[ac0y,field_part_list{f}],...
              'for',levels{i},'and',['GRACE-',sats{s}]))
          end
        end
      end
      %don't need the old fields no more
      clear fields

      %add gaps
      for i=1:numel(levels)
        for j=1:numel(fields_out)
          for s=1:numel(sats)
            tmp=obj.data_get_scalar(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)]));
            end_arc_idx=[false;diff(tmp.y(:,1))==0];
                gap_idx=[diff(tmp.t)>seconds(1)+tmp.t_tol;false];
            gap_t=tmp.t(end_arc_idx & gap_idx)+seconds(1);
            %build timeseries with arc ends
            gaps=simpletimeseries(gap_t,nan(numel(gap_t),tmp.width),...
              'format','datetime',...
              'labels',tmp.labels,...
              'units',tmp.y_units,...
              'timesystem',tmp.timesystem,...
              'descriptor',['gaps for ',tmp.descriptor]...
            );
            %augment and save
            obj=obj.data_set(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)]),...
              tmp.augment(gaps,'old',true,'new',true)...
            );  
            %debug date report
            if ~isempty(p.Results.debugdate)
              rep_date=datetime(p.Results.debugdate);
              disp(['DEBUG DATE: with gaps:',levels{i},':',sats{s},':',fields_out{j},' @ ',datestr(rep_date)]);
              idx=ts_now.idx(rep_date);
              obj.data_get_scalar(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)])).peek((idx):(idx+20));
            end              
          end
        end
      end

      %loop over all sat and level to check Job IDs agreement across all output fields
      for i=1:numel(levels)
        for s=1:numel(sats)
          for j=1:numel(fields_out)-1
            dn1=product.dataname.set_field_path([levels(i),fields_out(j  ),sats(s)]);
            dn2=product.dataname.set_field_path([levels(i),fields_out(j+1),sats(s)]);
            d1=obj.data_get_scalar(dn1);
            d2=obj.data_get_scalar(dn2);
            [~,i1,i2]=intersect(d1.t,d2.t);
            bad_idx=find(...
              d1.y(i1,jobid_col) ~= d2.y(i2,jobid_col) & ...
              d1.mask(i1) & ...
              d2.mask(i2) ...
            );
            if ~isempty(bad_idx)
              n=numel(bad_idx);
              msg=cell(1,2*n+2);
              msg{1}=str.tablify([5,6,24],'found',numel(bad_idx),'Job ID inconsistencies:');
              msg{2}=str.tablify([30,6,20,12],'data name','idx','t','Job ID');
              for k=1:n
                idx=i1(bad_idx(k));
                msg{2*k+1}=str.tablify([30,6,20,12],dn1,idx,d1.t(idx),num2str(d1.y(idx,jobid_col),'%i'));
                idx=i2(bad_idx(k));
                msg{2*k+2}=str.tablify([30,6,20,12],dn2,idx,d2.t(idx),num2str(d2.y(idx,jobid_col),'%i'));
              end
              error([mfilename,':',strjoin(msg,'\n')])
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
            for j=1:numel(fields_out)
              %some fields do not have t0
              if ~any(fields_out{j}(end)=='DQ') || ~isempty(strfind(fields_out{j},'Y'))
                disp(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              disp(str.tablify([8,32],'Checking',product.str))
              %save time series into dedicated var
              ts_now=obj.data_get_scalar(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)]));
              %forget about epochs that have been artificially inserted to represent gaps and end of arcs
              %the 1e-6 parcel is needed to avoid artificially-inserted gaps that have round-off errors
              idx1=find(diff(ts_now.t)>seconds(1)+1e-6);
              %get arc lenths
              al=ts_now.y(idx1,arclen_col);
              %get consecutive time difference
              dt=[seconds(diff(ts_now.t));0]; dt=dt(idx1);
              %find arcs that span over time stamps
              bad_idx=find(al-dt>2); %no abs here!
              %report if any such epochs have been found
              csr.report(obj.debug,bad_idx,'Ilegal arc length in the data',[levels{i},'.',fields_out{j},'.',sats{s}],...
                {'global idx','arc init t','arc length','succ time diff','delta arc len'},...
                {idx1,ts_now.t(idx1),al,dt,al-dt}...
              ) %#ok<FNDSB>
            end
          case {'aak','accatt'}
            %this check ensures that the t0 value is the same as the start of the arc
            for j=1:numel(fields_out)
              %the Y parameter was constructed from multitple parameters and some fields do not have t0
              if ~any(fields_out{j}(end)=='DQ') || ~isempty(strfind(fields_out{j},'Y')) 
                disp(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              disp(str.tablify([8,32],'Checking',product.str))
              %save time series into dedicated var
              ts_now=obj.data_get_scalar(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)]));
              %forget about epochs that have been artificially inserted to represent forward steps
              idx1=find(diff(ts_now.t)>seconds(1));
              %get t0
              t0=simpletimeseries.utc2gps(datetime(ts_now.y(idx1,t0_col),'convertfrom','modifiedjuliandate'));
              %find arcs that have (much) t0 different than their first epoch
              bad_idx=find(...
                abs(ts_now.t(idx1)-t0)>seconds(1) & ...
                ts_now.mask(idx1) & ...                          %ignore gaps
                [true;diff(ts_now.y(idx1,jobid_col))~=0] ...     %ignore epochs inside the same arc
              );
              %report if any such epochs have been found
              csr.report(obj.debug,bad_idx,'Ilegal t0 in the data',[levels{i},'.',fields_out{j},'.',sats{s}],...
                {'global idx','arc init time','t0','delta time'},...
                {idx1,ts_now.t(idx1),t0,ts_now.t(idx1)-t0}...
              ) %#ok<FNDSB>bo
            end
          end
        end
      end
    end
    function import_calpar_debug_plots(debug,name)
      if ~exist('debug','var') || isempty(debug)
        debug=false;
      end
      if ~exist('name','var') || isempty(name)
        name='grace.calpar.csr';
      end
      %create dir for plots
      plot_dir=fullfile(dataproduct.default_list.plot_dir,'import_calpar_debug_plots',csr.gitversion);
      if isempty(dir(plot_dir)); mkdir(plot_dir); end

      %load calibration parameters
      a=datastorage('debug',debug).init(name,'plot_dir',plot_dir);
      %retrieve product info
      sats=a.product_get(name).mdget('sats');
      %define start/stop pairs and level
      i=0;ssl=struct([]);
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z'};
      ssl(i).start=datetime('2002-08-06 00:00:00');
      ssl(i).stop =datetime('2002-08-06 23:59:59');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z'};
      ssl(i).start=datetime('2002-08-06 00:00:00');
      ssl(i).stop =datetime('2002-08-06 23:59:59');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z'};
      ssl(i).start=datetime('2002-08-16 00:00:00');
      ssl(i).stop =datetime('2002-08-18 23:59:59');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z'};
      ssl(i).start=datetime('2003-01-12 00:00:00');
      ssl(i).stop =datetime('2003-01-15 00:00:00');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z'};
      ssl(i).start=datetime('2003-11-29 00:00:00');
      ssl(i).stop =datetime('2003-12-02 00:00:00');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z','AC0XD','AC0YD','AC0ZD','AC0XQ','AC0YQ','AC0ZQ'};
      ssl(i).start=datetime('2006-06-03 00:00:00');
      ssl(i).stop =datetime('2006-06-07 00:00:00');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z','AC0XD','AC0YD','AC0ZD','AC0XQ','AC0YQ','AC0ZQ'};
      ssl(i).start=datetime('2006-06-12 00:00:00');
      ssl(i).stop =datetime('2006-06-19 00:00:00');
      i=i+1; ssl(i).field={'AC0X','AC0Y','AC0Z'};
      ssl(i).start=datetime('2012-06-30 00:00:00');
      ssl(i).stop =datetime('2012-07-03 00:00:00');
      %loop over the data
      for i=1:numel(ssl)
        p=a.trim('start',ssl(i).start,'stop',ssl(i).stop);
        for f=1:numel(ssl(i).field)
          for s=1:numel(sats)
            p.plot(...
              datanames(name).set_field_path({'*',ssl(i).field{f},sats{s}}),...
              'plot_together',{'aak','accatt','estim'}...
            );
          end
        end
      end
    end
    function obj=compute_calmod(obj,product,varargin)
      % parse mandatory arguments
      p=inputParser;
      p.addRequired('product', @(i) isa(i,'dataproduct'));
      p.parse(product);
      %paranoid sanity
      if product.nr_sources~=2
        error([mfilename,': number of sources in product ',product.name,...
          ' is expected to be 2, not ',num2str(product.nr_sources),'.'])
      end
      
      keyboard
      
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
        acc=obj.data_get_scalar(l1baccp.dataname.set_field_path([l1baccp.dataname.level,l1baccp.dataname.field,sats{s}]));
        %loop over all levels
        for l=1:numel(levels)
          %handle exceptions (also deals with non-existing data)
          if ~isa(acc,'simpletimeseries')
            %patch nan calibration model
            calmod=simpletimeseries(...
              [obj.start;obj.stop],...
              nan(2,numel(coords)),...
              'descriptor',['calibration model ',levels{l},' GRACE-',upper(sats{s}),' (empty)']...
            );
            disp(['Skipping  the ',calmod.descriptor])
          else
            %init models container
            calmod=simpletimeseries(acc.t,zeros(acc.length,numel(coords))).copy_metadata(acc);
            calmod.descriptor=['calibration model ',levels{l},' GRACE-',upper(sats{s})];
            disp(['Computing the ',calmod.descriptor])
            for c=1:numel(coords)
              %retreive data
              cal=struct(...
                'ac0' ,obj.data_get_scalar(calparp.dataname.set_field_path({levels{l},['AC0',coords{c}    ],sats{s}})),...
                'ac0d',obj.data_get_scalar(calparp.dataname.set_field_path({levels{l},['AC0',coords{c},'D'],sats{s}})),...
                'ac0q',obj.data_get_scalar(calparp.dataname.set_field_path({levels{l},['AC0',coords{c},'Q'],sats{s}}))...
              );
              %get field names (too lazy to type)
              fields=fieldnames(cal);
              %loop over all fields
              for f=1:numel(fields)
                %interpolate to time domain of measurements
                cal.(fields{f})=cal.(fields{f}).interp(acc.t);
                %sanity
                assert(~isempty(cal.(fields{f})),[mfilename,...
                  ': ',(fields{f}),' data is not available to perform this operation.'])
                %remove transient calpars, this is indicative of a gap (interpolation is done blindly over any gap length)
                cal.(fields{f})=cal.(fields{f}).mask_and([cal.(fields{f}).mask(1);diff(cal.(fields{f}).y(:,1))==0]);
                %retrieve time domain (it is the same for all cal pars)
                if strcmp(levels{l},'estim')
                  %start of day
                  sod=dateshift(cal.(fields{f}).t,'start','day');
                  %second of arc
                  soa=sod+seconds(cal.(fields{f}).y(:,4));
                  %time domain for the calibration model: units are days and zero epoch is the start of the arc
                  t.(fields{f})=days(acc.t-soa);
                else
                  %TODO: é preciso arranjar isto, o start arc tem de ser definido algures (para aak e accatt)
                  t.(fields{f})=days(acc.t-simpletimeseries.ToDateTime(cal.(fields{f}).y(:,end),'modifiedjuliandate'));
                end
              end
              %paranoid sanity check
              good_idx=~isnan(t.ac0);
              if any(t.ac0(good_idx)~=t.ac0d(good_idx)) || any(t.ac0(good_idx)~=t.ac0q(good_idx))
                error([mfilename,': calibration time domain inconsistent between parameters, debug needed!'])
              end
              %build calibration model
              calmod=calmod.set_cols(c,...
                cal.ac0.cols( param_col                 )+...
                cal.ac0d.cols(param_col).times(t.ac0d   )+...
                cal.ac0q.cols(param_col).times(t.ac0q.^2)...
              );
              %make debug plots
              if product.mdget('debug_plot')
                font_size_args={...
                  'plot_fontsize_title',24,...
                  'plot_fontsize_axis',18, ...
                  'plot_fontsize_label',20 ...
                };
                filename=product.file('plot','start',obj.start,'stop',obj.stop,...
                  'use_storage_period',false,'timestamp',true,...
                  'suffix',[sats{s},'.',coords{c}]...
                );
                %add level
                [p,f,e]=fileparts(filename{1});
                filename=fullfile(p,levels{l},[f,e]);
                if ~exist(fileparts(filename),'dir'); mkdir(fileparts(filename)); end
                figure;
                subplot(2,2,1)
                for f=1:numel(fields)
                  plot(acc.t,t.(fields{f})), hold on
                end
                grid on
                title('time domain')
                product.enforce_plot(font_size_args{:})
                legend(fields)
                for f=1:numel(fields)
                  subplot(2,2,f+1)
                  plot(acc.t,cal.(fields{f}).cols(param_col).y)
                  title(cal.(fields{f}).labels(param_col))
                  grid on
                  product.enforce_plot(font_size_args{:})
                end
                disp(['Created plot ',filename])
                saveas(gcf,filename)
              end
            end
          end
          %propagate it
          obj=obj.data_set(product.dataname.set_field_path([product.dataname.level,levels(l),sats(s)]),calmod);
        end
      end
    end
    function obj=import_acc_l1b(obj,product,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',@(i) isa(i,'dataproduct'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(product,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
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
        %load (and save the data in mat format, as handled by simpletimeseries.import
        obj=obj.data_set(product.dataname.set_field_path({product.dataname.level,product.dataname.field,sats{s}}),...
          simpletimeseries.import(infile,'cut24hrs',false)...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
    function obj=import_acc_mod(obj,product,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',@(i) isa(i,'dataproduct'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(product,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
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
        obj=obj.data_set(product.dataname.set_field_path({product.dataname.level,product.dataname.field,sats{s}}),...
          simpletimeseries.import(infile,'cut24hrs',true)...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
    function obj=import_acc_resid_test(obj,product,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',@(i) isa(i,'dataproduct'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(product,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
      %retrieve relevant parameters
      sats  =product.mdget('sats');
      coords=product.mdget('coords');
      indir =product.mdget('import_dir');
      infile_template =product.mdget('filename');
      %gather list of daily data files
      [outfiles,timestamplist]=product.file('data',varargin{:},...
        'start',p.Results.start,...
        'stop', p.Results.stop...
      );
      %loop over the satellites
      for f=1:numel(outfiles)
        if ~exist(outfiles{f},'file')
          for c=1:numel(coords)
            for s=1:numel(sats)
              infile=file.wildcard(...
                strrep(strrep(strrep(infile_template,...
                  '<sat>',sats{s}),...
                  '<coord>',coords{c}),...
                  '<date>',datestr(timestamplist(f),'yy-mm-dd'))...
              );
              d=simpletimeseries.import(infile,'cut24hrs',false);
              simpletimeseries.import(infile,'cut24hrs',false)
            end
          end
        
          %save data
          s=obj.datatype_get(product.dataname.type); %#ok<*NASGU>
          save(outfiles{f},'s');
          clear s
        else
          %load data
          load(outfiles{f},'s');
          fields=fieldnames(s); 
          for i=1:numel(fields)
            obj=obj.data_set(product,s.(fields{i}));
          end
        end

          
      end
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
        obj=obj.data_set(product.dataname.set_field_path({product.dataname.level,product.dataname.field,sats{s}}),...
          simpletimeseries.import(infile,'cut24hrs',true)...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
    function rm_data(mode,varargin)
      if ~exist('mode','var') || isempty(mode)
        mode='all';
      end
      switch lower(mode)
      case 'all'
        %define all data to be removed
%           'grace.acc.l1b.nrtdm',...
%           'grace.acc.mod.nrtdm',...
%           'grace.acc.l1b.csr',...
%           'grace.acc.mod.csr',...
%           'graceacccal',...
        d={...
          'grace.acc.cal_csr',...
          'grace.calpar_csr',...
          'grace.calpar_csr_corr',...
          'grace.calpar_csr_stats'...
        };
        %recursive call
        for i=1:numel(d)
          csr.rm_data(d{i});
        end
      case 'graceacccal'
        datadir=fullfile(getenv('HOME'),'data','csr','GraceAccCal');
        system(['rm -fv ',fullfile(datadir,'to-delete','*.mat')])
        system(['mv -fv ',fullfile(datadir,'*.mat'),' ',fullfile(datadir,'to-delete')])
      case {'import_calpar_debug_plots','calpar_debug_plots'}
        plot_dir=fullfile('plot',lower(mode),csr.gitversion);
        disp(plot_dir)
        if ~isempty(dir(plot_dir))
          [status,result]=system(['rm -fvr "',plot_dir,'"']);
          assert(status==0,['error removing ',mode,'.'])
          disp(result)
        end
      otherwise
        dataproduct(mode,'metadata_dir',obj.metadata_dir).rm_data(varargin{:});
      end
    end
    function calpar_debug_plots(debug)
      if ~exist('debug','var') || isempty(debug)
        debug=false;
      end
      %create dir for plots
      plot_dir=fullfile('plot','calpar_debug_plots',csr.gitversion);
      if isempty(dir(plot_dir)); mkdir(plot_dir); end
      
      %define list of days to plot
      lod=datetime({...
        '2008-02-24',...
        '2008-02-25'...
%         '2014-09-13',...
%         '2014-09-14'...
%         '2014-08-01',...
%         '2014-08-02'...
%         '2013-02-26',...
%         '2013-02-25',...  
%         '2013-02-24',...
%         '2013-02-23',...
%         '2013-02-22',...
%         '2013-02-21'...
%         '2013-02-20',...
%         '2002-05-04',... %nominal
%         '2002-05-11',...
%         '2002-09-07',...
%         '2002-04-15',... %exception
%         '2002-05-16',...
%         '2002-09-28',...
%         '2002-09-30'...
%         '2002-04-27',... %exception
%         '2002-08-06',... %exception
%         '2002-08-16','2002-08-17',...
%         '2003-01-13','2003-01-14',...
%         '2003-11-20',...
%         '2006-06-14','2006-06-15',...
%         '2002-04-15',...
%         '2002-08-03',...
%         '2002-08-26',...
%         '2002-04-27',...
%         '2002-09-30'...
      });
      %loop over all requested days
      for i=1:numel(lod)
        start=lod(i);
        stop=lod(i)+days(1)-seconds(1);
        %plot it
        datastorage('debug',debug,'start',start,'stop',stop).init('grace.acc.cal_csr_plots','plot_dir',plot_dir);
      end
    end
    function out=test
      out=datastorage('debug',true).init('grace.calpar.csr');
    end
  end
end