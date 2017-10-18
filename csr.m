classdef csr
  %static
  properties(Constant)
    arc_timeline_file='/Users/teixeira/data/csr/GraceAccCal/arctimeline.GraceAccCal';
  end
  methods(Static)
    %% utilities
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
        fid = file.open(logname,'a');  
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
          switch numel(data{l})
          case 1
            msg_data{l}=data{l};
          case numel(idx)
            msg_data{l}=data{l}(k);
          otherwise
            msg_data{l}=data{l}(idx(k));
          end
        end
        msg{k+2}=str.tablify(20,ids,idx(k),msg_data{:});
      end
      if debug;disp(strjoin(msg(1:min([20,numel(msg)])),'\n')); else disp(msg{1}); end; 
      if log_flag;csr.log(msg);end
    end
    %% debugging
    function out=debug_parameters
      %define start/stop pairs and level
      i=0;out=struct([]);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-04-15 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-04-27 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-05-04 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-05-11 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-05-16 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-08-03 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-08-06 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-08-16 00:00:00'); out(i).stop=out(i).start+days(2)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-08-26 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-09-07 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-09-28 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2002-09-30 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2003-01-12 00:00:00'); out(i).stop=datetime('2003-01-15 23:59:59');
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2003-11-20 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2003-11-21 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
%       i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
%       out(i).start=datetime('2003-11-29 00:00:00'); out(i).stop=datetime('2003-12-02 23:59:59');
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z','AC0XD','AC0YD','AC0ZD','AC0XQ','AC0YQ','AC0ZQ'};
      out(i).start=datetime('2006-06-03 00:00:00'); out(i).stop=datetime('2006-06-07 23:59:59');
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z','AC0XD','AC0YD','AC0ZD','AC0XQ','AC0YQ','AC0ZQ'};
      out(i).start=datetime('2006-06-12 00:00:00'); out(i).stop=datetime('2006-06-19 23:59:59');
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2008-02-24 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2008-02-25 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2012-06-30 00:00:00'); out(i).stop=datetime('2012-07-03 23:59:59');
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2013-02-20 00:00:00'); out(i).stop=datetime('2013-02-26 23:59:59');
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2014-08-01 00:00:00'); out(i).stop=out(i).start+days(2)-seconds(1);
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2014-09-13 00:00:00'); out(i).stop=out(i).start+days(2)-seconds(1);
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2015-08-14 00:00:00'); out(i).stop=out(i).start+days(1)-seconds(1);
      i=i+1; out(i).field={'AC0X','AC0Y','AC0Z'};
      out(i).start=datetime('2017-01-07 00:00:00'); out(i).stop=out(i).start+days(3)-seconds(1);

    end
    function out=debug_plots(mode,debug,plot_dir)
      if ~exist('debug','var') || isempty(debug)
        debug=true;
      end
      if ~exist('mode','var') || isempty(mode)
         mode='import_calpar';
      end
      if ~exist('plot_dir','var') || isempty(plot_dir)
         plot_dir=fullfile(plot_dir,['debug_plots_',mode],csr.gitversion);
      end
      %create dir for plots
      if isempty(dir(plot_dir)); mkdir(plot_dir); end
      %get dates of plots
      ssl=csr.debug_parameters;
      %outputs
      out=cell(size(ssl));
      %additional arguments
      args={};
      %branch on mode
      switch mode
      case 'all'
        for i={'import_calpar','compute_calmod','calpar','calpar-gps'}
          csr.debug_plots(i{1});
          close all
        end
        %done
        return
      case 'import_calpar'
        name='grace.calpar.csr';
        %load calibration parameters
        out=datastorage('debug',debug).init(name,'plot_dir',plot_dir);
        %retrieve product info
        sats=out.product_get(name).mdget('sats');
        %loop over the data
        for i=1:numel(ssl)
          p=out.trim('start',ssl(i).start,'stop',ssl(i).stop);
          for f=1:numel(ssl(i).field)
            for s=1:numel(sats)
              p.plot(...
                datanames(name).set_field_path({'*',ssl(i).field{f},sats{s}}),...
                'plot_together',{'aak','accatt','estim'}...
              );
            end
          end
        end
        %done
        return
      case 'calpar-gps'
        start=datetime('2017-01-01');
        stop =datetime('2017-01-02');
        out=datastorage('start',start,'stop',stop,'debug',debug);
        for p={'grace.acc.l1b.csr','grace.acc.calmod.csr','grace.acc.cal.csr','grace.acc.mod.nrtdm','grace.acc.mod.csr'}
          out=out.init(p{1});
        end
        %done
        return
      case 'calpar';               name='grace.acc.cal.csr.plots';
      case 'compute_calmod';       name='grace.acc.calmod.csr';       args={'debug_plot',true};
      case 'estimate_poly_calmod'; name='grace.acc.cal.poly0.plots';
      case 'psd';                  name='grace.acc.psd.plots';
      otherwise
        error(['unknown mode ''',mode,'''.'])  
      end
      for i=1:numel(ssl)
        out{i}=datastorage('start',ssl(i).start,'stop',ssl(i).stop,'debug',debug).init(name,'plot_dir',plot_dir,args{:});
      end
    end
    function out=test(mode,start)
      if ~exist('mode','var') || isempty(mode)
        mode='poly_calmod';
        mode='calpar';
      end
      if ~exist('start','var') || isempty(start)
%         start=datetime('2008-02-26');
        start=datetime('2002-04-15');
      end
      stop =start+days(1)-seconds(1);
      %translate modes to metadata names
      switch mode
      case 'calpar';        name='grace.calpar.csr';
      case 'import_acc_l1b';name='grace.acc.l1b.csr';
      case 'calmod';        name='grace.acc.calmod.csr';
      case 'mod';           name='grace.acc.mod.csr';
      case 'nrtdm';         name='grace.acc.mod.nrtdm';
      case 'calacc';        name='grace.acc.cal.csr';
      case 'poly_calmod';   name='grace.acc.calmod.poly0';
      otherwise
        error(['unknown mode ''',mode,'''.'])  
      end
      out=datastorage('debug',true,'start',start,'stop',stop).init(name,'recompute',true);
      out.peek
    end
    %% long-term biases
    function out=ltb_data(filename)
      %original:
      %  52640.00
      % -0.54818E-06 -0.24827E-10
      %  0.86557E-05  0.17811E-08
      % -0.77612E-06  0.49781E-10
      %internal:
      % MJD        AC0X         AC0Y*        AC0Z
      % 0         -0.24827E-10  0.17811E-08  0.49781E-10  %linear term
      % 52640.00  -0.54818E-06  0.86557E-05 -0.77612E-06  %bias term
      out=flipud(transpose(dlmread(filename)));
    end
    function out=ltb(t,sat,field,scl,fn)
      %defaults
      if ~exist('field','var') || isempty(field)
        field='all';
      end
      if ~exist('scl','var') || isempty(scl)
        scl=1;
      end
      if ~exist('filename','var') || isempty(fn)
        switch sat
        case 'A'; fn=[getenv('HOME'),'/data/csr/corral-tacc/input/bsA2003'];
        case 'B'; fn=[getenv('HOME'),'/data/csr/corral-tacc/input/bsB2003'];
        otherwise; error(['Unknown sat ''',sat,'''.'])
        end
      elseif iscellstr(fn)
        switch sat
        case 'A'; fn=fn{1};
        case 'B'; fn=fn{2};
        otherwise; error(['Unknown sat ''',sat,'''.'])  
        end
      end
      %translate field to data index
      switch upper(field)
      case {'X','AC0X'}
        lbt_idx=2;
      case {'Y','AC0Y1','AC0Y2','AC0Y3','AC0Y4','AC0Y5','AC0Y6','AC0Y7','AC0Y8'}
        lbt_idx=3;
      case {'Z','AC0Z'}
        lbt_idx=4;
      case {'ALL','-ALL'}
        %if mode is '-all', then "remove" the long-term bias
        if strcmpi(field,'-ALL');scl=-scl;end
        x=csr.ltb(t,sat,'x',scl,fn);
        y=csr.ltb(t,sat,'y',scl,fn);
        z=csr.ltb(t,sat,'z',scl,fn)';
        out=x.glue(...
            y.glue(....
            z));
        out.descriptor=['LTB read from ',str.clean(fn,'file')];
        return
      otherwise
        out=[];
        return
      end
      %get the data
      ltb_data=csr.ltb_data(fn);
      %number of day since reference epoch (in the ltb file)
      d=time.mjd(t)-ltb_data(2,1);
      %create timeseries object
      out=simpletimeseries(...
        t,...
        scl*polyval(ltb_data(:,lbt_idx),d),...
        'format','datetime',...
        'timesystem','gps',...
        'labels',{field},...
        'units',{'m/s^2'},...
        'descriptor',[field,' LTB read from ',str.clean(fn,'file')]...
      );
    end
    function out=ltb_apply(product,in,sat,field)
      %get relevant parameters
        bias_files=product.mdget('bias_files');
      calpar_scale=product.mdget('calpar_scale','default',1);
         ltb_scale=product.mdget(   'ltb_scale','default',1);
      %get long-term biases
      ltb=csr.ltb(in.t,sat,field,ltb_scale,bias_files);
      %add long-term biases (unless ltb is empty, which means this field had no ltb)
      if isempty(ltb)
        out=in;
      else
        out=in.scale(calpar_scale)+ltb;
      end
    end
    %% constructors
    function obj=import_calpar(obj,product,varargin)
      %open log file
      csr.log
      %parse optional arguments
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('debugdate', [], @(i) ischar(i) || isempty(i));
      p.parse(varargin{:});
      %get names of parameters and levels
      levels     =product.mdget('levels');
      fields     =product.mdget('fields');
      fields_out =product.mdget('fields_out');
      sats       =product.mdget('sats');
      param_col  =product.mdget('param_col');
      jobid_col  =product.mdget('jobid_col');
      arclen_col =product.mdget('arclen_col');
      t0_col     =product.mdget('t0_col');

      %load data
      for i=1:numel(levels)
        for j=1:numel(fields)
          tmp=struct('A',[],'B',[]);
          for s=1:numel(sats)
            %read L1B data
            f=fullfile(product.mdget('import_dir'),['gr',sats{s},'.',fields{j},'.',levels{i},'.GraceAccCal']);
            tmp.(sats{s})=simpletimeseries.import(f,'cut24hrs',false);
            %apply long-term bias
            tmp.(sats{s})=tmp.(sats{s}).set_cols(param_col,...
              csr.ltb_apply(product,tmp.(sats{s}).get_cols(param_col),sats{s},fields{j})...
            );
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
          str.say(str.tablify([15,6,3,6],'loaded data for',levels{i},'and',fields{j}))
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
            %rename the relevant object fields to remove the '1'
            ts_now.labels=strrep(ts_now.labels,field,[ac0y,field_part_list{f}]);
            ts_now.descriptor=strrep(ts_now.descriptor,field,[ac0y,field_part_list{f},'[1-8]']);
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
                str.say('DEBUG DATE: merge AC0Y*:',levels{i},':',sats{s},':',field,' @ ',datestr(rep_date));
                idx=ts_now.idx(rep_date);
                ts_now.peek((idx-10):(idx+10));
              end
            end
            %save the data
            obj=obj.data_set(product.dataname.set_field_path([levels(i),{[ac0y,field_part_list{f}]},sats(s)]),ts_now);
            %user feedback
            str.say(str.tablify([29,5,3,6,3,7],'merged cross-track parameter',[ac0y,field_part_list{f}],...
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
            ts_now=obj.data_get_scalar(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)]));
            %debug date report
            if ~isempty(p.Results.debugdate)
              rep_date=datetime(p.Results.debugdate);
              str.say('DEBUG DATE: w/out gaps:',levels{i},':',sats{s},':',fields_out{j},' @ ',datestr(rep_date));
              idx=ts_now.idx(rep_date);
              ts_now.peek((idx-10):(idx+10));
            end
            %get end of arcs and non-consecutive time indexes
            end_arc_idx=[false;diff(ts_now.y(:,1))==0];
                gap_idx=[diff(ts_now.t)>seconds(1)+ts_now.t_tol;false];
            %extend calibration parameters into the gap
            gap_start_idx=find(end_arc_idx & gap_idx);
            gap_stop_idx=gap_start_idx+1;
%             ext_len=minutes(10);
%             extension=min( ts_now.t(gap_stop_idx)-ts_now.t(gap_start_idx),ext_len*ones(size(gap_start_idx))*2 )/2;
            extension=ts_now.t(gap_stop_idx)-ts_now.t(gap_start_idx);
            ts_now.t(gap_start_idx)=ts_now.t(gap_start_idx)+time.round_seconds(extension/2);
            ts_now.t(gap_stop_idx )=ts_now.t(gap_start_idx)+seconds(1);
            %build timeseries with arc ends
            gap_t=ts_now.t(gap_start_idx)+seconds(1);
            gaps=simpletimeseries(gap_t,nan(numel(gap_t),ts_now.width),...
              'format','datetime',...
              'labels',ts_now.labels,...
              'units',ts_now.y_units,...
              'timesystem',ts_now.timesystem,...
              'descriptor',['gaps for ',ts_now.descriptor]...
            );
            %augment (keep it separate from saving, so that date report works as expected)
            ts_now=ts_now.augment(gaps,'old',true,'new',true);
            %debug date report
            if ~isempty(p.Results.debugdate)
              rep_date=datetime(p.Results.debugdate);
              str.say('DEBUG DATE: with gaps:',levels{i},':',sats{s},':',fields_out{j},' @ ',datestr(rep_date));
              idx=ts_now.idx(rep_date);
              ts_now.peek((idx-10):(idx+10));
            end              
            %save
            obj=obj.data_set(product.dataname.set_field_path([levels(i),fields_out(j),sats(s)]),ts_now);  
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
                str.say(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              str.say(str.tablify([8,32],'Checking',product.str))
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
                str.say(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              str.say(str.tablify([8,32],'Checking',product.str))
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
    function obj=import_acc_l1b(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      %retrieve relevant parameters
      sats =product.mdget('sats');
      indir=product.mdget('import_dir');
      version=product.mdget('version');
      %gather list of daily data files
      [~,timestamplist]=product.file('data',varargin{:},...
        'start',obj.start,...
        'stop', obj.stop...
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
        %load the data in mat format, as handled by simpletimeseries.import
        out=simpletimeseries.import(infile,'cut24hrs',false);
        %save the data, if not empty
        if ~isempty(out)
          obj=obj.data_set(product.dataname.set_field_path(sats(s)),out);
        else
          str.say('No data in file(s):',[char(10),strjoin(infile,char(10))])
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function outfile=export_acc(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      %parse optional arguments
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('format','ACC1B', @(i) ischar(i));
      p.addParameter('rm_ltb',true,    @(i) islogical(i));
      p.parse(varargin{:});
      %branch on output format
      switch p.Results.format
        case 'ACC1B'
          outdir_default=fullfile('YY','acc','asc');
          filename='ACC1B_YYYY-MM-DD_SAT_VERSION.asc';
          sat_name='GRACE SAT';
        case 'msodp'
          outdir_default=fullfile('YY','acc','dat');
          filename='grace.YYYY-MM-DD_SAT_VERSION.acc.pre';
          sat_name='GRACESAT';
        otherwise
          error([mfilename,': cannot handle format ''',p.Results.format,'''.'])
      end
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      %retrieve relevant parameters
      sats =product.mdget('sats');
      outdir=product.mdget('export_dir','default',outdir_default);
      version=num2str(product.mdget('version'));
      %gather list of daily data files
      [~,startlist,stoplist]=product.file('data',varargin{:},...
        'start',obj.start,...
        'stop', obj.stop...
      );
      %outputs
      outfile=cell(size(startlist)*2);
      %loop over the satellites
      for s=1:numel(sats)
        %save the data, if not empty
        out=obj.data_get_scalar(product.dataname.append_field_leaf(sats{s}));
        %loop over all dates
        for i=1:numel(startlist)
          idx=(s-1)*numel(startlist)+i;
          %build input data filename
          outfile{idx}=str.rep(fullfile(outdir,filename),...
            'YYYY',datestr(startlist(i),'yyyy'),...
            'YY',  datestr(startlist(i),'yy'  ),...
            'MM',  datestr(startlist(i),'mm'  ),...
            'DD',  datestr(startlist(i),'dd'  ),...
            'SAT', sats{s},...
            'VERSION',version ...
          );
          %trim only today
          tmp=out.trim(startlist(i),stoplist(i));
          if ~isempty(tmp)
            %remove long-term biases, if requested
            if p.Results.rm_ltb
              tmp=csr.ltb_apply(product,tmp,sats{s},'-all');
            end
            %save the data in ACC1B format, as handled by simpletimeseries.export
            tmp.export(...
              outfile{idx},p.Results.format,'sat_name',strrep(sat_name,'SAT',sats{s}),...
              varargin{:}...
            );
          end
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function obj=compute_calmod(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      %sanity
      assert(product.nr_sources==2,...
        ['number of sources in product ',product.str,...
        ' is expected to be 2, not ',num2str(product.nr_sources),'.'])
      %get sources metadata
      calparp=obj.product_get(product.sources(1));
      l1baccp=obj.product_get(product.sources(2));
      %retrieve relevant parameters
      param_col =calparp.mdget('param_col');
      coords    =calparp.mdget('coords');
      sats      =product.mdget('sats');
       t_scale  =1/str.num(product.mdget('t' ,'default',1)); %str.num also handles 'inf'
      t2_scale  =1/str.num(product.mdget('t2','default',1)); %str.num also handles 'inf'

      for s=1:numel(sats)
        %gather quantities, l1b acc data only contains the satellite names (located at the leafs of cal par data)
        acc=obj.data_get_scalar(l1baccp.dataname.set_field_path(sats(s)));
        %handle exceptions (also deals with non-existing data)
        if ~isa(acc,'simpletimeseries')
          %patch nan calibration model
          calmod=simpletimeseries(...
            [obj.start;obj.stop],...
            nan(2,numel(coords)),...
            'descriptor',['calibration model ',product.str,', GRACE-',sats{s},' (empty)']...
          );
          str.say('Skipping  the ',calmod.descriptor)
          continue
        end
        %define field names (this is arbitrary but must have length 3)
        fields={'bias','linear','quad'};
        %gather calibratiom parameters
        for c=1:numel(coords)
          %retreive data
          cal.(coords{c})=struct(...
fields{1},obj.data_get_scalar(calparp.dataname.set_field_path([product.dataname.field_path,{['AC0',coords{c}    ]},sats(s)])),...
fields{2},obj.data_get_scalar(calparp.dataname.set_field_path([product.dataname.field_path,{['AC0',coords{c},'D']},sats(s)])),...
fields{3},obj.data_get_scalar(calparp.dataname.set_field_path([product.dataname.field_path,{['AC0',coords{c},'Q']},sats(s)])) ...
          );
        end
        %init models container
        calmod=simpletimeseries(acc.t,zeros(acc.length,numel(coords))).copy_metadata(acc);
        calmod.descriptor=['calibration model ',product.str,', GRACE-',sats{s}];
        str.say('Computing the ',calmod.descriptor)
        for c=1:numel(coords)
          obj.log('@','iter','product',product,'sat',sats{s},'coord',coords{c})
          %loop over all fields
          for f=1:numel(fields)
            %interpolate to time domain of measurements
            cal.(coords{c}).(fields{f})=cal.(coords{c}).(fields{f}).interp(acc.t);
            obj.log('@','iter-interp',...
              'coord',coords{c},...
              [fields{f},' size'],cal.(coords{c}).(fields{f}).size,...
              [fields{f},' gaps'],cal.(coords{c}).(fields{f}).nr_gaps...
            )
            %sanity
            assert(~isempty(cal.(coords{c}).(fields{f})),[mfilename,...
              ': ',(fields{f}),' data is not available to perform this operation.'])
            %remove transient calpars, this is indicative of a gap (interpolation is done blindly over any gap length)
            cal.(coords{c}).(fields{f})=cal.(coords{c}).(fields{f}).mask_and([...
                   cal.(coords{c}).(fields{f}).mask(1);...
              diff(cal.(coords{c}).(fields{f}).y(:,1))==0]...
            );
            obj.log('@','iter-fixed ',...
              'coord',coords{c},...
              [fields{f},' size'],cal.(coords{c}).(fields{f}).size,...
              [fields{f},' gaps'],cal.(coords{c}).(fields{f}).nr_gaps...
            )
            %retrieve time domain (it is the same for all cal pars)
            if any(strcmp(product.dataname.field_path,'estim'))
              %start of day
              sod=dateshift(cal.(coords{c}).(fields{f}).t,'start','day');
              %second of arc
              soa=sod+seconds(cal.(coords{c}).(fields{f}).y(:,4));
              %time domain for the calibration model: units are days and zero epoch is the start of the arc
              t.(fields{f})=days(acc.t-soa);
            else
              %TODO: ? preciso arranjar isto, o start arc tem de ser definido algures (para aak e accatt)
              t.(fields{f})=days(acc.t-simpletimeseries.ToDateTime(cal.(coords{c}).(fields{f}).y(:,end),'modifiedjuliandate'));
            end
          end
          %paranoid sanity check
          good_idx=~isnan(t.(fields{1}));
          assert(...
            ~any(t.(fields{1})(good_idx)~=t.(fields{2})(good_idx)) && ...
            ~any(t.(fields{1})(good_idx)~=t.(fields{3})(good_idx)),...
            'calibration time domain inconsistent between parameters, debug needed!')
          %build calibration model
          calmod=calmod.set_cols(c,...
            cal.(coords{c}).(fields{1}).cols(param_col                                 )+...
            cal.(coords{c}).(fields{2}).cols(param_col).times(t.(fields{2})   * t_scale)+...
            cal.(coords{c}).(fields{3}).cols(param_col).times(t.(fields{3}).^2*t2_scale)...
          );
          %make debug plots
          if product.debug_plot
            font_size_args={...
              'plot_fontsize_title',24,...
              'plot_fontsize_axis',18, ...
              'plot_fontsize_label',20 ...
            };
            %patch 
            filename=product.file('plot',...
              'start',obj.start,...
              'stop',obj.stop,...
              'use_storage_period',false,...
              'timestamp',true,...
              'suffix',[sats{s},'.',coords{c}]...
            );
            assert(numel(filename)==1,['expecting only one filename, not ',num2str(numel(filename)),'.'])
            if  ~exist(filename{1},'file')
              figure('visible',product.metadata.plot_visible);
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
                cal.(coords{c}).(fields{f}).plot(...
                  'columns',param_col,...
                  'zeromean',true,...
                  'title',cal.(coords{c}).(fields{f}).labels{param_col})
                h_label=get(gca,'YLabel');
                set(h_label,'String',strrep(get(h_label,'String'),cal.(coords{c}).(fields{f}).labels{param_col},''))
                grid on
                product.enforce_plot(font_size_args{:})
              end
              saveas(gcf,filename{1})
              str.say('Created plot ',filename{1})
            else
              str.say('Skipped plot ',filename{1})
            end
          end
        end
        %propagate it
        obj=obj.data_set(product.dataname.append_field_leaf(sats{s}),calmod);
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function obj=estimate_poly_calmod(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      %sanity
      assert(product.nr_sources==2,...
        ['number of sources in product ',product.str,...
        ' is expected to be 1, not ',num2str(product.nr_sources),'.'])
      %get sources metadata
      modaccp=obj.product_get(product.sources(1));
      l1baccp=obj.product_get(product.sources(2));
      %retrieve relevant parameters
      coords          =product.mdget('coords');
      sats            =product.mdget('sats');
      arcs_per_day    =cell2mat(product.mdget('arcs_per_day'));
      arcs_poly_order =cell2mat(product.mdget('arcs_poly_order'));
      l1b_median_order=product.mdget('l1b_median_order');
      arc_scheme      =product.mdget('arc_scheme');
      atl_file        =product.mdget('arc_timeline_file');
      %check if the product satellites include the one in the current product
      sat_now=product.dataname.field_path;
      if ~any(strcmp(sats,sat_now))
        obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop,'sats',sats,'sat from product',sat_now)
        %nothing to do
        return
      end
      %ignore this warning
      warning_id='MATLAB:polyfit:RepeatedPointsOrRescale';
      %gather list of days to process
      if dateshift(obj.stop,'start','day')==obj.stop
        [t0,t1]=time.day_list(obj.start,obj.stop-days(1));
      else
        [t0,t1]=time.day_list(obj.start,obj.stop);
      end
      switch arc_scheme
      case 'blind'
        startstop=[t0,t1];
      case 'csr'
        startstop=[];
        for i=1:numel(t0)
          startstop=[startstop;...
            cells.c2m(...
              csr.arctimeline_get(...
                'arc_timeline_file',atl_file,...
                'year',  year(t0(i)),...
                'month',month(t0(i)),...
                'day',    day(t0(i)),...
                'out',{'start','stop'}...
              )...
            )...
          ]; %#ok<AGROW>
        end
        %need to be sure all epochs are within the object start/stop
        %(csr.arctimeline_get returns all entries which satisfy the logical union of the inputs)
        startstop(startstop(:,1)<obj.start,1)=obj.start;
        startstop(startstop(:,2)>obj.stop ,2)=obj.stop;
      otherwise
        error(['Cannot handle ''arc_scheme'' with value ''',arc_scheme,'''.'])
      end
      for d=1:size(startstop,1)
        obj.log('@','startstop feedback',...
          ['start(',num2str(d),')'],startstop(d,1),...
          ['stop(',num2str(d),')'],startstop(d,2))
      end
      %gather quantities
      l1b=obj.data_get_scalar(l1baccp.dataname.set_field_path(sat_now));
      mod=obj.data_get_scalar(modaccp.dataname.set_field_path(sat_now));
      %handle exceptions (also deals with non-existing data)
      if ~isa(l1b,'simpletimeseries')
        %patch nan calibration model
        calmod=simpletimeseries(...
          [obj.start;obj.stop],...
          nan(2,numel(coords)),...
          'descriptor',['calibration model ',product.str,' (empty)']...
        );
        str.say('Skipping  the ',calmod.descriptor)
      else
        str.say('Computing the calibration model ',product)
        %loop over all coordinates
        for c=1:numel(coords)
          col_idx=find(cells.iscellstrfind({'X','Y','Z'},coords{c}));
          %loop over all relevant (sub-) arcs
          for d=1:size(startstop,1)
            %retrieve (sub-) arcs
            t=csr.subarc_startstop(startstop(d,:),arcs_per_day(c));
            %loop over the number of arcs per day
            for a=1:size(t,1)
              obj.log('@','iter',...
                'product',product,'coord',coords{c},...
                'arc',startstop(d,:),...
                'sub-arc',t(a,:),'sub-arc',a)
              %get trimmed time series
              arc_mod=mod.trim(t(a,1),t(a,2)).cols(c);
              arc_l1b=l1b.trim(t(a,1),t(a,2)).cols(c).median(l1b_median_order).interp(arc_mod.t);
              
              %get calibration curve to fit
              arc_cur=arc_mod+arc_l1b;
              %compute calibration model
              warning('off',warning_id)
              arc_calmod=arc_cur.scale(-1).polyfit(arcs_poly_order(c));
              warning('on',warning_id)
              
%               %get calibration curve to fit
%               arc_cur1=arc_mod-arc_l1b;
%               arc_cur2=arc_mod+arc_l1b;
%               %compute calibration model
%               warning('off',warning_msg)
%               [arc_calmod1,S1]=arc_cur1.polyfit(arcs_poly_order(c));
%               [arc_calmod2,S2]=arc_cur2.polyfit(arcs_poly_order(c));
%               warning('on',warning_msg)
%               %find which calibration curve works best
%               if S1.normr>S2.normr
%                 factor=-1;
%                 arc_calmod=arc_calmod2;
%               else
%                 factor=1;
%                 arc_calmod=arc_calmod1;
%               end
                            
              %create/append daily data
              if a==1; day_calmod=arc_calmod;
              else     day_calmod=day_calmod.append(arc_calmod);
              end
            end
            %create/append complete time domain
            if d==1; coord_calmod=day_calmod;
            else     coord_calmod=coord_calmod.append(day_calmod);
            end
          end
          %create/glue all coordinates
          if c==1; calmod=coord_calmod;
          else     calmod=calmod.glue(coord_calmod);
          end          
        end
        calmod.descriptor=['calibration model ',product.str];
      end
      %propagate it
      obj=obj.data_set(product,calmod);
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function obj=import_acc_mod(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      %retrieve relevant parameters
      sats =product.mdget('sats');
      indir=product.mdget('import_dir');
      acc_version =product.mdget('acc_version' );
      gps_version =product.mdget('gps_version' );
      grav_version=product.mdget('grav_version');
      %gather list of daily data files
      [~,timestamplist]=product.file('data',varargin{:},'start',obj.start,'stop', obj.stop);
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
        obj.log('@','iter','sat',sats{s},'files to load',infile)
        %load (and save the data in mat format, as handled by simpletimeseries.import)
        obj=obj.data_set(product.dataname.set_field_path(sats(s)),...
          simpletimeseries.import(infile,'cut24hrs',true)...
        );
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    %% atomic loading of calibration parameters
    function out=estim_translate_param(in)
      if ~isempty(strfind(in,'AC0'))
        out=in;
      elseif ~isempty(strfind(in,'AC'))
        out=[in(1:2),'0',in(3:4)];
      elseif ~isempty(strfind(in,'C0Y'))
        out=['AC',in(3:4),in(1)];
      elseif ~isempty(strfind(in,'CYD')) || ...
             ~isempty(strfind(in,'CYQ'))
        out=['AC0',in(3:4),in(1)];
      elseif ~isempty(strfind(in,'S-T')) || ...
             ~isempty(strfind(in,'C-T')) || ...
             ~isempty(strfind(in,'S-N')) || ...
             ~isempty(strfind(in,'C-N'))
        switch numel(in)
        case 4; out=[in(2),in(4),'0',in(1  )];
        case 5; out=[in(3),in(5),    in(1:2)];
        otherwise; error('singularity, debug needed!')
        end
      else
        out=in;
      end
    end
    function out=estim_filename(varargin)
      %parse optional arguments
      p=inputParser;p.KeepUnmatched=true;
      p.addParameter('estimdir', '', @(i) ischar(i)    && ~isempty(i));
      p.addParameter('arc',      [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i));
      p.parse(varargin{:});
      out=file.wildcard(fullfile(p.Results.estimdir,...
        ['R.common.arc',num2str(p.Results.arc),'.to.arc',num2str(p.Results.arc),'.*.estim']));
      assert(numel(out)==1,['Can only handle one estim file at a time, not ',num2str(numel(out)),':',10,strjoin(out(:),char(10))])
      out=out{1};
    end
    function out=estim_load(filename)
      if exist([filename,'.mat'],'file')
        load([filename,'.mat'])
      elseif strcmp(filename(end-3:end),'.mat')
        load(filename)
      else
        %init the outputs
        out=struct('A',[],'B',[]);
        %open the file
        fid=file.open(filename);
        %load the data
        dat=textscan(fid,'%5s%5s %f %f %f');
        %close the file
        fclose(fid);
        %parse the data
        for i=1:numel(dat{1})
          fn=csr.estim_translate_param(dat{2}{i});
          if ~isvarname(fn)
            fn=['f',strrep(fn,'-','_')];
          end
          out.(strrep(dat{1}{i},'GRC-','')).(fn)=dat{5}(i);
        end
        save([filename,'.mat'],'out')
      end
    end
    function T=arctimeline_load(filename)
      %init the outputs
      start=datetime;
      stop =datetime;
      %open the file
      fid=file.open(filename);
      %load the data
      % 1  2  3  4     5       6       7
      %18 12/17/03 52990 49700.0 36700.0
      dat=textscan(fid,'%f %f/%f/%f %f %f %f');
      %close the file
      fclose(fid);
      %table names
      arc=dat{1};
      month=dat{2};
      day=dat{3};
      year=2000+dat{4};
      date=datetime([year,month,day]);
      mjd=dat{5};
      assert(all(mjd==time.mjd(date)),['Discrepancy between date and MJD in the file ''',filename,'''.'])
      start=date+ seconds(dat{6});
      stop =start+seconds(dat{7});
      %create the table
      T=table(year,month,day,arc,start,stop);
    end
    function out=arctimeline_get(varargin)
      %parse optional arguments
      p=inputParser;
      p.addParameter('arc_timeline_file',csr.arc_timeline_file, @(i) ischar(i) && ~isempty(i));
      p.addParameter('year', [], @(i) iscell(i) || isnumeric(i)  ||  isempty(i));
      p.addParameter('month',[], @(i) iscell(i) || isnumeric(i)  ||  isempty(i));
      p.addParameter('arc',  [], @(i) iscell(i) || isnumeric(i)  ||  isempty(i));
      p.addParameter('day',  [], @(i) iscell(i) || isnumeric(i)  ||  isempty(i));
      p.addParameter('start',[], @(i) iscell(i) || isdatetime(i) ||  isempty(i));
      p.addParameter('stop', [], @(i) iscell(i) || isdatetime(i) ||  isempty(i));
      p.addParameter('out',  {}, @(i) iscell(i)                  && ~isempty(i));
      p.parse(varargin{:});
      %load the arc timeline
      arctimeline=csr.arctimeline_load(p.Results.arc_timeline_file);
      %may need to change parameters, so need to get them out of the parser object
      par=p.Results;
      %init logic key
      key=true(height(arctimeline),1);
      %build logic key
      for i={'year','month','arc','day','start','stop'}
        if ~isempty(par.(i{1}))
          %make cells into matrices
          par.(i{1})=cells.c2m(par.(i{1}));
          key_now=false(size(key));
          for j=1:numel(par.(i{1}))
            key_now=key_now| (arctimeline.(i{1})==par.(i{1})(j));
          end
          key=key&key_now;
        end
      end
      %warn if nothing is returned
      if ~any(key)
        %try to patch things up
        if ~isempty(par.year) && ~isempty(par.month) && ~isempty(par.day)
          out=cell(size(par.out));
          for i=1:numel(par.out)
            switch par.out{i}
            case 'start'; out{i}=datetime(par.year,par.month,par.day,0,0,0);
            case 'stop';  out{i}=datetime(par.year,par.month,par.day,23,59,59);
            case 'arc';   out{i}=[];
            end
          end
          str.say('Patched default arc for',str.show(varargin))
        else
          out=[];
          str.say('No arc found for',str.show(varargin))
        end
        return
      end
      relevant=arctimeline(key,par.out);
      %reduce from table to fields requested in par.out
      out=cell(size(par.out));
      for i=1:numel(par.out)
        out{i}=relevant.(par.out{i});
        %special handling
        switch par.out{i}
        case 'stop'
          %shave 1 second from stop time
          out{i}=out{i}-seconds(1);
        end
      end
    end
    function t=subarc_startstop(startstop,nr_sub_arcs)
      %short-circuit trivial calls
      if nr_sub_arcs==1
        t=startstop;
      end
      %vector mode
      if size(startstop,1)>1
        assert(size(startstop,2)==2,'In vector mode, can only handle input ''startstop'' with two columns.')
        t=[];
        for i=1:size(startstop,1)
          t=[t;csr.subarc_startstop(startstop(i,:),nr_sub_arcs)]; %#ok<AGROW>
        end
        return
      end
      %some sanity
      assert(numel(startstop)==2,'In scalar mode, can only handle input ''startstop'' with two entries.')
      %there are two instances where the arcs overshoot mid-night, fix those
      if dateshift(startstop(1),'end','day')==dateshift(startstop(2),'start','day') && ...
        startstop(2)-dateshift(startstop(1),'end','day')<minutes(10)
        startstop(2)=dateshift(startstop(1),'end','day')-seconds(1);
      end
      %more sanity, MORE!
      assert(...
        dateshift(startstop(1),'start','day') ==  dateshift(startstop(2),'start','day') && ...
        dateshift(startstop(1),'end',  'day') ==  dateshift(startstop(2),'end',  'day'),...
        'Cannot handle input ''startstop'' that refers to different days')
      %init loop parameters
      t=[];
      %trivial call
      if startstop(1)==startstop(2)
        t=[];
        return
      end
      %determine sub-arc period
      periodicity=days(1)/nr_sub_arcs;
      %loop over all sub-arcs
      for i=1:nr_sub_arcs
        %get sub-arc starts/ends
        t_now=[startstop(1)+periodicity*(i-1),...
          min([startstop(1)+periodicity*(i  )-seconds(1),startstop(2)])...
        ];
        %singularity
        if diff(t_now)<=0
          return
        end
        %append
        t=[t;t_now]; %#ok<AGROW>
      end
    end
    function arc=arc_timeseries(startstop,calpars)
      %get persistent vars 
      persistent prefix suffix units coords nr_sub_arcs
      if isempty(prefix)
        prefix='AC0';
        suffix={'','D','Q'};
        units={'m/s^2','m/s^2','m/s^2'}; %these units are wrong but they need to agree with each other to build the cal mod2
        coords={'x','y','z'};
        nr_sub_arcs=[1,8,1];
      end
      %get sats
      sats=fieldnames(calpars);
      %loop over all sats
      for s=1:numel(sats)
        %loop over all components
        for c=1:numel(coords)
          %loop over all suffixes
          for i=1:numel(suffix)
            %build field name
            field=[prefix,upper(coords{c}),suffix{i}];
            %build ordinate: get sub-arc start/stop
            t=csr.subarc_startstop(startstop,nr_sub_arcs(c));
            %loop over all sub-arcs
            for day_loc=1:size(t,1)
              %retrieve (sub-) arc data
              if nr_sub_arcs(c)==1
                y=ones(2,1)*calpars.(sats{s}).(field);
              else
                y=ones(2,1)*calpars.(sats{s}).([field,num2str(day_loc)]);
              end
              %build timeseries object
              ts=simpletimeseries(transpose(t(day_loc,:)),y,...
                'format','datetime',...
                'labels',{field},...
                'units',units(i),...
                'timesystem','gps',...
                'descriptor',field...
              );
              %init/append output
              if day_loc==1
                arc.(sats{s}).(field)=ts;
              else
                arc.(sats{s}).(field)=arc.(sats{s}).(field).append(ts);
              end
            end
          end
        end
      end
    end
    function calmod=arc_build(product,varargin)
      %parse optional arguments
      p=inputParser;p.KeepUnmatched=true;
      p.addParameter('estimdir', '', @(i) ischar(i)    && ~isempty(i));
      p.addParameter('year',     [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i));
      p.addParameter('month',    [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i));
      p.addParameter('arc',      [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i));
      p.parse(varargin{:});
      %get names of parameters and levels
      coords    =product.mdget('coords');
      units=cell(size(coords));units(:)={'m/s^2'};
      sats      =product.mdget('sats');
      atl_file  =product.mdget('arc_timeline_file');
      %get start/stop of current arc
      startstop=cells.c2m(...
          csr.arctimeline_get(...
            'arc_timeline_file',atl_file,...
            'year', p.Results.year,...
            'month',p.Results.month,...
            'arc',  p.Results.arc,...
            'out',{'start','stop'}...
        )...
      );
      %load calpars
      calpars=csr.estim_load(csr.estim_filename(varargin{:}));
      %build calibration parameters time series
      arc=csr.arc_timeseries(startstop,calpars);
      %define field names (this is arbitrary but must have length 3)
      fields={'bias','linear','quad'};
      %set calibration model time domain
      t=startstop(1):seconds(1):startstop(2);
      %set the time domain in units compatible with the calibration parameters:
      %units are days and zero epoch is the start of the arc
      tc=days(t-startstop(1));
      %loop over all satellites
      for s=1:numel(sats)
        %gather calibratiom parameters
        for c=1:numel(coords)
          %retreive data
          cal.(coords{c})=struct(...
            fields{1},arc.(sats{s}).(['AC0',coords{c}    ]).interp(t),...
            fields{2},arc.(sats{s}).(['AC0',coords{c},'D']).interp(t),...
            fields{3},arc.(sats{s}).(['AC0',coords{c},'Q']).interp(t) ...
          );
          %loop over all fields
          for f=1:numel(fields)
            %remove transient calpars, this is indicative of a gap (interpolation is done blindly over any gap length)
            cal.(coords{c}).(fields{f})=cal.(coords{c}).(fields{f}).mask_and([...
                   cal.(coords{c}).(fields{f}).mask(1);...
              diff(cal.(coords{c}).(fields{f}).y(:,1))==0]...
            );
          end
        end
        %init calibration model time series
        calmod.(sats{s})=simpletimeseries(t,zeros(numel(t),numel(coords)),...
          'format','datetime',...
          'labels',coords,...
          'units',units,...
          'timesystem','gps',...
          'descriptor',['calibration model ',product.str,', GRACE-',sats{s}]...
        );
        %build calibration model
        for c=1:numel(coords)
          calmod.(sats{s})=calmod.(sats{s}).set_cols(c,...
            cal.(coords{c}).(fields{1}).cols(1             )+...
            cal.(coords{c}).(fields{2}).cols(1).times(tc   )+...
            cal.(coords{c}).(fields{3}).cols(1).times(tc.^2)...
          );
        end
        %apply long-term biases
        calmod.(sats{s})=csr.ltb_apply(product,calmod.(sats{s}),sats{s},'all');
      end
    end
    function out=arcs_build(product,varargin)
      %parse optional arguments
      p=inputParser;
      p.addParameter('estimdir',{}, @(i)                 ischar(i) && ~isempty(i));
      p.addParameter('year',    {}, @(i) iscell(i) || isnumeric(i) && ~isempty(i));
      p.addParameter('month',   {}, @(i) iscell(i) || isnumeric(i) && ~isempty(i));
      p.addParameter('day',     {}, @(i) iscell(i) || isnumeric(i) || isempty(i));
      p.parse(varargin{:});
      %translate days into arcs (estim files are arc-wise)
      if ~isempty(p.Results.day)
        %get arcs of current day
        arcs=cells.c2m(csr.arctimeline_get(...
          'arc_timeline_file',product.mdget('arc_timeline_file'),...
          'year', p.Results.year,...
          'month',p.Results.month,...
          'day',  p.Results.day,...
          'out',{'arc'}...
        ));
      else
        arcs=[];
      end
      %add arcs, which now default to the arcs of the given day (if given)
      p.addParameter('arc',arcs, @(i) iscell(i) || isnumeric(i) || isempty(i));
      p.parse(varargin{:});
      %handle scalar inputs
      n=max([...
        numel(p.Results.year),...
        numel(p.Results.month),...
        numel(p.Results.day),...
        numel(p.Results.arc)...
      ]);
      for i={'year','month','day','arc'}
        par.(i{1})=cells.deal(p.Results.(i{1}),[n,1]);
      end
      %sanity: all inputs must have the same length
      assert(all(cellfun(@(i) numel(par.(i)),{'year','month','day','arc'})==n),...
        'All time parameters must have the same length, debug needed!')
      %loop over all inputs
      s.msg='Assembling arcs';s.n=n;
      for i=1:n
        tmp=csr.arc_build(product,...
          'estimdir', p.Results.estimdir,...
          'year',     par.year{i},...
          'month',    par.month{i},...
          'day',      par.day{i},...
          'arc',      par.arc{i}...
        );
        if i==1
          out=tmp;
        else
          out=structs.objmethod2('append',out,tmp);
        end
        s=time.progress(s);
      end
    end
    %% legacy (needs checking)
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
      case {'debug_plots_import_calpar','debug_plots_calpar'}
        plot_dir=fullfile('plot',lower(mode),csr.gitversion);
        if ~isempty(dir(plot_dir))
          [status,result]=system(['rm -fvr "',plot_dir,'"']);
          assert(status==0,['error removing ',mode,': ',result])
        end
      otherwise
        dataproduct(mode,'metadata_dir',obj.metadata_dir).rm_data(varargin{:});
      end
    end
  end
end