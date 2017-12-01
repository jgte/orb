classdef csr
  %static
  properties(Constant)
    arc_timeline_file=[getenv('HOME'),'/data/csr/GraceAccCal/arctimeline.GraceAccCal'];
    estim_dir_root=[getenv('HOME'),'/data/csr/EstimData/RL05'];
    estimdir_file=[csr.estim_dir_root,'/EstimDirs_RL05'];
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
      case 'calpar';      name='grace.calpar.csr';
      case 'import_l1b';  name='grace.acc.l1b.csr';
      case 'calmod';      name='grace.acc.calmod.csr';
      case 'mod';         name='grace.acc.mod.csr';
      case 'nrtdm';       name='grace.acc.mod.nrtdm';
      case 'calacc';      name='grace.acc.cal.csr';
      case 'poly_calmod'; name='grace.acc.calmod.poly0';
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
    function out=ltb_args(i)
      out=varargs({
        'field',    'all', @(i) ischar(i)    &&                ~isempty(i);...
        'ltb_scale',    1, @(i) isnumeric(i) && isscalar(i) && ~isempty(i);...
        'calpar_scale', 1, @(i) isnumeric(i) && isscalar(i) && ~isempty(i);...
        'sat',        'A', @(i) ischar(i)                   && ~isempty(i);...
        'bias_files', {...
          [getenv('HOME'),'/data/csr/corral-tacc/input/bsA2003'],...
          [getenv('HOME'),'/data/csr/corral-tacc/input/bsB2003']...
        },                 @(i) iscellstr(i) && numel(i)==2; ...
      });
      if exist('i','var') && ~isempty(i)
        out=out.(i);
      end
    end
    function out=ltb(varargin)
      v=varargs.wrap('sources',{{...
        't',    datetime('now'),@(i) isdatetime(i) && ~isempty(i);...
      },csr.ltb_args},varargin{:});
      %resolve bias files
      switch v.sat
      case 'A'; bias_file=v.bias_files{1};
      case 'B'; bias_file=v.bias_files{2};
      otherwise; error(['Unknown sat ''',v.sat,'''.'])  
      end
      %translate field to data index
      switch upper(v.field)
      case {'X','AC0X'}
        lbt_idx=2;
      case {'Y','AC0Y1','AC0Y2','AC0Y3','AC0Y4','AC0Y5','AC0Y6','AC0Y7','AC0Y8'}
        lbt_idx=3;
      case {'Z','AC0Z'}
        lbt_idx=4;
      case {'ALL','-ALL'}
        %if mode is '-all', then "remove" the long-term bias
        if strcmpi(v.field,'-ALL');v.ltb_scale=-v.ltb_scale;end
        x=csr.ltb(v.varargin{:},'field','x');
        y=csr.ltb(v.varargin{:},'field','y');
        z=csr.ltb(v.varargin{:},'field','z')';
        out=x.glue(...
            y.glue(....
            z));
        out.descriptor=['LTB read from ',str.clean(bias_file,'file')];
        return
      otherwise
        out=[];
        return
      end
      %get the data
      ltb_data=csr.ltb_data(bias_file);
      %number of day since reference epoch (in the ltb file)
      d=time.mjd(v.t)-ltb_data(2,1);
      %create timeseries object
      out=simpletimeseries(...
        v.t,...
        v.ltb_scale*polyval(ltb_data(:,lbt_idx),d),...
        'format','datetime',...
        'timesystem','gps',...
        'labels',{v.field},...
        'units',{'m/s^2'},...
        'descriptor',[v.field,' LTB read from ',str.clean(bias_file,'file')]...
      );
    end
    function out=ltb_apply(varargin)
      v=varargs.wrap('sources',{{...
        'ts',    [], @(i) isa(i,'simpletimeseries') && ~isempty(i);...
      },csr.ltb_args},varargin{:});
      %get long-term biases
      ltb=csr.ltb('t',v.ts.t,v.varargin{:});
      %add long-term biases (unless ltb is empty, which means this field had no ltb)
      if isempty(ltb)
        out=v.ts;
      else
        out=v.ts.scale(v.calpar_scale)+ltb;
      end
    end
    %% importers
    function obj=import_calpar_slow(obj,product,varargin)
      %open log file
      csr.log
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'debugdate',   [], @(i) ischar(i)   || isempty(i);
      },csr.ltb_args,product.metadata},varargin{:});
      %load data
      for i=1:numel(v.levels)
        for j=1:numel(v.calpars)
          tmp=struct('A',[],'B',[]);
          for s=1:numel(v.sats)
            %read L1B data
            f=fullfile(v.import_dir,['gr',v.sats{s},'.',v.calpars{j},'.',v.levels{i},'.GraceAccCal']);
            tmp.(v.sats{s})=simpletimeseries.import(f,'cut24hrs',false);
            %apply long-term bias
            tmp.(v.sats{s})=tmp.(v.sats{s}).set_cols(v.param_col,...
              csr.ltb_apply(v.varargin{:},...
                'ts',tmp.(v.sats{s}).get_cols(v.param_col),...
                'sat',v.sats{s},...
                'field',v.calpars{s}...
              )...
            );
            %additional processing: add end of arcs
            switch v.levels{i}
            case {'aak','accatt'}
              error('Needs implementation')
            case 'estim'
              %get arc stars
              arc_starts=tmp.(v.sats{s}).t;
              %build arc ends (arc duration given explicitly)
              arc_ends=arc_starts+seconds(tmp.(v.sats{s}).y(:,v.arclen_col))-seconds(1);
              %patch missing arc durations
              idx=find(isnat(arc_ends));
              %report edge cases
              csr.report(obj.debug,idx,'Arcs without arc length',f,...
                {'arc start','arc duraction'},...
                {arc_starts,  tmp.(v.sats{s}).y(:,v.arclen_col)}...
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
            if ~isempty(strfind(v.calpars{j},'AC0Y'))
              %there are 8 segments per day
              periodicity=days(1)/8;
              %get day location for this parameter
              day_loc=str2double(v.calpars{j}(end));
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
            arc_start_y=tmp.(v.sats{s}).y;
              arc_end_y=tmp.(v.sats{s}).y;

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
            if ~isempty(v.debugdate)
              rep_date=datetime(v.debugdate);
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
              'labels',tmp.(v.sats{s}).labels,...
              'units',tmp.(v.sats{s}).y_units,...
              'timesystem',tmp.(v.sats{s}).timesystem,...
              'descriptor',tmp.(v.sats{s}).descriptor...
            );
            %build timeseries with arc ends
            arc_end_ts=simpletimeseries(arc_ends,arc_end_y,...
              'format','datetime',...
              'labels',    tmp.(v.sats{s}).labels,...
              'units',     tmp.(v.sats{s}).y_units,...
              'timesystem',tmp.(v.sats{s}).timesystem,...
              'descriptor',['end of arcs for ',tmp.(v.sats{s}).descriptor]...
            );

            %augment arc starts with arc ends (only new data)
            tmp.(v.sats{s})=arc_start_ts.augment(arc_end_ts,'old',true);

          end

          %propagate data to object
          for s=1:numel(v.sats)
            obj=obj.data_set(product.dataname.set_field_path([v.levels(i),v.calpars(j),v.sats(s)]),tmp.(v.sats{s}));     
          end
          %user feedback
          str.say(str.tablify([15,6,3,6],'loaded data for',v.levels{i},'and',v.calpars{j}))
        end
      end

      %merge cross-track accelerations together
      ac0y='AC0Y';
      field_part_list={'','D','Q'};
      %loop over all levels and sats
      for i=1:numel(v.levels)
        for s=1:numel(v.sats)
          for f=1:numel(field_part_list)
            %start with first field
            calpar=[ac0y,field_part_list{f},'1'];
            ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),calpar,v.sats(s)]));
            %rename the relevant object fields to remove the '1'
            ts_now.labels=strrep(ts_now.labels,calpar,[ac0y,field_part_list{f}]);
            ts_now.descriptor=strrep(ts_now.descriptor,calpar,[ac0y,field_part_list{f},'[1-8]']);
            %loop over all other calpars
            for fpl=2:8
              calpar=[ac0y,field_part_list{f},num2str(fpl)];
              ts_now=ts_now.augment(...
                obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),calpar,v.sats(s)])),...
                'quiet',true,...
                'old',true,...
                'skip_gaps',true...
              );
              %debug date report
              if ~isempty(v.debugdate)
                rep_date=datetime(v.debugdate);
                str.say('DEBUG DATE: merge AC0Y*:',v.levels{i},':',v.sats{s},':',calpar,' @ ',datestr(rep_date));
                idx=ts_now.idx(rep_date);
                ts_now.peek((idx-10):(idx+10));
              end
            end
            %save the data
            obj=obj.data_set(product.dataname.set_field_path([v.levels(i),{[ac0y,field_part_list{f}]},v.sats(s)]),ts_now);
            %user feedback
            str.say(str.tablify([29,5,3,6,3,7],'merged cross-track parameter',[ac0y,field_part_list{f}],...
              'for',v.levels{i},'and',['GRACE-',v.sats{s}]))
          end
        end
      end

      %add gaps
      for i=1:numel(v.levels)
        for j=1:numel(v.calpars_out)
          for s=1:numel(v.sats)
            ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]));
            %debug date report
            if ~isempty(v.debugdate)
              rep_date=datetime(v.debugdate);
              str.say('DEBUG DATE: w/out gaps:',v.levels{i},':',v.sats{s},':',v.calpars_out{j},' @ ',datestr(rep_date));
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
            if ~isempty(v.debugdate)
              rep_date=datetime(v.debugdate);
              str.say('DEBUG DATE: with gaps:',v.levels{i},':',v.sats{s},':',v.calpars_out{j},' @ ',datestr(rep_date));
              idx=ts_now.idx(rep_date);
              ts_now.peek((idx-10):(idx+10));
            end              
            %save
            obj=obj.data_set(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]),ts_now);  
          end
        end
      end

      %loop over all sat and level to check Job IDs agreement across all output calpars
      for i=1:numel(v.levels)
        for s=1:numel(v.sats)
          for j=1:numel(v.calpars_out)-1
            dn1=product.dataname.set_field_path([v.levels(i),v.calpars_out(j  ),v.sats(s)]);
            dn2=product.dataname.set_field_path([v.levels(i),v.calpars_out(j+1),v.sats(s)]);
            d1=obj.data_get_scalar(dn1);
            d2=obj.data_get_scalar(dn2);
            [~,i1,i2]=intersect(d1.t,d2.t);
            bad_idx=find(...
              d1.y(i1,v.jobid_col) ~= d2.y(i2,v.jobid_col) & ...
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
                msg{2*k+1}=str.tablify([30,6,20,12],dn1,idx,d1.t(idx),num2str(d1.y(idx,v.jobid_col),'%i'));
                idx=i2(bad_idx(k));
                msg{2*k+2}=str.tablify([30,6,20,12],dn2,idx,d2.t(idx),num2str(d2.y(idx,v.jobid_col),'%i'));
              end
              error([mfilename,':',strjoin(msg,'\n')])
            end
          end
        end
      end

      %loop over all sats, levels and calpars to:
      % - in case of estim: ensure that there are no arcs with lenghts longer than consecutive time stamps
      % - in case of aak and accatt: ensure that the t0 value is the same as the start of the arc
      for s=1:numel(v.sats)
        %loop over all required levels
        for i=1:numel(v.levels)
          switch v.levels{i}
          case 'estim'
            %this check ensures that there are no arcs with lenghts longer than consecutive time stamps
            for j=1:numel(v.calpars_out)
              %some calpars do not have t0
              if ~any(v.calpars_out{j}(end)=='DQ') || ~isempty(strfind(v.calpars_out{j},'Y'))
                str.say(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              str.say(str.tablify([8,32],'Checking',product.str))
              %save time series into dedicated var
              ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]));
              %forget about epochs that have been artificially inserted to represent gaps and end of arcs
              %the 1e-6 parcel is needed to avoid artificially-inserted gaps that have round-off errors
              idx1=find(diff(ts_now.t)>seconds(1)+1e-6);
              %get arc lenths
              al=ts_now.y(idx1,v.arclen_col);
              %get consecutive time difference
              dt=[seconds(diff(ts_now.t));0]; dt=dt(idx1);
              %find arcs that span over time stamps
              bad_idx=find(al-dt>2); %no abs here!
              %report if any such epochs have been found
              csr.report(obj.debug,bad_idx,'Ilegal arc length in the data',[v.levels{i},'.',v.calpars_out{j},'.',v.sats{s}],...
                {'global idx','arc init t','arc length','succ time diff','delta arc len'},...
                {idx1,ts_now.t(idx1),al,dt,al-dt}...
              ) %#ok<FNDSB>
            end
          case {'aak','accatt'}
            %this check ensures that the t0 value is the same as the start of the arc
            for j=1:numel(v.calpars_out)
              %the Y parameter was constructed from multitple parameters and some calpars do not have t0
              if ~any(v.calpars_out{j}(end)=='DQ') || ~isempty(strfind(v.calpars_out{j},'Y')) 
                str.say(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              str.say(str.tablify([8,32],'Checking',product.str))
              %save time series into dedicated var
              ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]));
              %forget about epochs that have been artificially inserted to represent forward steps
              idx1=find(diff(ts_now.t)>seconds(1));
              %get t0
              t0=simpletimeseries.utc2gps(datetime(ts_now.y(idx1,v.t0_col),'convertfrom','modifiedjuliandate'));
              %find arcs that have (much) t0 different than their first epoch
              bad_idx=find(...
                abs(ts_now.t(idx1)-t0)>seconds(1) & ...
                ts_now.mask(idx1) & ...                          %ignore gaps
                [true;diff(ts_now.y(idx1,v.jobid_col))~=0] ...     %ignore epochs inside the same arc
              );
              %report if any such epochs have been found
              csr.report(obj.debug,bad_idx,'Ilegal t0 in the data',[v.levels{i},'.',v.calpars_out{j},'.',v.sats{s}],...
                {'global idx','arc init time','t0','delta time'},...
                {idx1,ts_now.t(idx1),t0,ts_now.t(idx1)-t0}...
              ) %#ok<FNDSB>bo
            end
          end
        end
      end
    end
    function obj=import_calpar(obj,product,varargin)
      %open log file
      csr.log
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'debugdate',   [], @(i) ischar(i)   || isempty(i);
      },csr.ltb_args,product.metadata},varargin{:});
      %load data
      for i=1:numel(v.levels)
        for j=1:numel(v.calpars)
          tmp=struct('A',[],'B',[]);
          for s=1:numel(v.sats)
            %read L1B data
            f=fullfile(v.import_dir,['gr',v.sats{s},'.',v.calpars{j},'.',v.levels{i},'.GraceAccCal']);
            tmp.(v.sats{s})=simpletimeseries.import(f,'cut24hrs',false);
            %apply long-term bias
            tmp.(v.sats{s})=tmp.(v.sats{s}).set_cols(v.param_col,...
              csr.ltb_apply(v.varargin{:},...
                'ts',tmp.(v.sats{s}).get_cols(v.param_col),...
                'sat',v.sats{s},...
                'field',v.calpars{s}...
              )...
            );
            %additional processing: add end of arcs
            switch v.levels{i}
            case {'aak','accatt'}
              %get arc stars
              arc_starts=tmp.(v.sats{s}).t;
              %build arc ends
              arc_ends=[arc_starts(2:end);dateshift(arc_starts(end),'end','day')]-seconds(1);
%                 %arc ends are at maximum 24 hours after arc starts (only for those arcs starting at mid-night)
%                 fix_idx=arc_ends-arc_starts>days(1) & ...
%                   seconds(arc_starts-dateshift(arc_starts,'start','day'))<tmp.(v.sats{s}).t_tol;
%                 arc_ends(fix_idx)=arc_starts(fix_idx)+days(1)-seconds(1);
            case 'estim'
              %get arc stars
              arc_starts=tmp.(v.sats{s}).t;
              %build arc ends (arc duration given explicitly)
              arc_ends=arc_starts+seconds(tmp.(v.sats{s}).y(:,v.arclen_col))-seconds(1);
              %patch missing arc durations
              idx=find(isnat(arc_ends));
              %report edge cases
              csr.report(obj.debug,idx,'Arcs without arc length',f,...
                {'arc start','arc duraction'},...
                {arc_starts,  tmp.(v.sats{s}).y(:,v.arclen_col)}...
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
            if ~isempty(strfind(v.calpars{j},'AC0Y'))
              %there are 8 segments per day
              periodicity=days(1)/8;
              %get day location for this parameter
              day_loc=str2double(v.calpars{j}(end));
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
            arc_start_y=tmp.(v.sats{s}).y;
              arc_end_y=tmp.(v.sats{s}).y;

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
            if ~isempty(v.debugdate)
              rep_date=datetime(v.debugdate);
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
              'labels',tmp.(v.sats{s}).labels,...
              'units',tmp.(v.sats{s}).y_units,...
              'timesystem',tmp.(v.sats{s}).timesystem,...
              'descriptor',tmp.(v.sats{s}).descriptor...
            );
            %build timeseries with arc ends
            arc_end_ts=simpletimeseries(arc_ends,arc_end_y,...
              'format','datetime',...
              'labels',    tmp.(v.sats{s}).labels,...
              'units',     tmp.(v.sats{s}).y_units,...
              'timesystem',tmp.(v.sats{s}).timesystem,...
              'descriptor',['end of arcs for ',tmp.(v.sats{s}).descriptor]...
            );

            %augment arc starts with arc ends (only new data)
            tmp.(v.sats{s})=arc_start_ts.augment(arc_end_ts,'old',true);

          end

          %propagate data to object
          for s=1:numel(v.sats)
            obj=obj.data_set(product.dataname.set_field_path([v.levels(i),v.calpars(j),v.sats(s)]),tmp.(v.sats{s}));     
          end
          %user feedback
          str.say(str.tablify([15,6,3,6],'loaded data for',v.levels{i},'and',v.calpars{j}))
        end
      end

      %merge cross-track accelerations together
      ac0y='AC0Y';
      field_part_list={'','D','Q'};
      %loop over all levels and sats
      for i=1:numel(v.levels)
        for s=1:numel(v.sats)
          for f=1:numel(field_part_list)
            %start with first field
            calpar=[ac0y,field_part_list{f},'1'];
            ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),calpar,v.sats(s)]));
            %rename the relevant object fields to remove the '1'
            ts_now.labels=strrep(ts_now.labels,calpar,[ac0y,field_part_list{f}]);
            ts_now.descriptor=strrep(ts_now.descriptor,calpar,[ac0y,field_part_list{f},'[1-8]']);
            %loop over all other calpars
            for fpl=2:8
              calpar=[ac0y,field_part_list{f},num2str(fpl)];
              ts_now=ts_now.augment(...
                obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),calpar,v.sats(s)])),...
                'quiet',true,...
                'old',true,...
                'skip_gaps',true...
              );
              %debug date report
              if ~isempty(v.debugdate)
                rep_date=datetime(v.debugdate);
                str.say('DEBUG DATE: merge AC0Y*:',v.levels{i},':',v.sats{s},':',calpar,' @ ',datestr(rep_date));
                idx=ts_now.idx(rep_date);
                ts_now.peek((idx-10):(idx+10));
              end
            end
            %save the data
            obj=obj.data_set(product.dataname.set_field_path([v.levels(i),{[ac0y,field_part_list{f}]},v.sats(s)]),ts_now);
            %user feedback
            str.say(str.tablify([29,5,3,6,3,7],'merged cross-track parameter',[ac0y,field_part_list{f}],...
              'for',v.levels{i},'and',['GRACE-',v.sats{s}]))
          end
        end
      end

      %add gaps
      for i=1:numel(v.levels)
        for j=1:numel(v.calpars_out)
          for s=1:numel(v.sats)
            ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]));
            %debug date report
            if ~isempty(v.debugdate)
              rep_date=datetime(v.debugdate);
              str.say('DEBUG DATE: w/out gaps:',v.levels{i},':',v.sats{s},':',v.calpars_out{j},' @ ',datestr(rep_date));
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
            if ~isempty(v.debugdate)
              rep_date=datetime(v.debugdate);
              str.say('DEBUG DATE: with gaps:',v.levels{i},':',v.sats{s},':',v.calpars_out{j},' @ ',datestr(rep_date));
              idx=ts_now.idx(rep_date);
              ts_now.peek((idx-10):(idx+10));
            end              
            %save
            obj=obj.data_set(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]),ts_now);  
          end
        end
      end

      %loop over all sat and level to check Job IDs agreement across all output calpars
      for i=1:numel(v.levels)
        for s=1:numel(v.sats)
          for j=1:numel(v.calpars_out)-1
            dn1=product.dataname.set_field_path([v.levels(i),v.calpars_out(j  ),v.sats(s)]);
            dn2=product.dataname.set_field_path([v.levels(i),v.calpars_out(j+1),v.sats(s)]);
            d1=obj.data_get_scalar(dn1);
            d2=obj.data_get_scalar(dn2);
            [~,i1,i2]=intersect(d1.t,d2.t);
            bad_idx=find(...
              d1.y(i1,v.jobid_col) ~= d2.y(i2,v.jobid_col) & ...
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
                msg{2*k+1}=str.tablify([30,6,20,12],dn1,idx,d1.t(idx),num2str(d1.y(idx,v.jobid_col),'%i'));
                idx=i2(bad_idx(k));
                msg{2*k+2}=str.tablify([30,6,20,12],dn2,idx,d2.t(idx),num2str(d2.y(idx,v.jobid_col),'%i'));
              end
              error([mfilename,':',strjoin(msg,'\n')])
            end
          end
        end
      end

      %loop over all sats, levels and calpars to:
      % - in case of estim: ensure that there are no arcs with lenghts longer than consecutive time stamps
      % - in case of aak and accatt: ensure that the t0 value is the same as the start of the arc
      for s=1:numel(v.sats)
        %loop over all required levels
        for i=1:numel(v.levels)
          switch v.levels{i}
          case 'estim'
            %this check ensures that there are no arcs with lenghts longer than consecutive time stamps
            for j=1:numel(v.calpars_out)
              %some calpars do not have t0
              if ~any(v.calpars_out{j}(end)=='DQ') || ~isempty(strfind(v.calpars_out{j},'Y'))
                str.say(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              str.say(str.tablify([8,32],'Checking',product.str))
              %save time series into dedicated var
              ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]));
              %forget about epochs that have been artificially inserted to represent gaps and end of arcs
              %the 1e-6 parcel is needed to avoid artificially-inserted gaps that have round-off errors
              idx1=find(diff(ts_now.t)>seconds(1)+1e-6);
              %get arc lenths
              al=ts_now.y(idx1,v.arclen_col);
              %get consecutive time difference
              dt=[seconds(diff(ts_now.t));0]; dt=dt(idx1);
              %find arcs that span over time stamps
              bad_idx=find(al-dt>2); %no abs here!
              %report if any such epochs have been found
              csr.report(obj.debug,bad_idx,'Ilegal arc length in the data',[v.levels{i},'.',v.calpars_out{j},'.',v.sats{s}],...
                {'global idx','arc init t','arc length','succ time diff','delta arc len'},...
                {idx1,ts_now.t(idx1),al,dt,al-dt}...
              ) %#ok<FNDSB>
            end
          case {'aak','accatt'}
            %this check ensures that the t0 value is the same as the start of the arc
            for j=1:numel(v.calpars_out)
              %the Y parameter was constructed from multitple parameters and some calpars do not have t0
              if ~any(v.calpars_out{j}(end)=='DQ') || ~isempty(strfind(v.calpars_out{j},'Y')) 
                str.say(str.tablify([8,32],'Skipping',product.str))
                continue
              end
              str.say(str.tablify([8,32],'Checking',product.str))
              %save time series into dedicated var
              ts_now=obj.data_get_scalar(product.dataname.set_field_path([v.levels(i),v.calpars_out(j),v.sats(s)]));
              %forget about epochs that have been artificially inserted to represent forward steps
              idx1=find(diff(ts_now.t)>seconds(1));
              %get t0
              t0=simpletimeseries.utc2gps(datetime(ts_now.y(idx1,v.t0_col),'convertfrom','modifiedjuliandate'));
              %find arcs that have (much) t0 different than their first epoch
              bad_idx=find(...
                abs(ts_now.t(idx1)-t0)>seconds(1) & ...
                ts_now.mask(idx1) & ...                          %ignore gaps
                [true;diff(ts_now.y(idx1,v.jobid_col))~=0] ...     %ignore epochs inside the same arc
              );
              %report if any such epochs have been found
              csr.report(obj.debug,bad_idx,'Ilegal t0 in the data',[v.levels{i},'.',v.calpars_out{j},'.',v.sats{s}],...
                {'global idx','arc init time','t0','delta time'},...
                {idx1,ts_now.t(idx1),t0,ts_now.t(idx1)-t0}...
              ) %#ok<FNDSB>bo
            end
          end
        end
      end
    end
    function obj=import_l1b(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'import_subdir','',@(i) ischar(i);...
        'filename',     '',@(i) ischar(i);...
        'cat_command',  '',@(i) ischar(i);...
      },product.metadata},varargin{:});
      %branch on input format to set format-dependent parameters
      switch upper(v.format)
      case 'ACC-ASCII'
        v.import_subdir=fullfile('YY','acc','asc');
        v.filename='ACC1B_YYYY-MM-DD_SAT_VERSION.asc';
        v.format='ACC1B';
      case 'ACC-MSODP'
        v.import_subdir=fullfile('YY','acc','dat');
        v.filename='grace.YYYY-MM-DD_SAT_VERSION.acc.pre';
        v.format=lower(v.format);
      case {'ACC1B','AHK1B','GNV1B','KBR1B','MAS1B','SCA1B','THR1B','CLK1B','GPS1B','IHK1B','MAG1B','TIM1B','TNK1B','USO1B','VSL1B'}
        v.cat_command=[fullfile(getenv('HOME'),'data','grace','cat-l1b.sh'),' YYYYMMDD ',upper(v.format),' SAT'];
      otherwise
        error([mfilename,': cannot handle data ''',v.format,'''.'])
      end
      %gather list of days
      [~,startlist]=product.file('data',varargin{:},...
        'start',obj.start,...
        'stop', obj.stop...
      );
      %loop over the satellites
      for s=1:numel(v.sats)
        %branch on binary (requires dedicated cat tool) or ascii data
        if isempty(v.cat_command)
          infile=cell(size(startlist));
          %loop over all dates
          for i=1:numel(startlist)
            %build input data filename
            infile{i}=str.rep(fullfile(v.import_dir,v.import_subdir,v.filename),...
              'YYYY',datestr(startlist(i),'yyyy'),...
              'YY',  datestr(startlist(i),'yy'  ),...
              'MM',  datestr(startlist(i),'mm'  ),...
              'DD',  datestr(startlist(i),'dd'  ),...
              'SAT', v.sats{s},...
              'VERSION',v.version ...
            );
          end
          %load the data, as handled by simpletimeseries.import
          out=simpletimeseries.import(infile,'cut24hrs',false,'format',v.format);
        else
          %loop over all dates
          for i=1:numel(startlist)
            %define temporary file
            infile{i}=fullfile(filesep,'tmp',...
              ['GRACE-',v.sats{s},'.',datestr(startlist(i),'yyyy-mm-dd'),'.',upper(v.format),'.',str.rand(6),'.tmp']...
            );
            %define the cat tool command
            com=[str.rep(v.cat_command,...
              'YYYY',datestr(startlist(i),'yyyy'),...
              'YY',  datestr(startlist(i),'yy'  ),...
              'MM',  datestr(startlist(i),'mm'  ),...
              'DD',  datestr(startlist(i),'dd'  ),...
              'SAT', v.sats{s}...
            ),' > ',infile{i}];
            %call cat tool
            [status,result]=system(com);
            %make sure that worked
            assert(status==0,['The following command failed with status ',num2str(status),':',10,...
              com,10,...
              'with the following screen output:',10,...
              result...
            ])
          end
          %load the data, as handled by simpletimeseries.import
          out=simpletimeseries.import(infile,'cut24hrs',false,'format',v.format,'save_mat',false);
          %clean up the tmp dir
          for i=1:numel(infile); delete(infile{i}); end
        end
        %save the data, if not empty
        if ~isempty(out)
          obj=obj.data_set(product.dataname.set_field_path(v.sats(s)),out);
        else
          str.say('No data in file(s):',[char(10),strjoin(infile,char(10))])
        end
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function obj=import_acc_mod(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{product.metadata},varargin{:});
      %gather list of daily data files
      [~,timestamplist]=product.file('data',varargin{:},'start',obj.start,'stop', obj.stop);
      %loop over the satellites
      for s=1:numel(v.sats)
        infile=cell(size(timestamplist));
        %loop over all dates
        for i=1:numel(timestamplist)
          %build input data filename
          infile{i}=fullfile(v.import_dir,datestr(timestamplist(i),'yy'),datestr(timestamplist(i),'mm'),'gps_orb_l',...
            ['grc',v.sats{s},'_gps_orb_',datestr(timestamplist(i),'yyyy-mm-dd'),...
            '_RL',v.acc_version,'_GPSRL',v.gps_version,'_RL',v.grav_version,'.*.acc']...
          );
        end
        obj.log('@','iter','sat',v.sats{s},'files to load',infile)
        %load (and save the data in mat format, as handled by simpletimeseries.import)
        obj=obj.data_set(product.dataname.set_field_path(v.sats(s)),...
          simpletimeseries.import(infile,'cut24hrs',true)...
        );
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    %% exporters
    function outfile=export_acc(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'format', 'ACC1B', @(i) ischar(i);...
        'rm_ltb', true,    @(i) islogical(i);...
      },csr.ltb_args,product.metadata},varargin{:});
      %branch on output format
      switch v.format
        case 'ACC1B'
          v.export_subdir=fullfile('YY','acc','asc');
          filename='ACC1B_YYYY-MM-DD_SAT_VERSION.asc';
          sat_name='GRACE SAT';
        case 'msodp'
          v.export_subdir=fullfile('YY','acc','dat');
          filename='grace.YYYY-MM-DD_SAT_VERSION.acc.pre';
          sat_name='GRACESAT';
        otherwise
          error([mfilename,': cannot handle format ''',v.format,'''.'])
      end
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      %gather list of daily data files
      [~,startlist,stoplist]=product.file('data',varargin{:},...
        'start',obj.start,...
        'stop', obj.stop...
      );
      %outputs
      outfile=cell(size(startlist)*2);
      %loop over the satellites
      for s=1:numel(v.sats)
        %save the data, if not empty
        out=obj.data_get_scalar(product.dataname.append_field_leaf(v.sats{s}));
        %loop over all dates
        for i=1:numel(startlist)
          idx=(s-1)*numel(startlist)+i;
          %build output data filename
          outfile{idx}=str.rep(fullfile(v.export_dir,v.export_subdir,filename),...
            'YYYY',datestr(startlist(i),'yyyy'),...
            'YY',  datestr(startlist(i),'yy'  ),...
            'MM',  datestr(startlist(i),'mm'  ),...
            'DD',  datestr(startlist(i),'dd'  ),...
            'SAT', v.sats{s},...
            'VERSION',v.version ...
          );
          %trim only today
          tmp=out.trim(startlist(i),stoplist(i));
          if ~isempty(tmp)
            %remove long-term biases, if requested
            if v.rm_ltb
              tmp=csr.ltb_apply(v.varargin{:},'ts',tmp,'sat',v.sats{s},'field','-all');
            end
            %save the data in ACC1B format, as handled by simpletimeseries.export
            tmp.export(...
              outfile{idx},v.format,'sat_name',strrep(sat_name,'SAT',v.sats{s}),...
              varargin{:}...
            );
          end
        end
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
    function out=estim_dir_default(varargin)
      persistent estimdata
      %parse optional arguments
      p=inputParser;p.KeepUnmatched=true;
      p.addParameter('estim_dir',     '',  @(i) ischar(i)                   || isempty(i));
      p.addParameter('year',          [],  @(i) isnumeric(i) && isscalar(i) || isempty(i));
      p.addParameter('month',         [],  @(i) isnumeric(i) && isscalar(i) || isempty(i));
      p.addParameter('release',   'RL05',  @(i) ischar(i));
      p.addParameter('estim_dir_root',csr.estim_dir_root,@(i) ischar(i));
      p.parse(varargin{:});
      %honour given estim_dir
      if ~isempty(p.Results.estim_dir)
        out=p.Results.estim_dir;
      else
        assert(~isempty(p.Results.year) && ~isempty(p.Results.month),['If argument ''estim_dir'' is empty, ',...
          'then need both arguments ''year'' and ''month.'])
        if isempty(estimdata)
          %load estimdir file
          fid=file.open(csr.estimdir_file);
%       0 2005 February  200502_180.GEO /corral-tacc/utexas/csr/byaa705/grace/grav/RL05_05-02/a/iter/ ./solution ./postfit ./L2 53402 53430
          estimdata=textscan(fid,'%1d %4d %s %s %s %s %s %s %d %d','MultipleDelimsAsOne',true);
          fclose(fid);
        end
        %get solution index
        date_particle=[num2str(p.Results.year-2000,'%02d'),'-',num2str(p.Results.month,'%02d')];
        idx=cells.iscellstrempty(estimdata{5},[p.Results.release, '_',date_particle]) | ...
            cells.iscellstrempty(estimdata{5},[p.Results.release,'b_',date_particle]); 
        %ensure it exists
        assert(any(idx),['There is no GRACE ',p.Results.release,' solution for 20',date_particle])
        %ensure it is unique
        assert(sum(idx)==1,['There are multiple GRACE ',p.Results.release,' solutions for 20',date_particle])
        %pluck solution and split it into its directory path
        estimdir=strsplit(estimdata{5}{idx},'/');
        %append subdir to local estim_dir_root
        out=fullfile(p.Results.estim_dir_root,estimdir{8:10},'solution');
      end
    end
    function out=estim_filename(varargin)
      %parse optional arguments
      p=inputParser;p.KeepUnmatched=true;
      p.addParameter('estim_dir',     csr.estim_dir_default(varargin{:}), @(i) ischar(i));
      p.addParameter('arc',           [],  @(i) ((isnumeric(i) && ~isempty(i) && isscalar(i)) || strcmp(i,'*')) );
      p.addParameter('enforce_scalar',true,@(i)   islogical(i) && ~isempty(i) && isscalar(i));
      p.parse(varargin{:});
      out=file.wildcard(fullfile(p.Results.estim_dir,...
        ['R.common.arc',num2str(p.Results.arc),'.to.arc',num2str(p.Results.arc),'.*.estim']));
      if p.Results.enforce_scalar
        assert(numel(out)==1,['Can only handle one estim file at a time, not ',num2str(numel(out)),':',10,strjoin(out(:),char(10))])
        out=out{1};
      end
    end
    function out=estim_load(filename)
      if exist([filename,'.mat'],'file')
        load([filename,'.mat'])
      elseif strcmp(filename(end-3:end),'.mat')
        load(filename)
      else
        %init the outputs
        out=struct('A',[],'B',[]);
        %allocate record array
        dat=cell(file.length(filename),5);
        %open the file
        fid=file.open(filename);
        %load the data
        for i=1:size(dat,1)
          l=fgetl(fid);
%         00000000011111111112222222222333333333344444444445555555555666666666677777777778888888
%         12345678901234567890123456789012345678901234567890123456789012345678901234567890123456
%         GRC-A ZDOT      -1.883386301884850e+03  -9.058436043141965e-05  -1.883386392469210e+03
          dat{i,1    }=str.clean( l(1 :5  ),' ');
          dat{i,2    }=str.clean( l(6 :10 ),' ');
          dat(i,3:end)=cells.m2c(str.num(   l(11:end)    ));
        end
        %close the file
        fclose(fid);
        %parse the data
        for i=1:size(dat,1)
          fn=csr.estim_translate_param(dat{i,2});
          if ~isvarname(fn)
            fn=['f',strrep(fn,'-','_')];
          end
          out.(strrep(dat{i,1},'GRC-','')).(fn)=dat{i,5};
        end
        save([filename,'.mat'],'out')
      end
    end
    function T=arctimeline_load(filename)
      persistent arctimeline_table arctimeline_filename
      if isempty(arctimeline_table) || ~strcmp(arctimeline_filename,filename)
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
        arctimeline_table=table(year,month,day,arc,start,stop);
        %save filename of this arctimelines
        arctimeline_filename=filename;
      end
      T=arctimeline_table;
    end
    function out=arctimeline_get(varargin)
      %parse optional arguments
      p=inputParser;p.KeepUnmatched=true;
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
          str.say('Patched default arc for',...
            str.show(cellfun(@(i) {[i,':'],p.Results.(i),';'},{'year','month','arc','day','start','stop'},'UniformOutput',false))...
          )
        else
          out=[];
          str.say('No arc found for',...
            str.show(cellfun(@(i) {[i,':'],p.Results.(i),';'},{'year','month','arc','day','start','stop'},'UniformOutput',false))...
          )
        end
        return
      end
      relevant=arctimeline(key,par.out);
      %reduce from table the fields requested in par.out
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
    function [start,stop]=arctimeline_startstop(year,month)
      %call mother routine
      out=csr.arctimeline_get('year',year,'month',month,'out',{'start','stop'});
      %assign outputs
      start=out{1}(1);
      stop=out{2}(end);
    end
    function t=subarc_startstop(startstop,nr_sub_arcs)
      %zero sub-arcs need to be reduced to one
      if nr_sub_arcs==0;nr_sub_arcs=1;end
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
    function arc=arc_timeseries(startstop,calpars,varargin)
      %get persistent vars 
      persistent v units
      if isempty(v)
        % add input arguments and metadata to collection of parameters 'v'
        v=varargs.wrap('sources',{{...
          'calpar_prefix',         'AC0', @(i) ischar(i)    && ~isempty(i);...
          'calpar_suffixes',{'','D','Q'}, @(i) iscellstr(i) && ~isempty(i);...
          'calpar_coords'  {'X','Y','Z'}, @(i) iscellstr(i) && ~isempty(i);...
          'calpar_cyclesperday', [1,8,1], @(i) isnumeric(i) && ~isempty(i);...
          },csr.ltb_args},varargin{:});
        %these units are wrong but they need to agree with each other to build the cal mod (which is done elsewhere)
        units={'m/s^2','m/s^2','m/s^2'}; 
      end
      %get sats
      sats=fieldnames(calpars);
      %loop over all sats
      for s=1:numel(sats)
        %loop over all components
        for c=1:numel(v.calpar_coords)
          %loop over all suffixes
          for i=1:numel(v.calpar_suffixes)
            %build field name
            field=[v.calpar_prefix,upper(v.calpar_coords{c}),v.calpar_suffixes{i}];
            %build ordinate: get sub-arc start/stop
            t=csr.subarc_startstop(startstop,v.calpar_cyclesperday(c));
            %add a gap at the end of each arc (it's prettier)
            t=[t(:,1),t(:,2)-seconds(1),t(:,2)];
            %loop over all sub-arcs
            for day_loc=1:size(t,1)
              %retrieve (sub-) arc data
              if v.calpar_cyclesperday(c)>1
                field_now=[field,num2str(day_loc)];
              else
                field_now=field;
              end
              %skip if this calpar does not exists
              if ~isfield(calpars.(sats{s}),field_now)
                continue
              end
              y=ones(2,1)*calpars.(sats{s}).(field_now);
              %build timeseries object (with gap at the last entry)
              ts=simpletimeseries(transpose(t(day_loc,:)),[y;NaN],...
                'format','datetime',...
                'labels',{field_now},...
                'units',units(i),...
                'timesystem','gps',...
                'descriptor',field_now...
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
    function calmod=arc_build(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'year',                     [], @(i) isnumeric(i) && isscalar(i);...
        'month',                    [], @(i) isnumeric(i) && isscalar(i);...
        'arc',                      [], @(i) isnumeric(i) && isscalar(i);...
        'sats',              {'A','B'}, @(i) iscellstr(i) && ~isempty(i);...
        'desc',              'unknown', @(i) ischar(i)    && ~isempty(i);...
        'timestep',         seconds(5), @(i) isduration(i)&& ~isempty(i);...
        'calpar_coords'  {'X','Y','Z'}, @(i) iscellstr(i) && ~isempty(i);...
        'calpar_prefix',         'AC0', @(i) ischar(i)    && ~isempty(i);...
        'calpar_suffixes',{'','D','Q'}, @(i) iscellstr(i) && ~isempty(i);...
        'calpar_meaning', {'bias','linear','quad'}, @(i) iscellstr(i) && ~isempty(i);...
      },csr.ltb_args},varargin{:});
      %calpar_meaning is used to connect with bias, linear, quadratic time functions
      %therefore, it have the same length as calpar_suffixes
      assert(numel(v.calpar_meaning)==numel(v.calpar_suffixes),['The number of elements in ''calpar_suffixes'' must be ',...
        'the same as in ''calpar_meaning''.'])
      %derived parameters
      units=cell(size(v.calpar_coords));units(:)={'m/s^2'};
      %get start/stop of current arc(s)
      startstop=cells.c2m(...
          csr.arctimeline_get(varargin{:},...
            'year', v.year,...
            'month',v.month,...
            'arc',  v.arc,...
            'out',{'start','stop'}...
        )...
      );
      %check if this arc is part of the final solution (if not, then csr.arctimeline_get returned empty startstop)
      if isempty(startstop)
        calmod=[];
        return
      end
      %load calpars
      calpars=csr.estim_load(csr.estim_filename(varargin{:}));
      %build calibration parameters time series
      arc=csr.arc_timeseries(startstop,calpars);
      %define field names (this is arbitrary but must have the same length as calpar_suffixes)
      fields=v.calpar_meaning ;
      %set calibration model time domain
      t=startstop(1):seconds(1):startstop(2);
      %set the time domain in units compatible with the calibration parameters:
      %units are days and zero epoch is the start of the arc
      tc=days(t-startstop(1));
      %loop over all satellites
      for s=1:numel(v.sats)
        %gather calibration parameters
        for c=1:numel(v.calpar_coords)
          %init this calpar
          cal.(v.calpar_coords{c})=struct();
          %loop over all calpars for this coordinate and sat
          for f=1:numel(fields)
            %build calpar field name (e.g. AC0X or AC0YD1)
            field_now=[v.calpar_prefix,v.calpar_coords{c},v.calpar_suffixes{f}];
            %make sure this calpar is available
            if isfield(arc.(v.sats{s}),field_now)
              %interpolate blindly over the whole time domain
              cal.(v.calpar_coords{c}).(fields{f})=arc.(v.sats{s}).(field_now).interp(t);
              %remove transient calpars, this is indicative of a gap (interpolation is done blindly over any gap length)
              cal.(v.calpar_coords{c}).(fields{f})=cal.(v.calpar_coords{c}).(fields{f}).mask_and([...
                     cal.(v.calpar_coords{c}).(fields{f}).mask(1);...
                diff(cal.(v.calpar_coords{c}).(fields{f}).y(:,1))==0]...
              );
            end
          end
        end
        %init calibration model time series
        calmod.(v.sats{s})=simpletimeseries(t,zeros(numel(t),numel(v.calpar_coords)),...
          'format','datetime',...
          'labels',v.calpar_coords,...
          'units',units,...
          'timesystem','gps',...
          'descriptor',['calibration model ',v.desc,', GRACE-',v.sats{s}]...
        );
        %build calibration model
        for c=1:numel(v.calpar_coords)
          for f=1:numel(fields)
            if isfield(cal.(v.calpar_coords{c}),fields{f})
              switch fields{f}
              case 'bias';   increment=cal.(v.calpar_coords{c}).(fields{f}).cols(1);
              case 'linear'; increment=cal.(v.calpar_coords{c}).(fields{f}).cols(1).times(tc   );
              case 'quad';   increment=cal.(v.calpar_coords{c}).(fields{f}).cols(1).times(tc.^2);
              otherwise
                error(['Calibration parameter with meaning ''',fields{f},''' is not implemented.'])
              end
              calmod.(v.sats{s})=calmod.(v.sats{s}).set_cols(c,calmod.(v.sats{s}).get_cols(c)+increment);
            end
          end
        end
        %apply long-term biases
        calmod.(v.sats{s})=csr.ltb_apply(v.varargin{:},...
          'ts',calmod.(v.sats{s}),...
          'sat',v.sats{s},...
          'field','all'...
        );
      end
    end
    function out=arcs_build(varargin)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'year',     [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i);...
        'month',    [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i);....
        'day',      [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i);....
        'arc',      [], @(i) isnumeric(i) && ~isempty(i) && isscalar(i);...
      }},varargin{:});
      %translate days into arcs (estim files are arc-wise)
      if isempty(v.arc)
        if ~isempty(v.day)
          %get arcs of current day
          v.arc=cells.c2m(csr.arctimeline_get(...
            'arc_timeline_file',v.arc_timeline_file,...
            'year', v.year,...
            'month',v.month,...
            'day',  v.day,...
            'out',{'arc'}...
          ));
        %if no day or arc is given, then need to find them all
        else
          %get all arcs in estim_dir
          estim_files=csr.estim_filename(varargin{:},...
            'arc','*',...
            'enforce_scalar',false...
          );
          %parse filenames to retrieve arc numbers
          arcs=cell(size(estim_files));
          for i=1:numel(arcs)
            [~,f]=fileparts(estim_files{i});
            arcs{i}=strsplit(strrep(f,'R.common.arc',''),'.to.arc');
            arcs{i}=str2double(arcs{i}{1});
          end
          v.arc=sort(cells.c2m(arcs));
        end
      end
      %handle scalar inputs
      n=max([...
        numel(v.year),...
        numel(v.month),...
        numel(v.day),...
        numel(v.arc)...
      ]);
      for i={'year','month','day','arc'}
        v.(i{1})=cells.deal(v.(i{1}),[n,1]);
      end
      %sanity: all inputs must have the same length
      assert(all(cellfun(@(i) numel(v.(i)),{'year','month','day','arc'})==n),...
        'All time parameters must have the same length, debug needed!')
      %init outputs
      out=[];
      %loop over all inputs
      s.msg='Assembling arcs';s.n=n;
      for i=1:n
        tmp=csr.arc_build(varargin{:},...
          'year',     v.year{i},...
          'month',    v.month{i},...
          'day',      v.day{i},...
          'arc',      v.arc{i}...
        );
        %save if this arc is not empty
        if ~isempty(tmp)
          if isempty(out)
            out=tmp;
          else
            out=structs.objmethod2('append',out,tmp);
          end
        end
        s=time.progress(s);
      end
    end
    function obj=import_estim(obj,product,varargin)
      %This contructor needs the following arguments:
      % - 'estim_dir' (if empty, it retrieves the RL05 estim data)
      % - 'year'
      % - 'month'
      %
      %A complete collection of estim files relevant to 'year' and 'month' 
      %should be sitting obidiently inside 'estim_dir'.
      %
      %NOTICE: it is important that these values are in agreement between
      %each other, there is no obvious way of making this internally consistent.
      %
      %The mandatory argument 'level' serves as a sort of descriptor, to 
      %differentiate different types of estim files.
      
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'desc',  product.str, @(i) ischar(i);...
        'year',  year( obj.start+(obj.stop-obj.start)/2), @(i) isnumeric(i) && isscalar(i);...
        'month', month(obj.start+(obj.stop-obj.start)/2), @(i) isnumeric(i) && isscalar(i);...
      },product.metadata},varargin{:});
      %build the arc time series 
      arcs=csr.arcs_build(v.varargin{:});
      %save it
      obj=obj.data_set(product,arcs);
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    %% operators
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
      
      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        't'   ,1,        @(i) isnumeric(i) && ~isempty(i);...
        't2'  ,1,        @(i) isnumeric(i) && ~isempty(i);...
        'sats',{'A','B'},@(i) iscellstr(i) && ~isempty(i);...
      },product.metadata},varargin{:});

      %loop over the sats
      for s=1:numel(v.sats)
        %gather quantities, l1b acc data only contains the satellite names (located at the leafs of cal par data)
        acc=obj.data_get_scalar(l1baccp.dataname.set_field_path(v.sats(s)));
        %handle exceptions (also deals with non-existing data)
        if ~isa(acc,'simpletimeseries')
          %patch nan calibration model
          calmod=simpletimeseries(...
            [obj.start;obj.stop],...
            nan(2,numel(coords)),...
            'descriptor',['calibration model ',product.str,', GRACE-',v.sats{s},' (empty)']...
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
fields{1},obj.data_get_scalar(calparp.dataname.set_field_path([product.dataname.field_path,{['AC0',coords{c}    ]},v.sats(s)])),...
fields{2},obj.data_get_scalar(calparp.dataname.set_field_path([product.dataname.field_path,{['AC0',coords{c},'D']},v.sats(s)])),...
fields{3},obj.data_get_scalar(calparp.dataname.set_field_path([product.dataname.field_path,{['AC0',coords{c},'Q']},v.sats(s)])) ...
          );
        end
        %init models container
        calmod=simpletimeseries(acc.t,zeros(acc.length,numel(coords))).copy_metadata(acc);
        calmod.descriptor=['calibration model ',product.str,', GRACE-',v.sats{s}];
        str.say('Computing the ',calmod.descriptor)
        for c=1:numel(coords)
          obj.log('@','iter','product',product,'sat',v.sats{s},'coord',coords{c})
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
            cal.(coords{c}).(fields{1}).cols(param_col                              )+...
            cal.(coords{c}).(fields{2}).cols(param_col).times(t.(fields{2})   * v.t )+...
            cal.(coords{c}).(fields{3}).cols(param_col).times(t.(fields{3}).^2* v.t2)...
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
              'suffix',[v.sats{s},'.',coords{c}]...
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
        obj=obj.data_set(product.dataname.append_field_leaf(v.sats{s}),calmod);
      end
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function obj=estimate_poly_calmod(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)

      %sanity
      assert(product.nr_sources==2,...
        ['number of sources in product ',product.str,...
        ' is expected to be 2, not ',num2str(product.nr_sources),'.'])
      %get sources metadata
      modaccp=obj.product_get(product.sources(1));
      l1baccp=obj.product_get(product.sources(2));

      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{product.metadata,{...
        'coords',    {'X','Y','Z'}, @(i) iscellstr(i) && ~isempty(i);...
        'sats',          {'A','B'}, @(i) iscellstr(i) && ~isempty(i);...
        'calpar_cyclesperday',    [1,8,1], @(i) isnumeric(i) && ~isempty(i);...
        'calpar_poly_order', [2,2,2], @(i) isnumeric(i) && ~isempty(i);...
        'l1b_median_order',     10, @(i) isnumeric(i) && ~isempty(i) && isscalar(i);...
      }},varargin{:});
      
      %check if the product satellites include the one in the current product
      sat_now=product.dataname.field_path;
      if ~any(strcmp(v.sats,sat_now))
        obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop,'sats',v.sats,'sat from product',sat_now)
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
      switch v.arc_scheme
      case 'blind'
        startstop=[t0,t1];
      case 'csr'
        startstop=[];
        for i=1:numel(t0)
          startstop=[startstop;...
            cells.c2m(...
              csr.arctimeline_get(...
                'arc_timeline_file',v.atl_file,...
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
        error(['Cannot handle ''arc_scheme'' with value ''',v.arc_scheme,'''.'])
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
          nan(2,numel(v.coords)),...
          'descriptor',['calibration model ',product.str,' (empty)']...
        );
        str.say('Skipping  the ',calmod.descriptor)
      else
        str.say('Computing the calibration model ',product)
        %loop over all coordinates
        for c=1:numel(v.coords)
          col_idx=cells.cellstrfind({'X','Y','Z'},v.coords{c});
          %loop over all relevant (sub-) arcs
          for d=1:size(startstop,1)
            %retrieve (sub-) arcs
            t=csr.subarc_startstop(startstop(d,:),v.calpar_cyclesperday(c));
            %loop over the number of arcs per day
            for a=1:size(t,1)
              obj.log('@','iter',...
                'product',product,'coord',v.coords{c},...
                'arc',startstop(d,:),...
                'sub-arc',t(a,:),'sub-arc',a)
              %get trimmed time series
              arc_mod=mod.trim(t(a,1),t(a,2)).cols(c);
              arc_l1b=l1b.trim(t(a,1),t(a,2)).cols(c).median(v.l1b_median_order).interp(arc_mod.t);
              
              %get calibration curve to fit
              arc_cur=arc_mod+arc_l1b;
              %compute calibration model
              warning('off',warning_id)
              arc_calmod=arc_cur.scale(-1).polyfit(v.calpar_poly_order(c));
              warning('on',warning_id)
              
%               %get calibration curve to fit
%               arc_cur1=arc_mod-arc_l1b;
%               arc_cur2=arc_mod+arc_l1b;
%               %compute calibration model
%               warning('off',warning_msg)
%               [arc_calmod1,S1]=arc_cur1.polyfit(v.calpar_poly_order(c));
%               [arc_calmod2,S2]=arc_cur2.polyfit(v.calpar_poly_order(c));
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
    function obj=smooth(obj,product,varargin)
      %branch on implicit mode
      if isempty(product)
        % direct mode: the data is given in obj
        dat=obj;
        % use this function to show messages
        dispfun=@(i) str.say('stack_delta',1,i{:});
        % define dummy variable v
        v={};
      else
        obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
        %sanity
        assert(product.nr_sources==1,...
          ['number of sources in product ',product.str,...
          ' is expected to be 1, not ',num2str(product.nr_sources),'.'])
        % get the sats
        v=varargs.wrap('sources',{product.metadata,{...
          'sats',          {'A','B'}, @(i) iscellstr(i) && ~isempty(i);...
        }},varargin{:});
        %check if the product satellites include the one in the current product
        sat_now=product.dataname.field_path;
        if ~any(strcmp(v.sats,sat_now))
          obj.log('@','out (skipped)','product',product,'start',obj.start,'stop', obj.stop,'sats',v.sats,'sat from product',sat_now)
          %nothing to do
          return
        end
        % get the data
        dat=obj.data_get_scalar(product.sources(1));
        % use this function to show messages
        dispfun=@(i)  obj.log('@',i{:});
        % parse the metadata
        v=varargs.wrap('sources',{product.metadata},varargin{:});
      end
      % parse the resample period
      v=varargs.wrap('sources',{v,{...
        'smooth',        0, @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'median',        0, @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'resample',      0, @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'bandpass',[0,inf], @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'trim_edges', true, @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
      }},varargin{:});
      % define the operations (they are applied in this order!)
      ops={'resample','median','smooth','bandpass'};
      % operate the data
      for i=1:numel(ops)
        % convert to duration
        if ~isduration(v.(ops{i}));v.(ops{i})=seconds(v.(ops{i}));end
        % only operate if the time step is smaller than the requested period
        if any(v.(ops{i})(isfinite(v.(ops{i}))) > 0)
          dispfun({ops{i},'t_span',v.(ops{i})})
          dat=dat.(ops{i})(v.(ops{i}));
        else
          dispfun({ops{i},'t_span',v.(ops{i}),'skipped because data step is',dat.step})
        end
      end
      if v.trim_edges
        %find trim length
        trim_len=max([v.smooth,v.median]);
        if ~isduration(trim_len);trim_len=seconds(trim_len);end
        if ~all(v.bandpass~=[0,inf])
          trim_len=max([trim_len,(dat.stop-dat.start)*0.05]);
        end
        % trim edges
        dat=dat.trim(dat.start+trim_len,dat.stop-trim_len);
      end
      %branch on implicit mode
      if isempty(product)
        % direct mode: back propagate
        obj=dat;
      else
        %propagate it
        obj=obj.data_set(product,dat);
        obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
      end
    end
    %% temperature correction
    function obj=estimate_temp_corr(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop', obj.stop)
      %some internal parameters
      smooth_default=hours(6);
      col=3;
      col_ref=1;
      xn={'k_R','k_C','dt'};

      %sanity
      assert(product.nr_sources==4,...
        ['number of sources in product ',product.str,...
        ' is expected to be 4, not ',num2str(product.nr_sources),'.'])

      % add input arguments and metadata to collection of parameters 'v'
      v=varargs.wrap('sources',{{...
        'coords',        {'X','Y','Z'}, @(i) iscellstr(i) && ~isempty(i);...
        'sats',              {'A','B'}, @(i) iscellstr(i) && ~isempty(i);...
        'bandpass',            [inf,0], @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'smooth',             hours(3), @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'median',                    0, @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'resample',         minutes(1), @(i)(isnumeric(i) || isduration(i)) && ~isempty(i);...
        'searchrange',     [100,100,2], @(i) isnumeric(i)  && isscalar(i) && ~isempty(i);...
        'searchlen',               100, @(i) isnumeric(i)  && isscalar(i) && ~isempty(i);...
        'RelTol',                 1e-4, @(i) isnumeric(i)  && isscalar(i) && ~isempty(i);...
        'AbsTol',                 1e-4, @(i) isnumeric(i)  && isscalar(i) && ~isempty(i);...
        'x0',                 [1,1,20], @(i) isnumeric(i)  && numel(i)==3;...
        'xs',           [1e-14,1e-6,1], @(i) isnumeric(i)  && numel(i)==3;...
        'ahk_col',                   1, @(i) isnumeric(i)  && isscalar(i);...
        'test',                  false, @(i) islogical(i)  && isscalar(i) && ~isempty(i);...
        'test_plot',             true, @(i) islogical(i)  && isscalar(i) && ~isempty(i);...
      },product.metadata},varargin{:});
    
      %check if the product satellites include the one in the current product
      sat_now=product.dataname.field_path;
      if ~any(strcmp(v.sats,sat_now))
        obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop,'sats',v.sats,'sat from product',sat_now)
        %nothing to do
        return
      end
      % loop over all data
      ts=cell(1,product.nr_sources);
      for i=1:product.nr_sources
        % load and smooth the data
        ts{i}=csr.smooth(obj.data_get_scalar(product.sources(i)),[],v.varargin_for_wrap{:});
      end
      %set clearer labels
      ts{1}.labels={'L1B','L1B','L1B'};
      ts{2}.labels={'mod','mod','mod'};
      ts{3}.labels={'calmod','calmod','calmod'};
      %enforce common gaps
      ts=simpledata.op_multiple({'merge','mask_match'},ts,'estimate_temp_corr');
      %make room for outputs
      tsout=simpletimeseries.unitc(ts{1}.t,ts{1}.width,...
        'units',ts{1}.y_units,'descriptor',product.name,'timesystem',ts{1}.timesystem);
      % get forcing temperature
      tsf=ts{4}.cols(v.ahk_col);
      % init internal records (2 comes from the two target strategies)
      rnorm=zeros(1,2);sT=rnorm;bT=rnorm;
      x=cell(size(rnorm));tsc_try=x;tst_try=x;y_pop=x;x_pop=x;
      % loop over all coordinates
      for c=1:numel(v.coords)
        %get index of this coordinate
        ci=cells.cellstrfind({'x','y','z'},lower(v.coords{c}));
        for i=1:2
          switch i
          case 1
            %create target ts 1
            tst_try{i}=ts{1}.cols(ci)-ts{2}.cols(ci);            
          case 2
            %create target ts 2
            tst_try{i}=ts{3}.cols(ci);
            switch lower(v.coords{c})
              case 'x'; %do nothing
              case 'y'; tst_try{i}=tst_try{i}.scale(-1);
              case 'z'; %do nothing
            end
          end
          %"solve"
          [x{i},rnorm(i),tsc_try{i},sT(i),bT(i),y_pop{i},x_pop{i}]=csr.temp_corr_solve(tst_try{i},tsf,v);
          %apply empirical bias and scales
          tsc_try{i}=(tsc_try{i}+bT(i)).*sT(i);
        end
        
        %pick the target ts with the smallest residuals (i.e. the best try index)
        bti=find(rnorm==min(rnorm));
        %make nice confusing plot, if requested
        if v.test_plot
          %build plot filename
          plot_name=[strjoin({...
            product.codename,...
            [datestr(obj.start,'yyyymmdd'),'T',datestr(obj.stop,'yyyymmdd')],...
            v.coords{c},...
            ['Ta',num2str(v.ahk_col)],...
            structs.str(v,{{'bandpass'},{'smooth'},{'median'},{'resample'}},'_'),...
            structs.str(v,{{'RelTol'},{'AbsTol'}},'_'),...
            structs.str(v,{{'searchrange'},{'searchlen'}},'_'),...
          },'-'),'.png'];
          %check if it exists
          if ~exist(plot_name,'file')
            %build info text 
            text_str=cell(numel(xn)+3,1);
            for ti=1:numel(xn)
              text_str{ti}=[xn{ti},': ',str.show(x{bti}(:,ti),'',', '),'\times10^{',num2str(log10(v.xs(ti))),'}'];
            end
            text_str{ti+1}=['rnorm: ',str.show(rnorm(bti),'%.1f',', ')];
            text_str{ti+2}=['n: ',str.show(numel(y_pop{bti}),'%d',', ')];
            text_str{ti+3}=tst_try{bti}.labels{1};
            %plot ot
            plotting.figure(v.varargin_for_wrap{:});
            h{1}=tsf.plot(           'zeromean',true,'normalize',true);
            h{2}=tst_try{1}.plot(    'zeromean',true,'normalize',true);
            h{3}=tsc_try{1}.plot(    'zeromean',true,'normalize',true);
            h{4}=tst_try{2}.plot(    'zeromean',true,'normalize',true);
            h{5}=tsc_try{2}.plot(    'zeromean',true,'normalize',true);
            plotting.enforce(...
              'plot_legend',cellfun(@(i)i.legend{1},h,'UniformOutput', false),...
              'plot_legend_location','northeast',...
              'plot_fontsize_legend',14,...
              'plot_xdate',true,...
              'plot_xdateformat','dd',...
              'plot_xlabel',datestr(obj.start+(obj.stop-obj.start)/2,'mmm yyyy'),...
              'plot_title',[v.coords{c},'-axis GRACE-',product.dataname.field_path{1}]...
            );
            %add text box
            a=axis;
            text(a(1)+0.01*(a(2)-a(1)),a(3)+0.15*(a(4)-a(3)),strjoin(text_str,char(10)),'fontsize',14)
            %save plots
            saveas(gcf,plot_name)
          end
          %plot histograms of parameters and rnorm
          if ~isempty(y_pop)
            %redefine plot name
            plot_name=strrep(plot_name,'.png','.hist.png');
            %check if it exists
            if ~exist(plot_name,'file')
              %plotting histograms
              plotting.dhist(log10(x_pop{bti}(:,1)),log10(x_pop{bti}(:,2)),log10(y_pop{bti}),...
                v.varargin_for_wrap{:},...
                'xlabel',['log(',xn{1},')'],...
                'ylabel',['log(',xn{2},')'],...
                'zlabel','log(rnorm)')
              %save plots
              saveas(gcf,plot_name)
            end
          end
        end
        %save this columns
        tsout=tsout.set_cols(ci,tsc_try{bti});
      end
      %propagate it
      obj=obj.data_set(product,tsout);
      obj.log('@','out','product',product,'start',obj.start,'stop', obj.stop)
    end
    function [x,rnorm,tsc,sT,bT,y_pop,x_pop]=temp_corr_solve(ts_target,ts_forcing,v)
      %compute temperature-acceleration scale
      sT=ts_target.stats('mode','ampl')/ts_forcing.stats('mode','ampl');
      %scale target ts
      ts_target=ts_target.scale(1/sT);
      %compute temperature-acceleration bias
      bT=ts_target.stats('mode','mean')-ts_forcing.stats('mode','mean');
      %scale target ts
      ts_target=ts_target-bT;
      %assign easier names
      Tt=ts_target.y_masked;
      tt=ts_target.x_masked;
      Ta=ts_forcing.y_masked;
      ta=ts_forcing.x_masked;
      assert(all(abs(tt-ta)<1e-6),'time domain discrepancy')
      if v.test
        x=v.x0;
        rnorm=csr.temp_corr_fitness(ta,Ta,Tt,x,v);
        sT=zeros(1,2);bT=sT;y_pop=cell(1,2);x_pop=y_pop;
      else
        %call "solver"
        [x,y_pop,x_pop]=num.param_brute(...
          @(i) csr.temp_corr_fitness(ta,Ta,Tt,i,v),...
          v.x0./v.searchrange,...
          v.x0.*v.searchrange,...
          v.varargin{:},...
          'searchspace',    'log',...
          'vectorize',       true ...
        );
        rnorm=min(y_pop);
      end
      %get corrected accelerations
      [T,t]=csr.temp_corr_Tb(ta,Ta,x,v);
      %build ts object
      tsc=simpletimeseries(ts_forcing.x2t(t),T,'labels',{['temp eff in ',ts_target.labels{1}]});
    end
    function out=temp_corr_dTb(t,T,ta,Ta,x,xs)
      Ta_now=interp1(ta,Ta,t);
      out=x(:,1)*xs(1).*( (Ta_now).^4 - (T+x(:,3)*xs(3)).^4 ) + x(:,2)*xs(2).*( (Ta_now) - (T+x(:,3)*xs(3)) );
    end
    function [T,t]=temp_corr_Tb(ta,Ta,x,v)
      options = odeset('Vectorized','on','RelTol',v.RelTol,'AbsTol',v.AbsTol);
      Ta=Ta+273.15;
      [t,T] = ode45(@(ti,Ti) csr.temp_corr_dTb(ti,Ti,ta,Ta,x,v.xs) , ta, Ta(1)-x(:,3)*v.xs(3) , options );
      T=T-273.15;
    end
    function out=temp_corr_fitness(ta,Ta,Tt,x,v)
      Tb=csr.temp_corr_Tb(ta,Ta,x,v);
      res=Tb-Tt(:)*ones(1,size(x,1));
      res=res-ones(size(res,1),1)*mean(Tb);
      res=res+mean(Tt);
      out=sum(res.^2);
    end
    %% manual interfaces
    function obj=plot(start,stop,product,plot_dir)
      switch start
      case 'odd'
        date_list={...
          '2002-04-27';...
          '2002-05-16';...
          '2002-08-05';...
          '2002-08-06';...
          '2002-10-07';...
          '2002-11-15';...
          '2002-11-27';...
          '2002-12-04';...
          '2002-12-05';...
          '2002-12-13';...
          '2002-12-18';...
          '2003-01-01';...
          '2003-01-14';...
          '2003-02-11';...
          '2003-02-22';...
          '2004-04-05';...
          '2004-12-01';...
          '2005-02-24';...
          '2005-11-15';...
          '2006-01-26';...
          '2006-06-15';...
          '2011-03-28';...
          '2011-04-07';...
          '2011-08-25';...
          '2013-01-21';...
          '2014-01-08';...
          '2014-04-14';...
          '2014-04-18';...
          '2014-04-22';...
          '2014-04-25';...
          '2014-08-16';...
          '2014-09-19';...
          '2014-09-25';...
          '2015-02-10';...
          '2015-08-14';...
          '2015-12-11';...
          '2016-05-18';...
          '2016-05-25';...
          '2016-11-30';...
          '2017-01-13';...
        };
        for i=1:numel(date_list)
          tmp=csr.plot(date_list{i},'','odd');
          clear tmp
        end
      otherwise 
        if ~exist('stop','var') || isempty(stop)
          stop=dateshift(datetime(start),'end','day')-seconds(1);
        end
        product_default=dataproduct('grace.acc.cal.csr.plots/estim');
        if ~exist('product','var') || isempty(product)
          product=product_default;
        elseif strcmp(product,'odd')
          product=product_default;
          plot_dir=['/Volumes/Users/Teixeira/2017-11-17.Oddities/',datestr(start,'yyyy-mm-dd')];
        end
        if ~exist('plot_dir','var') || isempty(plot_dir)
          plot_dir=product.plot_dir;
        end
        product.plot_dir=plot_dir;
        obj=datastorage('start',datetime(start),'stop',datetime(stop),'debug',true).init(product);
        if ~product.mdget('plot_visible')
          close all
        end
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