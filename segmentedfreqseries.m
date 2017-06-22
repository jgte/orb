classdef segmentedfreqseries < simplefreqseries
  %static
  properties(Constant,GetAccess=private)
    %default value of input parameters
    parameter_list=struct(...
      'seg_length', struct('default',days(1), 'validation',@(i) isduration(i) && isscalar(i)),...
      'seg_overlap',struct('default',hours(3),'validation',@(i) isduration(i) && isscalar(i))...
    );
  end
  %read only
  properties(SetAccess=private)
    seg_length
    seg_overlap
    seg
  end
  %private (visible only to this object)
  properties(GetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
    y_merged
  end
  methods(Static)
    function out=parameters
      out=fieldnames(segmentedfreqseries.parameter_list);
    end
    function [out,idx]=time_segmented(time,seg_length,seg_overlap)
      if ~isdatetime(time)
        error([mfilename,'input ''time'' must be of class ''datetime'', not ''',class(time),'''.'])
      end
      if ~isduration(seg_length)
        error([mfilename,'input ''seg_length'' must be of class ''duration'', not ''',class(seg_length),'''.'])
      end
      if ~isduration(seg_overlap)
        error([mfilename,'input ''seg_overlap'' must be of class ''duration'', not ''',class(seg_overlap),'''.'])
      end
      if seg_overlap>=seg_length
        error([mfilename,'input ''seg_overlap'' (',num2str(seg_overlap),') must be smaller than input ''seg_length'' (',num2str(seg_length),').'])
      end
      %handle infinite segment length
      if ~isfinite(seg_length)
        seg_length=time(end)-time(1);
        seg_overlap=seconds(0);
      end
      %guess the number of segments
      n=ceil((time(end)-time(1))/seg_length*2);
      %init outputs
      out=cell(1,n);
      idx=cell(1,n);
      %init loop vars
      start=time(1);
      stop=time(1);
      c=0;
      %loop it
      while stop < time(end)
        %increment counter
        c=c+1;
        %determine stop time
        stop=start+seg_length;
        %find indexes
        idx{c}=[...
               find(time>=start,1,'first'),...
          min([find(time<=stop, 1,'last' ),numel(time)])...
        ];
        %save segment time domain
        out{c}=time(idx{c}(1):idx{c}(2));
        %set start time
        start=stop-seg_overlap;
      end
      %remove empty entries
      empty_idx=cells.isempty(out);
      out=out(~empty_idx);
      idx=idx(~empty_idx);
    end
    %general test for the current object
    function out=test_parameters(field,l,w)
      switch field
      case 'timestep'
        out=minutes(10);
      case 'time'
        timestep=segmentedfreqseries.test_parameters('timestep',l,w);
        out=datetime('2015-01-01')+(-l/2*timestep:timestep:(l/2-1)*timestep);
      case 'y-dyn-sin'
        t=juliandate(segmentedfreqseries.test_parameters('time',l,w),'modifiedjuliandate');
        t=t(:);
        % number of frequencies in the signal
        nf=1;
        %build frequency content of signal, starting and ending frequencies
        f=[0.1;0.3]*(l./logspace(0,1,nf));
        %build time series of instantaneous frequencies
        ft=interp1([t(1);t(end)],f,t);
        %make room for output
        out=zeros(l,w);
        %build signal
        for i=1:nf
          out=out+ones(l,1)*randn(1,w)*10; %random bias
          out=out+sin(2*pi*(t.*ft(:,i)))*(1+10*rand(1,w)); %sinusoidal signal (random amplitudes)
          out=out+randn(l,w)*0.1; %noise
        end
      otherwise
        out=simplefreqseries.test_parameters(field,l,w);
      end
    end
    function test(l,w)

      if ~exist('l','var') || isempty(l)
        l=1000;
      end
      if ~exist('w','var') || isempty(w)
        w=3;
      end

      %test current object
      args=segmentedfreqseries.test_parameters('args',l,w);
      columns={'columns',1};
      line_seg={'line',{'-o'}};
      line_y  ={'line',{'-+'}};

      %simple plot

      a=segmentedfreqseries(...
        segmentedfreqseries.test_parameters('time',l,w),...
        segmentedfreqseries.test_parameters('y-dyn-sin',l,w),...
        'format','datetime',...
        args{:}...
      );

%       figure
%       subplot(2,1,1)
%       a.plot(columns{:},line_y{:})
%       a.op('plot',columns{:},line_seg{:})
%
%       subplot(2,1,2)
%       a.plot_psd(columns{:},line_y{:})
%       a.op('plot_psd',columns{:},line_seg{:})

      %y_merged

      figure
      b=simplefreqseries(...
        a.t,...
        a.y_merged,...
        args{:},...
        'format','datetime'...
      );
      b.plot(columns{:},line_seg{:});
      a.plot(columns{:},line_y{:})
      c=a-b;
      if sum(c.y(:).^2)>1e-14
        error([mfilename,':BUG TRAP: discrepancy between merged and original data'])
      end


    end 
  end
  methods
    %% constructor
    function obj=segmentedfreqseries(t,y,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 't',     @(i)                 ~isscalar(i)); %this can be char, double or datetime
      p.addRequired( 'y',     @(i) isnumeric(i) && ~isscalar(i));
      %declare parameters
      for j=1:numel(segmentedfreqseries.parameters)
        %shorter names
        pn=segmentedfreqseries.parameters{j};
        %declare parameters
        p.addParameter(pn,segmentedfreqseries.parameter_list.(pn).default,segmentedfreqseries.parameter_list.(pn).validation)
      end
      % parse it
      p.parse(t,y,varargin{:});
      %call superclass
      obj=obj@simplefreqseries(p.Results.t,p.Results.y,varargin{:});
      % save parameters
      for i=1:numel(segmentedfreqseries.parameters)
        %shorter names
        pn=segmentedfreqseries.parameters{i};
        if ~isscalar(p.Results.(pn))
          %vectors are always lines (easier to handle strings)
          obj.(pn)=transpose(p.Results.(pn)(:));
        else
          obj.(pn)=p.Results.(pn);
        end
      end
      % remove mask from varargin
      args=simpledata.vararginclean(varargin,{'mask'});
      % cut into segments
      obj=obj.segmentate(p.Results.seg_length,p.Results.seg_overlap,args{:});
    end
    function obj=copy_metadata(obj,obj_in)
      %call superclass
      obj=copy_metadata@simplefreqseries(obj,obj_in);
      %propagate parameters of this object
      parameters=fields(segmentedfreqseries.parameter_list);
      for i=1:numel(parameters)
        if isprop(obj,parameters{i}) && isprop(obj_in,parameters{i})
          obj.(parameters{i})=obj_in.(parameters{i});
        end
      end
    end
    %% seg methods
    function obj=segmentate(obj,seg_length,seg_overlap,varargin)
      % separate time series into segments
      [ts,idx]=segmentedfreqseries.time_segmented(obj.t,...
        seg_length,...
        seg_overlap...
      );
      % initialize
      obj.seg=cell(size(ts));
      s.msg=[mfilename,': cutting into segments'];s.n=numel(ts);
      % propagate segments
      for i=1:numel(ts)
        %initialize
        obj.seg{i}=simplefreqseries(...
          ts{i},...
          obj.y(idx{i}(1):idx{i}(2),:),...
          'mask',obj.mask(idx{i}(1):idx{i}(2),:),...
          varargin{:}...
        );
        %copy metadata from unsegmented object
        obj.seg{i}=obj.seg{i}.copy_metadata(obj);
        %update descriptor
        obj.seg{i}.descriptor=['segment ',num2str(i),' of ',obj.descriptor];
        %inform
        s=time.progress(s,i);
      end
    end
    %return the segmend idx at time t
    function seg_idx=segment_idx(obj,t)
      seg_idx=[];
      for i=1:numel(obj.seg)
        if obj.seg{i}.t(1)<=t && t<=obj.seg{i}.t(end)
          seg_idx(end+1)=i; %#ok<AGROW>
        end
      end
    end
    %returns the segment(s) at time t
    function seg=segment(obj,t)
      seg_idx=segment_idx(obj,t);
      % outputs
      if isempty(seg_idx)
        seg=cell(0);
      else
        seg=obj.seg(seg_idx);
      end
    end
    %apply mask also in segments
    function obj=apply_mask(obj,start,stop)
      %call superclass
      obj=apply_mask@simpledata(obj,start,stop);
      %propagate to segments
      for i=1:numel(obj.seg)
        obj.seg{i}=obj.seg{i}.apply_mask(start,stop);
      end
    end
    %% y_merge
    function y=get.y_merged(obj)
      %make room for outputs
      y=nan(size(obj.y));
      %loop over all segments
      s.msg=[mfilename,': merging segments'];s.n=numel(obj.seg);
      for i=1:numel(obj.seg)

        %propagate overlapping values

        if i > 1
          %get starting and ending indexes of overlap relative to the
          %previous and current segments
          idx_curr_start=1;
          idx_prev_stop =numel(obj.seg{i-1}.t);
          idx_prev_start=obj.seg{i-1}.idx(obj.seg{i  }.t(idx_curr_start));
          idx_curr_stop =obj.seg{i  }.idx(obj.seg{i-1}.t(idx_prev_stop ));
          %bug trap
          if (idx_curr_stop-idx_curr_start) ~= (idx_prev_stop-idx_prev_start)
            error([mfilename,': BUG TRAP: discrepancy between the length of overlaping segments'])
          end
          %retrive previous and current overlaps
          y_prev=obj.seg{i-1}.y(idx_prev_start:idx_prev_stop,:);
          y_curr=obj.seg{i  }.y(idx_curr_start:idx_curr_stop,:);
          %overlap length
          n=size(y_prev,1);
          %build weights
          w_curr=transpose(linspace(0,1,n))*ones(1,obj.width);
          w_prev=flipud(w_curr);
          %bug trap
          if any(w_curr(:)+w_prev(:) ~= 1)
            error([mfilename,': BUG TRAP: weights do not add up to one'])
          end
          %get starting and ending indexes of overlap relative to the
          %global time domain
          idx_start=obj.idx(obj.seg{i  }.t(idx_curr_start));
          idx_stop =obj.idx(obj.seg{i-1}.t(idx_prev_stop ));
          %apply weights to each overlapping segment and add them together
          y(idx_start:idx_stop,:)=y_prev.*w_prev+y_curr.*w_curr;
        else
          %patch for the first segment
          idx_curr_stop=0;
        end

        %propagate non-overlapping section of current segment

        %get starting index of non-overlapping section relative to the
        %current segments
        idx_curr_start=idx_curr_stop+1;
        %get stopping index of non-overlapping section relative to the
        %current segments
        if i < numel(obj.seg)
          idx_curr_stop = obj.seg{i}.idx(obj.seg{i+1}.t(1))-1;
        else
          idx_curr_stop = numel(obj.seg{i}.t);
        end
        %get starting and ending indexes of non-overlapping section relative to the
        %global time domain
        idx_start=obj.idx(obj.seg{i}.t(idx_curr_start));
        idx_stop =obj.idx(obj.seg{i}.t(idx_curr_stop ));
        %copy data
        y(idx_start:idx_stop,:)=obj.seg{i}.y(idx_curr_start:idx_curr_stop,:);

        %user feedback
        s=time.progress(s,i);
      end
    end
    function obj=update_y_merged(obj)
      obj.y=obj.y_merged;
    end
    %% operate segment-wise
    function out=op(obj,operation,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'operation',                   @(i) ischar(i));
      p.addParameter('idx',        1:numel(obj.seg),@(i) isnumeric(i))
      p.addParameter('self_assign',false,           @(i) islogical(i) && isscalar(i))
      % parse it
      p.parse(operation,varargin{:});
      %init outputs
      out=cell(size(p.Results.idx));
      %loop over all segments
      for i=1:numel(p.Results.idx)
        out{i}=obj.seg{p.Results.idx(i)}.(p.Results.operation)(varargin{:});
      end
      %tweaking some super methods
      switch p.Results.operation
      case 'plot'
        xlim(datenum([obj.x(1) obj.x(end)]))
      end
      %check if outputs are of the same class as member 'seg'
      self_assign=true;
      for i=1:numel(out)
        self_assign=self_assign & strcmp(class(out{i}),class(obj.seg{p.Results.idx(i)}));
      end
      if self_assign && p.Results.self_assign
        %save self-assignement
        obj.seg=out;
        %update y
        obj.y=obj.y_merged;
        %update descriptor
        obj.descriptor=[obj.descriptor,' ',operation];
        %propagate
        out=obj;
      elseif ~self_assign && p.Results.self_assign
        error([mfilename,': cannot self-assign because there is a class discrepancy.'])
      end
      %clear outputs if not needed
      if nargout==0
        if p.Results.self_assign
          error([mfilename,': need one output argument if argument ''sefl_assign'' is true.'])
        end
        clear out
      end
    end
    %% overloading
    function obj=scale(obj,scale)
      obj=scale@simpledata(obj,scale);
      obj=obj.segmentate(obj.seg_length,obj.seg_overlap);
    end
    function obj=plus(obj,obj_new)
      obj=plus@simpledata(obj,obj_new);
      obj=obj.segmentate(obj.seg_length,obj.seg_overlap);
    end
    function obj=minus(obj,obj_new)
      obj=minus@simpledata(obj,obj_new);
      obj=obj.segmentate(obj.seg_length,obj.seg_overlap);
    end
    function obj=times(obj,obj_new)
      obj=times@simpledata(obj,obj_new);
      obj=obj.segmentate(obj.seg_length,obj.seg_overlap);
    end
    function obj=rdivide(obj,obj_new)
      obj=rdivide@simpledata(obj,obj_new);
      obj=obj.segmentate(obj.seg_length,obj.seg_overlap);
    end
    function obj=interp(obj,t_now,varargin)
      obj=interp@simplefreqseries(obj,t_now,varargin{:});
      obj=obj.segmentate(obj.seg_length,obj.seg_overlap,varargin{:});
    end
  end
end











