classdef attitude 
  %static
  properties(Constant,GetAccess=private)
    %list of data fields
    %NOTICE: needs updated when adding a new data type
    data_type_list={'quat','ang','angr','anga'};
    %default value of parameters
    %NOTICE: needs updated when adding a new parameter
    %NOTICE: it's not a bad idea to define attitude_dir in project.yaml, to have the proper value for data_dir.
    %        If there is no attitude_dir in project.yaml, the 'attitude' dir will need to exist in the same
    %        dir as this script.
    parameter_list={...
      'satname',   'unknown',@ischar;...
      'frame_from','srf',    @(i) simpletimeseries.isframe(i);...
      'frame_to',  'crf',    @(i) simpletimeseries.isframe(i);...
      'qsfirst',    true,    @islogical;...
      'verbose',   false,    @islogical;...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTICE: needs updated when adding a new data type (if relevant)
    compatible_parameter_list={'satname','frame_from','frame_to'};
  end
  %NOTICE: needs updated when adding a new parameter
  properties
    %parameters
    satname
    frame_from
    frame_to
    qsfirst
    verbose
  end
  %NOTICE: needs updated when adding a new data type
  properties(SetAccess=private)
    %data types
    quat
    angi
    angri
    angai
    initialized
  end
  %calculated only when asked for
  properties(Dependent)
    time
    ang
    angr
    anga
    vector_part
    scalar_part
  end
  methods(Static)
    %interface methods to object constants
    function out=data_types
      out=attitude.data_type_list;
    end
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(attitude.parameter_list); end
      out=v.picker(varargin{:});
    end
    %data source definitions
%     function out=nrtdm_product(in)
%       switch simpletimeseries.translatesat(in)
%         case 'ch'
%           out='CH_Basic/Orbit_CH-OG-3-RSO';
%         case 'ga'
%           out='GA_Basic/Orbit_NAVSOL';
%         case 'gb'
%           out='GB_Basic/Orbit_NAVSOL';
%         case 'sa'
%           out='SA_Basic/Orbit_L1B';
%         case 'sb'
%           out='SB_Basic/Orbit_L1B';
%         case 'sc'
%           out='SC_Basic/Orbit_L1B';
%         otherwise
%           error([mfilenane,': unknown NRTDM product for satellite ''',in,''', debug needed!'])
%       end
%     end
    %the function <format>_filename below define the filenames given sat, date and dir
    %for the data of different formats
    %NOTICE: data_dir is the top-most data dir, without specifying the satellite, data, etc
    function filename=grace_l1b_filename(satname,start,version,data_dir)
      filename=grace.grace_l1b_filename('SCA1B',satname,start,version,data_dir);
    end
    %NOTICE: add more <format>_filename wrapers here

    %wrapper for the <format>_filename routines, transparent for different <format>s
    %NOTICE: inputs satname and start can be arrays (cell array for satname)
    %NOTICE: the output is always a cell array
    function out=filename(format,satname,start,varargin)
      p=inputParser;
      p.addRequired( 'format',      @ischar);
      p.addRequired( 'satname',     @(i) ischar(i) || iscell(i));
      p.addRequired( 'start',       @(i) isdatetime(i) );
      p.addParameter('version','', @ischar);
      p.addParameter('data_dir',...
        simpletimeseries.parameters('value','data_dir'),...
        simpletimeseries.parameters('validation','data_dir'));
      p.parse(format,satname,start,varargin{:});
      %picking interface routine
      interface=str2func(['attitude.',format,'_filename']);
      %ensure cell satname
      satname=cells.scalar(satname,'set');
      %make room for outputs
      out=cell(1,numel(satname)*numel(start));
      %loop over vector inputs
      c=1;
      for s=1:numel(satname)
        for d=1:numel(start)
          %call interface routines
          out{c}=interface(satname{s},start(d),p.Results.version,p.Results.data_dir);
        end
      end
    end
    %import data according to format, satname, start and optionally data_dir:
    % - valid formats are according to the available <format>_filename routines
    % - valid satnames are according to simpletimeseries.translatesat
    % - start is datetime, from which the date in the filename is retrieved
    % - data_dir
    function obj=import(format,satname,varargin)
      % call parent
      sts=simpletimeseries.import(attitude.filename(format,satname,varargin{:}));
      %save loaded data in appropriate data type
      switch format
        case 'grace_l1b'
          %ensure quaternions are unitary
          sts=sts.assign_y(attitude.quat_unit(sts.y));
          args={'quat',sts};
          sat=simpletimeseries.translatesatname(satname);
          %TODO: check if this correct
          frame_from='srf';
          frame_to='crf';
          qsfirst=true;
        otherwise
          error(['Cannot handle format ''',format,'''.'])
      end
      if ~isempty(sts)
        obj=attitude(...
          args{:},...
          'satname',   sat,...
          'frame_from',frame_from,...
          'frame_to',  frame_to,...
          'qsfirst',   qsfirst...
        );
      else
        obj=[];
      end
    end
    %% math
    %sets norm of the quaternion equal to 1
    function q=quat_unit(q)
      scale=sqrt(sum(q.^2,2))*ones(1,4);
      if any(scale ~= 1)
        q=q./scale;
      end
    end
    %split scalar and vector parts
    function [s,v]=quat_split(q,qsfirst)
      if qsfirst
        s=q(:,1);
        v=q(:,2:4);
      else
        s=q(:,4);
        v=q(:,1:3);
      end
    end
    %joins scalar and vector parts
    function q=quat_join(s,v,qsfirst)
      if qsfirst
        q=[s,v];
      else
        q=[v,s];
      end
    end
    %scales the quaterion so it represents a rotation
    function q=quat_scale(q,qsfirst)
      [s,v]=attitude.quat_split(q,qsfirst);
      q=attitude.quat_join(s,v.*( abs(sin(acos(s))./sqrt(sum(v.^2,2)))*ones(1,3) ),qsfirst);
    end
    %compute the conjugate quaternion
    function q=quat_conj(q,qsfirst)
      [s,v]=attitude.quat_split(q,qsfirst);
      q=attitude.quat_join(s,-v,qsfirst);
    end
    %compute quaternion multiplication
    function q=quat_mult(q1,q2,qsfirst)
      %split
      [s1,v1]=attitude.quat_split(q1,qsfirst);
      [s2,v2]=attitude.quat_split(q2,qsfirst);
      %multiply the scale parts
      s=s1.*s2-sum(v1.*v2,2);
      %multiply the vector parts
      v=[s1.*v2(:,1),s1.*v2(:,2),s1.*v2(:,3)]+...
        [s2.*v1(:,1),s2.*v1(:,2),s2.*v1(:,3)]+...
        cross(v1,v2);
      %join
      q=attitude.quat_join(s,v,qsfirst);
    end
    %compute quaternion derivative
%     function dq=quat_diff(q)
%     end
    %% tests for the current object
    function out=test_parameters(field)
      switch lower(field)
      case 'format'
        out='grace_l1b';
      case 'satname'
        out='gracea';
      case 'start'
        out=datetime('2010-01-01');
      case 'version'
        out='02';
      otherwise
        error([mfilename,': unknown field ',field,'.'])
      end
    end
    function out=test(method)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      format= attitude.test_parameters('format');
      satname=attitude.test_parameters('satname');
      start=  attitude.test_parameters('start');
      switch(method)
        case 'all'
          for i={'filename','import','quat','ang','angr'}
            attitude.test(i{1});
          end
        case 'filename'
          out=attitude.filename(format,satname,start);
        case 'import'
          out=attitude.import(format,satname,start);
        case 'quat'
          out=attitude.test('import');
          figure
          out.quat.plot;
        case 'ang'
          out=attitude.test('import');
          figure
          out.ang.plot;
        case 'angr'
          out=attitude.test('import');
          figure
          out.angr.plot('zeromean',true);
        otherwise
          error(['Cannot handle test method ''',method,'''.'])
      end
    end
  end
  methods
    %% constructor
    function obj=attitude(varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      %parse the arguments with the names defined in attitude.data_type_list
      for j=1:numel(attitude.data_types)
        %shorter names
        dtn=attitude.data_types{j};
        %declare data types
        p.addParameter(dtn,[],(@(i) isa(i,'simpletimeseries')));
      end
      %declare parameters p
      [~,p,obj]=varargs.wrap('parser',p,'sinks',{obj},'sources',{attitude.parameters('obj')},varargin{:});
      %clean varargin
      varargin=cells.vararginclean(varargin,p.Parameters);
      % retrieve each data type
      for j=1:numel(attitude.data_types)
        %shorter names
        data_type=attitude.data_types{j};
        %skip if this data type is empty
        if ~isempty(p.Results.(data_type))
          %add new data type
          obj=obj.add_data_type(...
            data_type,p.Results.(data_type),...
            varargin{:}...
          );
        end
      end
      %propagate quaternion information to other data types
      obj=obj.update_ang;
      obj=obj.update_angr;
      %initialize internal records
      obj.initialized=true;
    end
    function obj=add_data_type(obj,data_type,data_value)
      %simplify things
      data_type=lower(data_type);
      %parse input
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('data_type' ,@(i) ischar(i) && cells.isincluded(attitude.data_types,i));
      p.addRequired('data_value',@(i) isa(i,'simpletimeseries'));
      % parse it
      p.parse(data_type,data_value);
      %sanity
      assert(isempty(obj.(data_type)),['data of type ''',data_type,''' has already been created. Use another method to append data.'])
      %propagate this data type
      obj.(data_type)=data_value;
    end
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      pn=[attitude.parameters('list');more_parameters(:)];
      for i=1:numel(pn)
        if isprop(obj,pn{i}) && isprop(obj_in,pn{i})
          obj.(pn{i})=obj_in.(pn{i});
        end
      end
      %propagate parameters of all non-empty data types
      for j=1:numel(attitude.data_types)
        %shorter names
        data_type=attitude.data_types{j};
        %sanity
        if xor(isempty(obj.(data_type)),isempty(obj_in.(data_type)))
          error([mfilename,': error propagating metadata of type ',data_type,': it does not exist in both objects.'])
        end
        %skip if data type is empty
        if ~isempty(obj.(data_type))
          obj.(data_type)=obj.(data_type).copy_metadata(obj_in.(data_type));
        end
      end
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      warning off MATLAB:structOnObject
      out=varargs(...
        structs.filter(struct(obj),[attitude.parameters('list');more_parameters(:)])...
      ).varargin;
      warning on MATLAB:structOnObject
    end
    function print(obj,tab)
      if ~exist('tab','var') || isempty(tab)
        tab=20;
      end
      disp(' --- Parameters --- ')
      for i=1:numel(attitude.parameters('list'))
        %shorter names
        p=attitude.parameters('value',i);
        disp([p,repmat(' ',1,tab-length(p)),' : ',str.show(obj.(p))])
      end
      d_list=attitude.data_types;
      for i=1:numel(d_list)
        %shorter names
        d=d_list{i};
        if ~isempty(obj.(d))
          disp([' --- ',d,' --- '])
          obj.(d).print
        end
      end
    end
    function msg(obj,in)
      if obj.verbose
        disp(in)
      end
    end
    %% properties of this object
    function obj=update_qsfirst(obj,qsfirst)
      %trivial call
      if (qsfirst==obj.qsfirst)
        return
      end
      %split scalar/vector parts
      [s,v]=attitude.quat_split(obj.quat.y,obj.qsfirst); 
      %swap scalar/vector parts
      if obj.qsfirst && ~qsfirst
        obj.quat=obj.quat.assign_y([v,s]); 
      elseif ~obj.qsfirst && qsfirst
        obj.quat=obj.quat.assign_y([s,v]); 
      end
    end
    function obj=set.qsfirst(obj,qsfirst)
      if ~isempty(obj.initialized) && obj.initialized %#ok<MCSUP>
        obj=obj.update_qsfirst(qsfirst);
      else
        obj.qsfirst=qsfirst;
      end
    end
    %% dependent properties of this object
    function s=get.scalar_part(obj)
      s=attitude.quat_split(obj.quat,obj.qsfirst);
    end
    function obj=set.scalar_part(obj,s)
      [~,v]=attitude.quat_split(obj.quat,obj.qsfirst);
      obj.quat=attitude.quat_join(s,v,obj.qsfirst);
    end
    function v=get.vector_part(obj)
      [~,v]=attitude.quat_split(obj.quat,obj.qsfirst);
    end
    function obj=set.vector_part(obj,v)
      s=attitude.quat_split(obj.quat,obj.qsfirst);
      obj.quat=attitude.quat_join(s,v,obj.qsfirst);
    end
    %% properties of data types
    function out=get_dtp(obj,propname)
      out=[];
      dt=attitude.data_types;
      for i=1:numel(dt)
        if ~isempty(obj.(dt{i}))
          out=obj.(dt{i}).(propname);
          break
        end
      end
      assert(~isempty(out),'all data types are empty.')
    end
    function obj=set_dtp(obj,propname,propvalue)
      odt=attitude.data_types;
      for i=1:numel(odt)
        if ~isempty(obj.(odt{i}))
          switch propname
            case {'time','t'}
              obj.(odt{i})=obj.(odt{i}).interp(t,'interp1_args',{'spline'});
            otherwise
              obj.(odt{i}).(propname)=propvalue;
          end
        end
      end
      obj.check_dtp(propname)
    end
    function check_dtp(obj,propname)
      odt=attitude.data_types;
      for i=2:numel(odt)
        if ~obj.isempty(odt{i})
          switch propname
            case {'time','t'}
              check=obj.(odt{1}).istequal(obj.(odt{i}));
            otherwise
              check=cells.isequal(...
                cells.set(obj.(odt{1}).(propname),'set'),...
                cells.set(obj.(odt{i}).(propname),'set'));
          end
          assert(check,...
            ['discrepancy in property ''',propname,''' between data types ''',...
            odt{1},''' and ''',odt{i},'''.']...
          )
        end
      end
    end
    function out=get.time(obj)
      out=obj.get_dtp('time');
    end
    function obj=set.time(obj,value)
      obj=obj.set_dtp('time',value);
    end
    %% euler angles
    function obj=update_ang(obj)
      assert(~isempty(obj.quat),'cannot compute data type ''ang'' because ''quat'' is empty')
      %getting euler angles
      y = SpinConv('QtoEA123', obj.update_qsfirst(false).quat.y);
      %building object
      obj.angi=simpletimeseries(obj.quat.t,y,...
        'format','datetime',...
        'units',{'deg', 'deg',  'deg'},...
        'labels', {'roll','pitch','yaw'},...
        'timesystem','gps',...
        'descriptor',obj.quat.descriptor...
      );
    end
    function obj=set.ang(obj,ang)
      %TODO: implement procedures that propagate to other datatypes (if needed)
      error('implementation needed')
    end
    function ang=get.ang(obj)
      assert(~isempty(obj.angi),'Need to call obj.update_ang before calling obj.ang')
      ang=obj.angi;
    end
    %% angular rates
    function obj=update_angr(obj)
      assert(~isempty(obj.quat),'cannot compute data type ''angr'' because ''quat'' is empty')
      %computing quaternion quantities
      qj=attitude.quat_conj(obj.quat.y,obj.qsfirst);
      %NOTICE: this calls simpledata.diff, but there is also simpletimeseries.deriv (unnecessary)
      dq=obj.quat.diff.y;
      %compute angular rates in frame_to
      y=2*attitude.quat_mult(qj,dq,obj.qsfirst);
      %%compute angular rates in frame_from
      %y=2*attitude.quat_mult(dq,qj,obj.qsfirst);
      if obj.qsfirst
          obj.msg([mfilename,': average of absolute value of scalar part of omega is ',...
              num2str(mean(abs(y(~isnan(y(:,1)),1))))])
          y=y(:,2:4);
      else
          obj.msg([mfilename,': average of absolute value of scalar part of omega is ',...
              num2str(mean(abs(y(~isnan(y(:,1)),4))))])
          y=y(:,1:3);
      end
      %building object
      obj.angri=simpletimeseries(obj.quat.t,y,...
        'format','datetime',...
        'units',{'deg/s', 'deg/s',  'deg/s'},...
        'labels', {'roll-rate','pitch-rate','yaw-rate'},...
        'timesystem','gps',...
        'descriptor',obj.quat.descriptor...
      );      
    end
    function obj=set.angr(obj,angr)
      %TODO: implement procedures that propagate to other datatypes (if needed)
      error('implementation needed')
    end
    function angr=get.angr(obj)
      assert(~isempty(obj.angi),'Need to call obj.update_angr before calling obj.angr')
      angr=obj.angri;
    end
    %% angular accelerations
    function out=get.anga(obj)
      if isempty(obj.anga)
        assert(~isempty(obj.quat),'cannot compute data type ''anga'' because ''quat'' is empty')
        %computing quaternion quantities
        q=obj.quat.y;
        qj=attitude.quat_conj(q,obj.qsfirst);
        dq=obj.quat.diff.y;
        %compute angular rates in frame_to
        y=2*attitude.quat_mult(qj,dq,obj.qsfirst);
        %%compute angular rates in frame_from
        %y=2*attitude.quat_mult(dq,qj,obj.qsfirst);
        if obj.qsfirst
            obj.msg([mfilename,': average of absolute value of scalar part of omega is ',...
                num2str(mean(abs(y(~isnan(y(:,1)),1))))])
            y=y(:,2:4);
        else
            obj.msg([mfilename,': average of absolute value of scalar part of omega is ',...
                num2str(mean(abs(y(~isnan(y(:,1)),4))))])
            y=y(:,1:3);
        end
        %building object
        obj.angr=simpletimeseries(obj.quat.t,y,...
          'format','datetime',...
          'units',{'deg/s', 'deg/s',  'deg/s'},...
          'labels', {'roll-rate','pitch-rate','yaw-rate'},...
          'timesystem','gps',...
          'descriptor',obj.quat.descriptor...
        );
      end
      out=obj.angr;
    end
    %% object properties
    function obj=set.satname(obj,in)
      obj.satname=simpletimeseries.translatesat(in);
    end
    function obj=set.frame_to(obj,in)
      obj.frame_to=simpletimeseries.translateframe(in);
    end
    function obj=set.frame_from(obj,in)
      obj.frame_from=simpletimeseries.translateframe(in);
    end
    %% general data_type scalar get method
    function out=get(obj,method,varargin)
      %get data types
      odt=attitude.data_types;
      %make room for outputs
      out=cell(size(odt));
      %loop over all data types
      for i=1:numel(odt)
        %check if this data type is not empty 
        if ~isempty(obj.(odt{i}))
          %check if this is a member
          if ismethod(obj.(odt{i}),method)
            %use the method on it, pass additional arguments
            out{i}=obj.(odt{i}).(method)(varargin{:});
          %check if this is a property
          elseif isprop(obj.(odt{i}),method)
            %call this member (varargin is ignored)
            out{i}=obj.(odt{i}).(method);
          end
        end
      end
      %remove dups and reduce to scalar if possible
      out=cells.scalar(cells.rm_duplicates(out),'get');
      %need to return something
      assert(~isempty(out),'all data types are empty.')
    end
    %% management
    function compatible(obj1,obj2,varargin)
      %This method checks if the objectives are referring to the same
      %type of data, i.e. the data length is not important.
      parameters=attitude.compatible_parameter_list;
      for i=1:numel(parameters)
        if ~isequal(obj1.(parameters{i}),obj2.(parameters{i}))
          error([mfilename,': discrepancy in parameter ',parameters{i},': ''',...
            obj1.(parameters{i}),''' ~= ''',obj2.(parameters{i}),'''.'])
        end
      end
      %check that all data type as compatible as well
      odt=attitude.data_types;
      for i=1:numel(odt)
        if ~isempty(obj1.(odt{i})) && ~isempty(obj2.(odt{i}))
          obj1.(odt{i}).compatible(obj2.(odt{i}),varargin{:})
        end
      end
    end
    %object obj1 will have the time domain of obj2 (interpolated if needed)
    function [obj1,obj2]=consolidate(obj1,obj2,varargin)
      %compatibility check
      compatible(obj1,obj2,varargin{:})
      %consolidate all data types
      counter=0;
      odt=attitude.data_types;
      for i=1:numel(odt)
        if ~isempty(obj1.(odt{i})) && ~isempty(obj2.(odt{i}))
          [obj1.(odt{i}),obj2.(odt{i})]=obj1.(odt{i}).consolidate(obj2.(odt{i}));
          counter=counter+1;
        end
      end
      if counter==0
        error([mfilename,': there were no common fields in the input objects.'])
      end
    end
    function out=isempty(obj,data_type)
      if ~exist('data_type','var') || isempty(data_type)
        odt=attitude.data_types;
        for i=1:numel(odt)
          if ~obj.isempty(odt{i});out=false;return;end
        end
        out=true;
      else
        out=isempty(obj.(data_type)) || all(obj.(data_type).y(:)==0);
      end
    end
    %% operator
    % uses a method from a superclass over all non-empty data types
    function obj=op(obj,operation,varargin)
      %operation counter
      counter=0;
      odt=attitude.data_types;
      %handle operations between two attitude sets
      if numel(varargin)==1 && isa(varargin{1},'attitude')
        err_msg='there are no common fields in the input objects';
        %shorter names
        obj1=varargin{1};
        %No need for consolidation, that is done inside the operation
        for i=1:numel(odt)
          %operate over all non-empty data types
          if ~isempty(obj.(odt{i})) && ~isempty(obj1.(odt{i}))
            %operate
            obj.(odt{i})=obj.(odt{i}).(operation)(obj1.(odt{i}));
            counter=counter+1;
          end
        end
      else
        err_msg='there are no non-empty fields in the input object';
        for i=1:numel(odt)
          %operate over all non-empty data types
          if ~isempty(obj.(odt{i}))
            if isprop(obj.(odt{i}),operation)
              %sanity
              if numel(varargin)>1
                error([mfilename,': when propagating data to field ',operation,...
                  ' can only handle one input argument, not ',num2str(numel(varargin)),'.'])
              end
              %propagate
              obj.(odt{i}).(operation)=varargin{1};
            else
              %operate
              obj.(odt{i})=obj.(odt{i}).(operation)(varargin{:});
            end
            counter=counter+1;
          end
        end
      end
      if counter==0
        error([mfilename,': ',err_msg,'.'])
      end
    end
    %% relative attitude (in the local frame)
    function rel=relative(obj1,obj2)
      %compute difference
      rel=obj1.op('minus',obj2);
    end
    %TODO: needs testing
    function out=periodic_stats(obj,period,varargin)
      % separate time series into segments
      [ts,idx]=segmentedfreqseries.time_segmented(obj.time,period,seconds(0));
      %get stats for all available stat types
      odt=attitude.data_types;
      for j=1:numel(odt)
        if ~isempty(obj.(odt{j}))
          % initialize
          s.msg=[mfilename,': cutting into segments data of type ''',odt{j},'''.'];s.n=numel(ts);
          clear tmp
          % propagate segments
          for i=1:numel(ts)
            %compute statistics
            tmp(i)=simpledata(...
              ts{i},...
              obj.(odt{j}).y(idx{i}(1):idx{i}(2),:),...
              'mask',obj.(odt{j}).mask(idx{i}(1):idx{i}(2),:),...
              varargin{:}...
            ).stats('struct'); %#ok<*AGROW>
            %inform
            s=simpledata.progress(s,i);
          end
          %propagate
          s_list=fields(tmp);
          for i=1:numel(s_list)
            stats.(odt{j}).(s_list{i})=transpose(reshape(...
              [tmp.(s_list{i})],...
              numel(tmp(1).(s_list{i})),...
              numel(ts)...
            ));
          end
        end
      end
      %build time domain
      t=datetime([],[],[]);
      for i=1:numel(ts)
        t(i)=mean(ts{i});
      end
      %build argument list
      o_list=fields(stats);
      args=cell(1,2*numel(o_list));
      %propagate all fields in statistics
      s_list=fields(tmp);
      for i=1:numel(s_list)
        %build argument list for this statistic
        for j=1:numel(o_list)
          args{2*j-1}=o_list{j};
          if size(stats.(o_list{j}).(s_list{i}),2)==1
            args{2*j  }=repmat(stats.(o_list{j}).(s_list{i}),1,attitude.data_type_list.(o_list{j}).size);
          elseif size(stats.(o_list{j}).(s_list{i}),2)==attitude.data_type_list.(o_list{j}).size
            args{2*j  }=stats.(o_list{j}).(s_list{i});
          else
            error([mfilename,': BUG TRAP: statistic with non-comformant number of columns. Debug needed!'])
          end
        end
        %build attitude object for this statistic
        out.(s_list{i})=attitude(t,args{:}).copy_metadata(obj);
      end

    end
  end
end
% 
% function euler=num_quaternion2euler(quat,scalar_first_flag,rotation_type)
% %computes the euler angles associated with a rotation represented
% %by the quaternion set <in>, which is assumed to be a list of
% %quaternions, where each quaternion is in one line (nr of columns is 4).
% 
%   %parameters
%   tol=1e-6;
%   ichk=0;
%   %sanity on inputs
%   if size(quat,2) ~= 4
%       error([mfilename,': input <quat> must have 4 columns.'])
%   end
%   %optionals
%   if ~exist('scalar_first_flag','var') || isempty(scalar_first_flag)
%       scalar_first_flag=false;
%   end
%   %branch on type of quaternion
%   if scalar_first_flag
%       euler=SpinCalc(['QtoEA',rotation_type],[quat(:,2:4),quat(:,1)],tol,ichk)*pi/180;
%   else
%       euler=SpinCalc(['QtoEA',rotation_type],quat,tol,ichk)*pi/180;
%   end
% end
% 

function OUTPUT = SpinConv(TYPES, INPUT, tol, ichk)
%SpinConv   Conversion from a rotation representation type to another
%
%   OUT = SpinConv(TYPES, IN, TOL, ICHK) converts a rotation representation
%   type (IN) to another (OUT). Supported conversion input/output types are
%   as follows:
%      1) Q      Rotation quaternions
%      2) EV     Euler vector and rotation angle (degrees)
%      3) DCM    Direction cosine matrix (a.k.a. rotation matrix)
%      4) EA###  Euler angles (12 possible sets) (degrees)
%   All representation types accepted as input (IN) and returned as output
%   (OUT) by SpinConv are meant to represent the rotation of a 3D
%   coordinate system (CS) relative to a rigid body or vector space
%   ("alias" transformation), rather than vice-versa ("alibi"
%   transformation).
%
%   OUT=SpinConv(TYPES,IN) is equivalent to OUT=SpinConv(TYPES,IN,10*eps,1)
%   OUT=SpinConv(TYPES,IN,TOL) is equiv. to OUT=SpinConv(TYPES,IN,TOL,1)
%
%   Input and output arguments:
%
%      TYPES - Single string value that specifies both the input type and
%            the desired output type. The allowed values are:
%
%               'DCMtoEA###'      'DCMtoEV'      'DCMtoQ'
%               'EA###toDCM'      'EA###toEV'    'EA###toQ'
%               'EVtoDCM'         'EVtoEA###'    'EVtoQ'
%               'QtoDCM'          'QtoEA###'     'QtoEV'
%               'EA###toEA###'
%
%            For cases that involve Euler angles, ### should be
%            replaced with the proper order desired. E.g., EA321
%            would be Z(yaw)-Y(pitch)-X(roll).
%
%      IN  - Array of N matrices or N vectors (N>0) corresponding to the
%            first entry in the TYPES string, formatted as follows:
%
%            DCM - (3×3×N) Array of rotation matrices. Each matrix R
%                  contains in its rows the versors of the rotated CS
%                  represented in the original CS, and in its columns the
%                  versors of the original CS represented in the rotated
%                  CS. This format is typically used when the column-vector
%                  convention is adopted: point coordinates are arranged in
%                  column vectors Vi, and the desired rotation is applied
%                  by pre-multiplying Vi by R (rotated Vi = R * Vi).
%            EA### - [psi,theta,phi] (N×3) row vector list containing, in
%                  each row, three Euler angles or Tait-Bryan angles.
%                  (degrees).
%            EV  - [m1,m2,m3,MU] (N×4) Row vector list containing, in each
%                  row, the components (m1, m2, m3) of an Euler rotation
%                  vector (represented in the original CS) and the Euler
%                  rotation angle about that vector (MU, in degrees).
%            Q   - [q1,q2,q3,q4] (N×4) Row vector list defining, in each
%                  row, a rotation quaternion. q4 = cos(MU/2), where MU is
%                  the Euler angle.
%
%      TOL - (Default value: TOL = 10 * eps) Tolerance value for deviations
%            from 1. Used to test determinant of rotation matrices or
%            length of unit vectors.
%      ICHK - (Default value: ICHK = 1) Flag controlling whether
%            near-singularity warnings are issued or not.
%            ICHK = 0 disables warnings.
%            ICHK = 1 enables them.
%      OUT - Array of N matrices or N vectors (N > 0) corresponding to the
%            second entry in the TYPES string, formatted as shown
%            above.
%
%   See also SpinCalc, degtorad, rad2deg.

% Version 2.2
% 2013 April 3
%
% Based on:
%    SpinCalc, Version 1.3 (MATLAB Central file #20696)
%    2009 June 30
% SpinCalc code by:
%    John Fuller
%    National Institute of Aerospace
%    Hampton, VA 23666
%    John.Fuller@nianet.org
% Debugged and optimized for speed by:
%    Paolo de Leva
%    University "Foro Italico"
%    Rome, Italy
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% Setting default values for missing input arguments
switch nargin
  case 2, tol = 10*eps; ichk = true;
  case 3, ichk = true;
  case 4
    if isequal(ichk, 0), ichk = false;
    else,                 ichk = true;
    end
  otherwise, narginchk(2, 4); % Allow 2 to 4 input arguments
end

% No TYPES string can be shorter than 4 or longer than 12 chars
len = length(TYPES);
if len>12 || len<4, error('Invalid entry for TYPES input string'); end

% Determine type of conversion from TYPES string
TYPES = upper(TYPES);
index = strfind(TYPES, 'TO');
TYPE.INPUT  = TYPES(1 : index-1);
TYPE.OUTPUT = TYPES(index+2 : len);
fields = {'INPUT', 'OUTPUT'}; % 1×2 cell
% Check validity of TYPES string, both for input and output
for f = 1:2
  IO = fields{f};
  type = TYPE.(IO);
  switch type
    case {'Q' 'EV' 'DCM'} % Valid TYPE
    otherwise
      % Check that TYPE is 'EA###'
      if length(type)~=5 || ~strcmp(type(1:2), 'EA')
        error('Invalid entry for TYPES input string')
      end
      TYPE.(IO) = 'EA';
      EAorder.(IO) = type(3:5);
      % Check that all characters in '###' are numbers between 1
      % and 3, and that no 2 consecutive characters are equal
      order31 = str2num(EAorder.(IO)'); %#ok<ST2NM> % 3×1 double
      if isempty(order31) || any ([order31<1; order31>3]) || ...
          order31(1)==order31(2) || ...
          order31(2)==order31(3)
        error('Invalid Euler angle order in TYPES string.')
      end
      % Type of EA sequence:
      %    1) Rotations about three distinct axes
      %    2) 1st and 3rd rotation about same axis
      if order31(1)==order31(3), EAtype.(IO) = 2;
      else                     , EAtype.(IO) = 1;
      end
  end
end

% Set N (number of rotations) and check INPUT size
s=size(INPUT);
switch numel(s)
  case 2; size1=s(1); size2=s(2);
  case 3; size1=s(1); size2=s(2); size3=s(3);
end
switch TYPE.INPUT
  case 'DCM' % (3×3×N) Direction cosine matrix
    N = size3;
    isnot_DCM = false;
    if ndims(INPUT)>3 || N==0 || size1~=3 || size2~=3
      error('Invalid INPUT size (INPUT must be 3×3×N for DCM type)')
    end
  case 'EA', v_length=3; Isize='N×3'; isnot_DCM=true;
  case 'Q',  v_length=4; Isize='N×4'; isnot_DCM=true;
  case 'EV', v_length=4; Isize='N×4'; isnot_DCM=true;
end
if isnot_DCM
  N = size1;
  if ~ismatrix(INPUT) || N==0 || size2~=v_length
    error(['Invalid INPUT size (INPUT must be ' ...
      Isize ' for ' TYPE.INPUT ' type)'])
  end
end

% Determine the quaternions that uniquely describe the rotation prescribed
% by INPUT. OUTPUT will be calculated in the second portion of the code
% from these quaternions.
switch TYPE.INPUT
  
  case 'DCM'
    % NOTE: Orthogonal matrixes may have determinant -1 or 1
    %       DCMs are special orthogonal matrices, with determinant 1
    improper  = false;
    DCM_not_1 = false;
    if N == 1
      % Computing deviation from orthogonality
      delta = INPUT * INPUT' - eye(3); % DCM*DCM' - I
      delta = delta(:); % 9×1 <-- 3×3
      % Checking determinant of DCM
      DET = det(INPUT);
      if DET<0, improper=true; end
      if ichk && abs(DET-1)>tol, DCM_not_1=true; end
      % Permuting INPUT
      INPUT = reshape(INPUT, [1 3 3]); % 1×3×3
    else
      % Computing deviation from orthogonality
      delta = multiprod(INPUT, multitransp(INPUT), [1 2]); % DCM*DCM'
      delta = bsxfun(@minus, delta, eye(3)); % Subtracting I
      delta = delta(:); % 9N×1 <-- 3×3×N
      % Checking determinant of DCMs
      DET = INPUT(1,1,:).*INPUT(2,2,:).*INPUT(3,3,:) -INPUT(1,1,:).*INPUT(2,3,:).*INPUT(3,2,:)...
        +INPUT(1,2,:).*INPUT(2,3,:).*INPUT(3,1,:) -INPUT(1,2,:).*INPUT(2,1,:).*INPUT(3,3,:)...
        +INPUT(1,3,:).*INPUT(2,1,:).*INPUT(3,2,:) -INPUT(1,3,:).*INPUT(2,2,:).*INPUT(3,1,:); % 1×1×N
      if any(DET<0), improper=true; end
      if ichk && any(abs(DET-1)>tol), DCM_not_1=true; end
      % Permuting INPUT
      INPUT = permute(INPUT, [3 1 2]); % N×3×3
    end
    % Issuing error messages or warnings
    if ichk && any(abs(delta)>tol)
      warning('Input DCM is not orthogonal.')
    end
    if improper, error('Improper input DCM'); end
    if DCM_not_1
      warning('Input DCM determinant off from 1 by more than tolerance.');
    end
    % Denominators for 4 distinct types of equivalent Q equations
    denom = [1 + INPUT(:,1,1) - INPUT(:,2,2) - INPUT(:,3,3),...
      1 - INPUT(:,1,1) + INPUT(:,2,2) - INPUT(:,3,3),...
      1 - INPUT(:,1,1) - INPUT(:,2,2) + INPUT(:,3,3),...
      1 + INPUT(:,1,1) + INPUT(:,2,2) + INPUT(:,3,3)];
    denom = 2 .* sqrt (denom); % N×4
    % Choosing for each DCM the equation which uses largest denominator
    [maxdenom, index] = max(denom, [], 2); % N×1
    clear delta DET denom
    Q = NaN(N,4); % N×4
    % EQUATION 1
    ii = (index==1); % (Logical vector) MAXDENOM==DENOM(:,1)
    if any(ii)
      Q(ii,:) = [                         0.25 .* maxdenom(ii,1),...
        (INPUT(ii,1,2)+INPUT(ii,2,1)) ./ maxdenom(ii,1),...
        (INPUT(ii,1,3)+INPUT(ii,3,1)) ./ maxdenom(ii,1),...
        (INPUT(ii,2,3)-INPUT(ii,3,2)) ./ maxdenom(ii,1)];
    end
    % EQUATION 2
    ii = (index==2); % (Logical vector) MAXDENOM==DENOM(:,2)
    if any(ii)
      Q(ii,:) = [(INPUT(ii,1,2)+INPUT(ii,2,1)) ./ maxdenom(ii,1),...
        0.25 .* maxdenom(ii,1),...
        (INPUT(ii,2,3)+INPUT(ii,3,2)) ./ maxdenom(ii,1),...
        (INPUT(ii,3,1)-INPUT(ii,1,3)) ./ maxdenom(ii,1)];
    end
    % EQUATION 3
    ii = (index==3); % (Logical vector) MAXDENOM==DENOM(:,3)
    if any(ii)
      Q(ii,:) = [(INPUT(ii,1,3)+INPUT(ii,3,1)) ./ maxdenom(ii,1),...
        (INPUT(ii,2,3)+INPUT(ii,3,2)) ./ maxdenom(ii,1),...
        0.25 .* maxdenom(ii,1),...
        (INPUT(ii,1,2)-INPUT(ii,2,1)) ./ maxdenom(ii,1)];
    end
    % EQUATION 4
    ii = (index==4); % (Logical vector) MAXDENOM==DENOM(:,4)
    if any(ii)
      Q(ii,:) = [(INPUT(ii,2,3)-INPUT(ii,3,2)) ./ maxdenom(ii,1),...
        (INPUT(ii,3,1)-INPUT(ii,1,3)) ./ maxdenom(ii,1),...
        (INPUT(ii,1,2)-INPUT(ii,2,1)) ./ maxdenom(ii,1),...
        0.25 .* maxdenom(ii)];
    end
    clear INPUT maxdenom index ii
    
  case 'EV'
    % Euler vector (EV) and angle MU in degrees
    EV = INPUT(:,1:3); % N×3
    halfMU = INPUT(:,4) * (pi/360); % (N×1) MU/2 in radians
    % Check that input m's constitute unit vector
    delta = sqrt(sum(EV.*EV, 2)) - 1; % N×1
    if any(abs(delta) > tol)
      error('(At least one of the) input Euler vector(s) is not a unit vector')
    end
    % Quaternion
    SIN = sin(halfMU); % (N×1)
    Q = [EV(:,1).*SIN, EV(:,2).*SIN, EV(:,3).*SIN, cos(halfMU)];
    clear EV delta halfMU SIN
    
  case 'EA'
    % Identify singularities (2nd Euler angle out of range)
    theta = INPUT(:, 2); % N×1
    if EAtype.INPUT == 1
      % Type 1 rotation (rotations about three distinct axes)
      if any(abs(theta)>=90)
        error('Second input Euler angle(s) outside -90 to 90 degree range')
      elseif ichk && any(abs(theta)>88)
        warning(['Second input Euler angle(s) near a '...
          'singularity (-90 or 90 degrees).'])
      end
    else
      % Type 2 rotation (1st and 3rd rotation about same axis)
      if any(theta<=0 | theta>=180)
        error('Second input Euler angle(s) outside 0 to 180 degree range')
      elseif ichk && any(theta<2 | theta>178)
        warning(['Second input Euler angle(s) near a '...
          'singularity (0 or 180 degrees).'])
      end
    end
    % Half angles in radians
    HALF = INPUT * (pi/360); % N×3
    Hpsi   = HALF(:,1); % N×1
    Htheta = HALF(:,2); % N×1
    Hphi   = HALF(:,3); % N×1
    % Pre-calculate cosines and sines of the half-angles for conversion.
    c1=cos(Hpsi); c2=cos(Htheta); c3=cos(Hphi);
    s1=sin(Hpsi); s2=sin(Htheta); s3=sin(Hphi);
    c13 =cos(Hpsi+Hphi);  s13 =sin(Hpsi+Hphi);
    c1_3=cos(Hpsi-Hphi);  s1_3=sin(Hpsi-Hphi);
    c3_1=cos(Hphi-Hpsi);  s3_1=sin(Hphi-Hpsi);
    clear HALF Hpsi Htheta Hphi
    switch EAorder.INPUT
      case '121', Q=[c2.*s13,  s2.*c1_3, s2.*s1_3, c2.*c13];
      case '232', Q=[s2.*s1_3, c2.*s13,  s2.*c1_3, c2.*c13];
      case '313', Q=[s2.*c1_3, s2.*s1_3, c2.*s13,  c2.*c13];
      case '131', Q=[c2.*s13,  s2.*s3_1, s2.*c3_1, c2.*c13];
      case '212', Q=[s2.*c3_1, c2.*s13,  s2.*s3_1, c2.*c13];
      case '323', Q=[s2.*s3_1, s2.*c3_1, c2.*s13,  c2.*c13];
      case '123', Q=[s1.*c2.*c3+c1.*s2.*s3, c1.*s2.*c3-s1.*c2.*s3, c1.*c2.*s3+s1.*s2.*c3, c1.*c2.*c3-s1.*s2.*s3];
      case '231', Q=[c1.*c2.*s3+s1.*s2.*c3, s1.*c2.*c3+c1.*s2.*s3, c1.*s2.*c3-s1.*c2.*s3, c1.*c2.*c3-s1.*s2.*s3];
      case '312', Q=[c1.*s2.*c3-s1.*c2.*s3, c1.*c2.*s3+s1.*s2.*c3, s1.*c2.*c3+c1.*s2.*s3, c1.*c2.*c3-s1.*s2.*s3];
      case '132', Q=[s1.*c2.*c3-c1.*s2.*s3, c1.*c2.*s3-s1.*s2.*c3, c1.*s2.*c3+s1.*c2.*s3, c1.*c2.*c3+s1.*s2.*s3];
      case '213', Q=[c1.*s2.*c3+s1.*c2.*s3, s1.*c2.*c3-c1.*s2.*s3, c1.*c2.*s3-s1.*s2.*c3, c1.*c2.*c3+s1.*s2.*s3];
      case '321', Q=[c1.*c2.*s3-s1.*s2.*c3, c1.*s2.*c3+s1.*c2.*s3, s1.*c2.*c3-c1.*s2.*s3, c1.*c2.*c3+s1.*s2.*s3];
      otherwise
        error('Invalid input Euler angle order (TYPES string)');
    end
    clear c1 c2 c3 s1 s2 s3 c13 s13 c1_3 s1_3 c3_1 s3_1
    
  case 'Q'
    if ichk && any(abs(sqrt(sum(INPUT.*INPUT, 2)) - 1) > tol)
      warning('(At least one of the) Input quaternion(s) is not a unit vector')
    end
    Q = INPUT;
end
clear TYPE.INPUT EAorder.INPUT

% Normalize quaternion(s) in case of deviation from unity.
% User has already been warned of deviation.
Qnorms = sqrt(sum(Q.*Q,2));
Q = [Q(:,1)./Qnorms, Q(:,2)./Qnorms, Q(:,3)./Qnorms, Q(:,4)./Qnorms]; % N×4

switch TYPE.OUTPUT
  
  case 'DCM'
    Q  = reshape(Q', [1 4 N]); % (1×4×N)
    SQ = Q.^2;
    OUTPUT = [   SQ(1,1,:)-SQ(1,2,:)-SQ(1,3,:)+SQ(1,4,:),  2.*(Q(1,1,:).*Q(1,2,:) +Q(1,3,:).*Q(1,4,:)), 2.*(Q(1,1,:).*Q(1,3,:) -Q(1,2,:).*Q(1,4,:));
      2.*(Q(1,1,:).*Q(1,2,:) -Q(1,3,:).*Q(1,4,:)),   -SQ(1,1,:)+SQ(1,2,:)-SQ(1,3,:)+SQ(1,4,:),  2.*(Q(1,2,:).*Q(1,3,:) +Q(1,1,:).*Q(1,4,:));
      2.*(Q(1,1,:).*Q(1,3,:) +Q(1,2,:).*Q(1,4,:)), 2.*(Q(1,2,:).*Q(1,3,:) -Q(1,1,:).*Q(1,4,:)),   -SQ(1,1,:)-SQ(1,2,:)+SQ(1,3,:)+SQ(1,4,:)];
    
  case 'EV'
    % Angle MU in radians and sine of MU/2
    halfMUrad = atan2( sqrt(sum(Q(:,1:3).*Q(:,1:3),2)), Q(:,4) ); % N×1
    SIN = sin(halfMUrad); % N×1
    index = (SIN==0); % (N×1) Logical index
    if any(index)
      % Initializing
      OUTPUT = zeros(N,4);
      % Singular cases (MU is zero degrees)
      OUTPUT(index, 1) = 1;
      % Non-singular cases
      SIN = SIN(~index, 1);
      OUTPUT(~index, :) = [Q(~index,1) ./ SIN, ...
        Q(~index,2) ./ SIN, ...
        Q(~index,3) ./ SIN, ...
        halfMUrad .* (360/pi)];
    else
      % Non-singular cases
      OUTPUT = [Q(:,1)./SIN, Q(:,2)./SIN, Q(:,3)./SIN, halfMUrad.*(360/pi)];
    end
    % MU greater than 180 degrees
    index = (OUTPUT(:,4) > 180); % (N×1) Logical index
    OUTPUT(index, :) = [-OUTPUT(index,1:3), 360-OUTPUT(index,4)];
    
  case 'EA'
    SQ = Q.^2;
    switch EAorder.OUTPUT
      case '121'
        OUTPUT = [atan2(Q(:,1).*Q(:,2) +Q(:,3).*Q(:,4), Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)),...
          acos(SQ(:,4)+SQ(:,1)-SQ(:,2)-SQ(:,3)),...
          atan2(Q(:,1).*Q(:,2) -Q(:,3).*Q(:,4), Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4))];
      case '232'
        OUTPUT = [atan2(Q(:,1).*Q(:,4) +Q(:,2).*Q(:,3), Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)),...
          acos(SQ(:,4)-SQ(:,1)+SQ(:,2)-SQ(:,3)),...
          atan2(Q(:,2).*Q(:,3) -Q(:,1).*Q(:,4), Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4))];
      case '313'
        OUTPUT = [atan2(Q(:,1).*Q(:,3) +Q(:,2).*Q(:,4), Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)),...
          acos(SQ(:,4)-SQ(:,1)-SQ(:,2)+SQ(:,3)),...
          atan2(Q(:,1).*Q(:,3) -Q(:,2).*Q(:,4), Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3))];
      case '131'
        OUTPUT = [atan2(Q(:,1).*Q(:,3) -Q(:,2).*Q(:,4), Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)),...
          acos(SQ(:,4)+SQ(:,1)-SQ(:,2)-SQ(:,3)),...
          atan2(Q(:,1).*Q(:,3) +Q(:,2).*Q(:,4), Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2))];
      case '212'
        OUTPUT = [atan2(Q(:,1).*Q(:,2) -Q(:,3).*Q(:,4), Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)),...
          acos(SQ(:,4)-SQ(:,1)+SQ(:,2)-SQ(:,3)),...
          atan2(Q(:,1).*Q(:,2) +Q(:,3).*Q(:,4), Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3))];
      case '323'
        OUTPUT = [atan2(Q(:,2).*Q(:,3) -Q(:,1).*Q(:,4), Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)),...
          acos(SQ(:,4)-SQ(:,1)-SQ(:,2)+SQ(:,3)),...
          atan2(Q(:,1).*Q(:,4) +Q(:,2).*Q(:,3), Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3))];
      case '123'
        OUTPUT = [atan2(2.*(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3)), SQ(:,4)-SQ(:,1)-SQ(:,2)+SQ(:,3)),...
          asin(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4))),...
          atan2(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)), SQ(:,4)+SQ(:,1)-SQ(:,2)-SQ(:,3))];
      case '231'
        OUTPUT = [atan2(2.*(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3)), SQ(:,4)+SQ(:,1)-SQ(:,2)-SQ(:,3)),...
          asin(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4))),...
          atan2(2.*(Q(:,1).*Q(:,4)-Q(:,3).*Q(:,2)), SQ(:,4)-SQ(:,1)+SQ(:,2)-SQ(:,3))];
      case '312'
        OUTPUT = [atan2(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2)), SQ(:,4)-SQ(:,1)+SQ(:,2)-SQ(:,3)),...
          asin(2.*(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3))),...
          atan2(2.*(Q(:,2).*Q(:,4)-Q(:,3).*Q(:,1)), SQ(:,4)-SQ(:,1)-SQ(:,2)+SQ(:,3))];
      case '132'
        OUTPUT = [atan2(2.*(Q(:,1).*Q(:,4)+Q(:,2).*Q(:,3)), SQ(:,4)-SQ(:,1)+SQ(:,2)-SQ(:,3)),...
          asin(2.*(Q(:,3).*Q(:,4)-Q(:,1).*Q(:,2))),...
          atan2(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)), SQ(:,4)+SQ(:,1)-SQ(:,2)-SQ(:,3))];
      case '213'
        OUTPUT = [atan2(2.*(Q(:,1).*Q(:,3)+Q(:,2).*Q(:,4)), SQ(:,4)-SQ(:,1)-SQ(:,2)+SQ(:,3)),...
          asin(2.*(Q(:,1).*Q(:,4)-Q(:,2).*Q(:,3))),...
          atan2(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)), SQ(:,4)-SQ(:,1)+SQ(:,2)-SQ(:,3))];
      case '321'
        OUTPUT = [atan2(2.*(Q(:,1).*Q(:,2)+Q(:,3).*Q(:,4)), SQ(:,4)+SQ(:,1)-SQ(:,2)-SQ(:,3)),...
          asin(2.*(Q(:,2).*Q(:,4)-Q(:,1).*Q(:,3))),...
          atan2(2.*(Q(:,1).*Q(:,4)+Q(:,3).*Q(:,2)), SQ(:,4)-SQ(:,1)-SQ(:,2)+SQ(:,3))];
      otherwise
        error('Invalid output Euler angle order (TYPES string).');
    end
    OUTPUT = OUTPUT * (180/pi); % (N×3) Euler angles in degrees
    theta  = OUTPUT(:,2);       % (N×1) Angle THETA in degrees
    % Check OUTPUT
    if any(~isreal( OUTPUT(:) ))
      error('SpinConv:Unreal', ...
        ['Unreal Euler output. Input resides too close to singularity.\n' ...
        'Please choose different output type.'])
    end
    % Type 1 rotation (rotations about three distinct axes)
    % THETA is computed using ASIN and ranges from -90 to 90 degrees
    if EAtype.OUTPUT == 1
      singularities = abs(theta) > 89.9; % (N×1) Logical index
      if any(singularities)
        firstsing = find(singularities, 1); % (1×1)
        error(['Input rotation # %s resides too close to Type 1 Euler singularity.\n' ...
          'Type 1 Euler singularity occurs when second angle is -90 or 90 degrees.\n' ...
          'Please choose different output type.'], num2str(firstsing));
      end
      % Type 2 rotation (1st and 3rd rotation about same axis)
      % THETA is computed using ACOS and ranges from 0 to 180 degrees
    else
      singularities = theta<0.1 | theta>179.9; % (N×1) Logical index
      if any(singularities)
        firstsing = find(singularities, 1); % (1×1)
        error(['Input rotation # %s resides too close to Type 2 Euler singularity.\n' ...
          'Type 2 Euler singularity occurs when second angle is 0 or 180 degrees.\n' ...
          'Please choose different output type.'], num2str(firstsing));
      end
    end
    
  case 'Q'
    OUTPUT = Q;
end
end



