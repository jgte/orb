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
    %import data according to format, satname, start and optionally data_dir:
    % - valid formats are according to swith loop in this routine
    % - valid satnames are according to simpletimeseries.translatesat
    % - start is datetime, from which the date in the filename is retrieved
    % - data_dir
    function obj=import(format,satname,start,varargin)
      %save loaded data in appropriate data type
      switch format
        case 'grace_l1b'
          %get datafile
          datafile=grace.grace_l1b_filename('SCA1B',satname,start,varargin{:});
          %load the data
          sts=simpletimeseries.import(datafile,'format','SCA1B');
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
        error(['unknown field ',field,'.'])
      end
    end
    function out=test(method)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      format= attitude.test_parameters('format');
      satname=attitude.test_parameters('satname');
      start=  attitude.test_parameters('start');
      test_list={'import','quat','ang','angr'};
      switch(method)
        case 'all'
          for i=1:numel(test_list)
            out{i}=attitude.test(test_list{i});
          end
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
      p=machinery.inputParser;
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
      p=machinery.inputParser;
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
          error(['error propagating metadata of type ',data_type,': it does not exist in both objects.'])
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
    function out=varargin(obj,more_parameters)
      out=varargs(obj.metadata(more_parameters)).varargin;
    end
    %% info methods
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
      obj.angi=ang;
    end
    function ang=get.ang(obj)
      if isempty(obj.angi)
        ang=[];
      else
        ang=obj.angi;
      end
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
          obj.msg([' average of absolute value of scalar part of omega is ',...
              num2str(mean(abs(y(~isnan(y(:,1)),1))))])
          y=y(:,2:4);
      else
          obj.msg([' average of absolute value of scalar part of omega is ',...
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
      obj.angri=angr;
    end
    function angr=get.angr(obj)
      if isempty(obj.angri)
        angr=[];
      else
        angr=obj.angri;
      end
    end
    %% angular accelerations
    function obj=update_anga(obj)
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
          obj.msg([' average of absolute value of scalar part of omega is ',...
              num2str(mean(abs(y(~isnan(y(:,1)),1))))])
          y=y(:,2:4);
      else
          obj.msg([' average of absolute value of scalar part of omega is ',...
              num2str(mean(abs(y(~isnan(y(:,1)),4))))])
          y=y(:,1:3);
      end
      %building object
      obj.angai=simpletimeseries(obj.quat.t,y,...
        'format','datetime',...
        'units',{'deg/s', 'deg/s',  'deg/s'},...
        'labels', {'roll-rate','pitch-rate','yaw-rate'},...
        'timesystem','gps',...
        'descriptor',obj.quat.descriptor...
      );
    end
    function obj=set.anga(obj,anga)
      obj.angai=anga;
    end
    function anga=get.anga(obj)
      if isempty(obj.angai)
        anga=[];
      else
        anga=obj.angai;
      end
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
          error(['discrepancy in parameter ',parameters{i},': ''',...
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
        error('there were no common fields in the input objects.')
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
    %% multiple object manipulation
    function[obj1,idx1,idx2]=append(obj1,obj2,varargin)
      %append all data types
      for i=1:numel(attitude.data_types)
        %simplify things
        data_type=lower(attitude.data_types{i});
        if obj1.isempty(data_type) || obj2.isempty(data_type)
          continue
        end
        %call upstream method
        [...
          obj1.(data_type),...
          idx1.(data_type),...
          idx2.(data_type)...
        ]=obj1.(data_type).append(...
          obj2.(data_type),...
          varargin{:}...
        );
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
                error(['when propagating data to field ',operation,...
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
        error(err_msg)
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
          s.msg=['cutting into segments data of type ''',odt{j},'''.'];s.n=numel(ts);
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
            error('BUG TRAP: statistic with non-comformant number of columns. Debug needed!')
          end
        end
        %build attitude object for this statistic
        out.(s_list{i})=attitude(t,args{:},obj.varargin{:});
      end

    end
  end
end

%% SpinConv
% https://nl.mathworks.com/matlabcentral/fileexchange/41562-spinconv
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
%% Multiple matrix multiplications, with array expansion enabled
% https://nl.mathworks.com/matlabcentral/fileexchange/8773-multiple-matrix-multiplications-with-array-expansion-enabled
function c = multiprod(a, b, idA, idB)
%MULTIPROD  Multiplying 1-D or 2-D subarrays contained in two N-D arrays.
%   C = MULTIPROD(A,B) is equivalent  to C = MULTIPROD(A,B,[1 2],[1 2])
%   C = MULTIPROD(A,B,[D1 D2]) is eq. to C = MULTIPROD(A,B,[D1 D2],[D1 D2])
%   C = MULTIPROD(A,B,D1) is equival. to C = MULTIPROD(A,B,D1,D1)
%
%   MULTIPROD performs multiple matrix products, with array expansion (AX)
%   enabled. Its first two arguments A and B are "block arrays" of any
%   size, containing one or more 1-D or 2-D subarrays, called "blocks" (*).
%   For instance, a 5�6�3 array may be viewed as an array containing five
%   6�3 blocks. In this case, its size is denoted by 5�(6�3). The 1 or 2
%   adjacent dimensions along which the blocks are contained are called the
%   "internal dimensions" (IDs) of the array (�).
%
%   1) 2-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], [DB1 DB2]) contains the products
%         of the P�Q matrices in A by the R�S matrices in B. [DA1 DA2] are
%         the IDs of A; [DB1 DB2] are the IDs of B.
%
%   2) 2-D by 1-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, [DA1 DA2], DB1) contains the products of the
%         P�Q matrices in A by the R-element vectors in B. The latter are
%         considered to be R�1 matrices. [DA1 DA2] are the IDs of A; DB1 is
%         the ID of B.
%
%   3) 1-D by 2-D BLOCK(S) (*)
%         C = MULTIPROD(A, B, DA1, [DB1 DB2]) contains the products of the 
%         Q-element vectors in A by the R�S matrices in B. The vectors in A
%         are considered to be 1�Q matrices. DA1 is the ID of A; [DB1 DB2]
%         are the IDs of B.
%
%   4) 1-D BY 1-D BLOCK(S) (*)
%      (a) If either SIZE(A, DA1) == 1 or SIZE(B, DB1) == 1, or both,
%             C = MULTIPROD(A, B, DA1, DB1) returns products of scalars by 
%             vectors, or vectors by scalars or scalars by scalars.
%      (b) If SIZE(A, DA1) == SIZE(B, DB1), 
%             C = MULTIPROD(A, B, [0 DA1], [DB1 0]) or 
%             C = MULTIPROD(A, B, DA1, DB1) virtually turns the vectors
%             contained in A and B into 1�P and P�1 matrices, respectively,
%             then returns their products, similar to scalar products.
%             Namely, C = DOT2(A, B, DA1, DB1) is equivalent to 
%             C = MULTIPROD(CONJ(A), B, [0 DA1], [DB1 0]).
%      (c) Without limitations on the length of the vectors in A and B,
%             C = MULTIPROD(A, B, [DA1 0], [0 DB1]) turns the vectors
%             contained in A and B into P�1 and 1�Q matrices, respectively,
%             then returns their products, similar to outer products.
%             Namely, C = OUTER(A, B, DA1, DB1) is equivalent to
%             C = MULTIPROD(CONJ(A), B, [DA1 0], [0 DB1]).
%
%   Common constraints for all syntaxes:
%      The external dimensions of A and B must either be identical or 
%      compatible with AX rules. The internal dimensions of each block
%      array must be adjacent (DA2 == DA1 + 1 and DB2 == DB1 + 1 are
%      required). DA1 and DB1 are allowed to be larger than NDIMS(A) and
%      NDIMS(B). In syntaxes 1, 2, and 3, Q == R is required, unless the
%      blocks in A or B are scalars. 
%
%   Array expansion (AX):
%      AX is a powerful generalization to N-D of the concept of scalar
%      expansion. Indeed, A and B may be scalars, vectors, matrices or
%      multi-dimensional arrays. Scalar expansion is the virtual
%      replication or annihilation of a scalar which allows you to combine
%      it, element by element, with an array X of any size (e.g. X+10,
%      X*10, or []-10). Similarly, in MULTIPROD, the purpose of AX is to
%      automatically match the size of the external dimensions (EDs) of A
%      and B, so that block-by-block products can be performed. ED matching
%      is achieved by means of a dimension shift followed by a singleton
%      expansion:
%      1) Dimension shift (see SHIFTDIM).
%            Whenever DA1 ~= DB1, a shift is applied to impose DA1 == DB1.
%            If DA1 > DB1, B is shifted to the right by DA1 - DB1 steps.
%            If DB1 > DA1, A is shifted to the right by DB1 - DA1 steps.
%      2) Singleton expansion (SX).
%            Whenever an ED of either A or B is singleton and the
%            corresponding ED of the other array is not, the mismatch is
%            fixed by virtually replicating the array (or diminishing it to
%            length 0) along that dimension.
% 
%   MULTIPROD is a generalization for N-D arrays of the matrix
%   multiplication function MTIMES, with AX enabled. Vector inner, outer,
%   and cross products generalized for N-D arrays and with AX enabled are
%   performed by DOT2, OUTER, and CROSS2 (MATLAB Central, file #8782).
%   Elementwise multiplications (see TIMES) and other elementwise binary
%   operations with AX enabled are performed by BAXFUN (MATLAB Central,
%   file #23084). Together, these functions make up the �ARRAYLAB toolbox�.
%
%   Input and output format:
%      The size of the EDs of C is determined by AX. Block size is
%      determined as follows, for each of the above-listed syntaxes:
%      1) C contains P�S matrices along IDs MAX([DA1 DA2], [DB1 DB2]).
%      2) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A         P�Q  (2-D)     [DA1 DA2]
%         B         R    (1-D)     DB1
%         C (a)     P    (1-D)     MAX(DA1, DB1)
%         C (b)     P�Q  (2-D)     MAX([DA1 DA2], [DB1 DB1+1])
%         ----------------------------------------------------
%         (a) The 1-D blocks in B are not scalars (R > 1).
%         (b) The 1-D blocks in B are scalars (R = 1).
%      3) Array     Block size     ID(s)
%         ----------------------------------------------------
%         A           Q  (1-D)     DA1
%         B         R�S  (2-D)     [DB1 DB2]
%         C (a)       S  (1-D)     MAX(DA1, DB1)
%         C (b)     R�S  (2-D)     MAX([DA1 DA1+1], [DB1 DB2])
%         ----------------------------------------------------
%         (a) The 1-D blocks in A are not scalars (Q > 1).
%         (b) The 1-D blocks in A are scalars (Q = 1).
%      4)     Array     Block size         ID(s)
%         --------------------------------------------------------------
%         (a) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         MAX(P,Q) (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (b) A         P        (1-D)     DA1
%             B         P        (1-D)     DB1
%             C         1        (1-D)     MAX(DA1, DB1)
%         --------------------------------------------------------------
%         (c) A         P        (1-D)     DA1
%             B         Q        (1-D)     DB1
%             C         P�Q      (2-D)     MAX([DA1 DA1+1], [DB1 DB1+1])
%         --------------------------------------------------------------
%
%   Terminological notes:
%   (*) 1-D and 2-D blocks are generically referred to as "vectors" and 
%       "matrices", respectively. However, both may be also called
%       �scalars� if they have a single element. Moreover, matrices with a
%       single row or column (e.g. 1�3 or 3�1) may be also called �row
%       vectors� or �column vectors�.
%   (�) Not to be confused with the "inner dimensions" of the two matrices
%       involved in a product X * Y, defined as the 2nd dimension of X and
%       the 1st of Y (DA2 and DB1 in syntaxes 1, 2, 3).
%
%   Examples:
%    1) If  A is .................... a 5�(6�3)�2 array,
%       and B is .................... a 5�(3�4)�2 array,
%       C = MULTIPROD(A, B, [2 3]) is a 5�(6�4)�2 array.
%
%       A single matrix A pre-multiplies each matrix in B
%       If  A is ........................... a (1�3)    single matrix,
%       and B is ........................... a 10�(3�4) 3-D array,
%       C = MULTIPROD(A, B, [1 2], [3 4]) is a 10�(1�4) 3-D array.
%
%       Each matrix in A pre-multiplies each matrix in B (all possible
%       combinations)
%       If  A is .................... a (6�3)�5   array,
%       and B is .................... a (3�4)�1�2 array,
%       C = MULTIPROD(A, B, [1 2]) is a (6�4)�5�2 array.
%
%   2a) If  A is ........................... a 5�(6�3)�2 4-D array,
%       and B is ........................... a 5�(3)�2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5�(6)�2   3-D array.
%
%   2b) If  A is ........................... a 5�(6�3)�2 4-D array,
%       and B is ........................... a 5�(1)�2   3-D array,
%       C = MULTIPROD(A, B, [2 3], [2]) is   a 5�(6�3)�2 4-D array.
%
%   4a) If both A and B are .................. 5�(6)�2   3-D arrays,
%       C = MULTIPROD(A, B, 2) is .......... a 5�(1)�2   3-D array, while
%   4b) C = MULTIPROD(A, B, [2 0], [0 2]) is a 5�(6�6)�2 4-D array
%
%   See also DOT2, OUTER, CROSS2, BAXFUN, MULTITRANSP.
% $ Version: 2.1 $
% CODE      by:            Paolo de Leva
%                          (Univ. of Rome, Foro Italico, IT)    2009 Jan 24
%           optimized by:  Paolo de Leva
%                          Jinhui Bai (Georgetown Univ., D.C.)  2009 Jan 24
% COMMENTS  by:            Paolo de Leva                        2009 Feb 24
% OUTPUT    tested by:     Paolo de Leva                        2009 Feb 24
% -------------------------------------------------------------------------
 narginchk(2, 4) ; % Allow 2 to 4 input arguments
switch nargin % Setting IDA and/or IDB
    case 2, idA = [1 2]; idB = [1 2];
    case 3, idB = idA;
end
% ESC 1 - Special simple case (both A and B are 2D), solved using C = A * B
     if ismatrix(a) && ismatrix(b) && ...
         isequal(idA,[1 2]) && isequal(idB,[1 2])
         c = a * b; return
     end
% MAIN 0 - Checking and evaluating array size, block size, and IDs
     sizeA0 = size(a);
     sizeB0 = size(b);
     [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
     squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
                                           sizeval(idA,idB, sizeA0,sizeB0);
% MAIN 1 - Applying dimension shift (first step of AX) and 
%          turning both A and B into arrays of either 1-D or 2-D blocks
     if sizeisnew(1), a = reshape(a, sizeA); end    
     if sizeisnew(2), b = reshape(b, sizeB); end
% MAIN 2 - Performing products with or without SX (second step of AX)
     if squashOK % SQUASH + MTIMES (fastest engine)
         c = squash2D_mtimes(a,b, idA,idB, sizeA,sizeB, squashOK); 
     elseif timesOK % TIMES (preferred w.r. to SX + TIMES)
         if sumOK, c = sum(a .* b, sumOK);
         else    , c =     a .* b; end
     elseif sxtimesOK % SX + TIMES
         if sumOK, c = sum(bsxfun(@times, a, b), sumOK);
         else    , c =     bsxfun(@times, a, b); end
     elseif mtimesOK % MTIMES (rarely used)
         c = a * b;
     end
% MAIN 3 - Reshaping C (by inserting or removing singleton dimensions)
     [sizeC, sizeCisnew] = adjustsize(size(c), shiftC, false, delC, false);
     if sizeCisnew, c = reshape(c, sizeC); end
end
function c = squash2D_mtimes(a, b, idA, idB, sizeA, sizeB, squashOK)
% SQUASH2D_MTIMES  Multiproduct with single-block expansion (SBX).
%    Actually, no expansion is performed. The multi-block array is
%    rearranged from N-D to 2-D, then MTIMES is applied, and eventually the
%    result is rearranged back to N-D. No additional memory is required.
%    One and only one of the two arrays must be single-block, and its IDs
%    must be [1 2] (MAIN 1 removes leading singletons). Both arrays
%    must contain 2-D blocks (MAIN 1 expands 1-D blocks to 2-D).
    if squashOK == 1 % A is multi-block, B is single-block (squashing A)
        % STEP 1 - Moving IDA(2) to last dimension
        nd = length(sizeA);
        d2 = idA(2);    
        order = [1:(d2-1) (d2+1):nd d2]; % Partial shifting
        a = permute(a, order); % ...�Q
        % STEP 2 - Squashing A from N-D to 2-D  
        q = sizeB(1);
        s = sizeB(2);
        lengthorder = length(order);
        collapsedsize = sizeA(order(1:lengthorder-1)); 
        n = prod(collapsedsize);
        a = reshape(a, [n, q]); % N�Q    
        fullsize = [collapsedsize s]; % Size to reshape C back to N-D
    else % B is multi-block, A is single-block (squashing B)
        % STEP 1 - Moving IDB(1) to first dimension
        nd = length(sizeB);
        d1 = idB(1);    
        order = [d1 1:(d1-1) (d1+1):nd]; % Partial shifting
        b = permute(b, order); % Q�...
        % STEP 2 - Squashing B from N-D to 2-D  
        p = sizeA(1);
        q = sizeA(2);
        lengthorder = length(order);
        collapsedsize = sizeB(order(2:lengthorder)); 
        n = prod(collapsedsize);
        b = reshape(b, [q, n]); % Q�N
        fullsize = [p collapsedsize]; % Size to reshape C back to N-D
    end
    % FINAL STEPS - Multiplication, reshape to N-D, inverse permutation
    invorder(order) = 1 : lengthorder;
    c = permute (reshape(a*b, fullsize), invorder);
end
function [sizeA, sizeB, shiftC, delC, sizeisnew, idA, idB, ...
          squashOK, sxtimesOK, timesOK, mtimesOK, sumOK] = ...
                                          sizeval(idA0,idB0, sizeA0,sizeB0)
%SIZEVAL   Evaluation of array size, block size, and IDs
%    Possible values for IDA and IDB:
%        [DA1 DA2], [DB1 DB2]
%        [DA1 DA2], [DB1]
%        [DA1],     [DB1 DB2]
%        [DA1],     [DB1]
%        [DA1 0],   [0 DB1]
%        [0 DA1],   [DB1 0]
%
%    sizeA/B     Equal to sizeA0/B0 if RESHAPE is not needed in MAIN 1
%    shiftC, delC    Variables controlling MAIN 3.
%    sizeisnew   1x2 logical array; activates reshaping of A and B.
%    idA/B       May change only if squashOK ~= 0
%    squashOK    If only A or B is a multi-block array (M-B) and the other
%                is single-block (1-B), it will be rearranged from N-D to
%                2-D. If both A and B are 1-B or M-B arrays, squashOK = 0.
%                If only A (or B) is a M-B array, squashOK = 1 (or 2).
%    sxtimesOK, timesOK, mtimesOK    Flags controlling MAIN 2 (TRUE/FALSE).
%    sumOK       Dimension along which SUM is performed. If SUM is not
%                needed, sumOK = 0.
% Initializing output arguments
    idA = idA0;
    idB = idB0;
     squashOK = 0;
    sxtimesOK = false;
      timesOK = false;
     mtimesOK = false;
        sumOK = 0;
    shiftC = 0;
    delC = 0;
% Checking for gross input errors
    NidA = numel(idA);
    NidB = numel(idB);
    idA1 = idA(1);
    idB1 = idB(1);
    if  NidA>2 || NidB>2 || NidA==0 || NidB==0 || ...
           ~isreal(idA1) ||    ~isreal(idB1)   || ...
        ~isnumeric(idA1) || ~isnumeric(idB1)   || ...
                 0>idA1  ||          0>idB1    || ... % negative 
         idA1~=fix(idA1) ||  idB1~=fix(idB1)   || ... % non-integer
         ~isfinite(idA1) ||  ~isfinite(idB1) % Inf or NaN               
        error('MULTIPROD:InvalidDimensionArgument', ...
        ['Internal-dimension arguments (e.g., [IDA1 IDA2]) must\n', ...
         'contain only one or two non-negative finite integers']);
    end
% Checking Syntaxes containing zeros (4b/c)
    declared_outer = false;
    idA2 = idA(NidA); % It may be IDA1 = IDA2 (1-D block)
    idB2 = idB(NidB);
    if any(idA==0) || any(idB==0)
        
        % "Inner products": C = MULTIPROD(A, B, [0 DA1], [DB1 0])
        if idA1==0 && idA2>0 && idB1>0 && idB2==0
            idA1 = idA2;
            idB2 = idB1;
        % "Outer products": C = MULTIPROD(A, B, [DA1 0], [0 DB1]) 
        elseif idA1>0 && idA2==0 && idB1==0 && idB2>0
            declared_outer = true;
            idA2 = idA1;
            idB1 = idB2;
        else
            error('MULTIPROD:InvalidDimensionArgument', ...
            ['Misused zeros in the internal-dimension arguments\n', ...
            '(see help heads 4b and 4c)']);
        end
        NidA = 1; 
        NidB = 1;
        idA = idA1;
        idB = idB1;
    elseif (NidA==2 && idA2~=idA1+1) || ...  % Non-adjacent IDs
           (NidB==2 && idB2~=idB1+1)
        error('MULTIPROD:InvalidDimensionArgument', ...
        ['If an array contains 2-D blocks, its two internal dimensions', ... 
        'must be adjacent (e.g. IDA2 == IDA1+1)']);
    end
% ESC - Case for which no reshaping is needed (both A and B are scalars)
    scalarA = isequal(sizeA0, [1 1]);
    scalarB = isequal(sizeB0, [1 1]);
    if scalarA && scalarB
        sizeA = sizeA0;
        sizeB = sizeB0;
        sizeisnew = [false false];
        timesOK = true; return
    end
% Computing and checking adjusted sizes
% The lengths of ADJSIZEA and ADJSIZEB must be >= IDA(END) and IDB(END)
    NsA = idA2 - length(sizeA0); % Number of added trailing singletons
    NsB = idB2 - length(sizeB0);
    adjsizeA = [sizeA0 ones(1,NsA)];
    adjsizeB = [sizeB0 ones(1,NsB)];
    extsizeA = adjsizeA([1:idA1-1, idA2+1:end]); % Size of EDs
    extsizeB = adjsizeB([1:idB1-1, idB2+1:end]);
    p = adjsizeA(idA1);
    q = adjsizeA(idA2);
    r = adjsizeB(idB1);
    s = adjsizeB(idB2);    
    scalarsinA = (p==1 && q==1);
    scalarsinB = (r==1 && s==1);
    singleA = all(extsizeA==1);
    singleB = all(extsizeB==1);
    if q~=r && ~scalarsinA && ~scalarsinB && ~declared_outer
       error('MULTIPROD:InnerDimensionsMismatch', ...
             'Inner matrix dimensions must agree.');
    end
% STEP 1/3 - DIMENSION SHIFTING (FIRST STEP OF AX)
%   Pipeline 1 (using TIMES) never needs left, and may need right shifting.
%   Pipeline 2 (using MTIMES) may need left shifting of A and right of B.
    shiftA = 0;
    shiftB = 0;
    diffBA = idB1 - idA1;    
    if scalarA % Do nothing
    elseif singleA && ~scalarsinB, shiftA = -idA1 + 1; %  Left shifting A
    elseif idB1 > idA1,            shiftA = diffBA;    % Right shifting A        
    end    
    if scalarB % Do nothing
    elseif singleB && ~scalarsinA, shiftB = -idB1 + 1; %  Left shifting B
    elseif idA1 > idB1,            shiftB = -diffBA;   % Right shifting B
    end
% STEP 2/3 - SELECTION OF PROPER ENGINE AND BLOCK SIZE ADJUSTMENTS
    addA  = 0; addB  = 0;
    delA  = 0; delB  = 0;
    swapA = 0; swapB = 0;
    idC1 = max(idA1, idB1);
    idC2 = idC1 + 1;
    checktimes = false;
    if (singleA||singleB) &&~scalarsinA &&~scalarsinB % Engine using MTIMES
        if singleA && singleB 
            mtimesOK = true;
            shiftC=idC1-1; % Right shifting C
            idC1=1; idC2=2;
        elseif singleA
            squashOK = 2;
            idB = [idB1, idB1+1] + shiftB;
        else % singleB
            squashOK = 1;
            idA = [idA1, idA1+1] + shiftA;
        end
        if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS
            % OK 
        elseif NidA==2        % 2) 2-D BLOCKS BY 1-D BLOCKS
            addB=idB1+1; delC=idC2;
        elseif NidB==2        % 3) 1-D BLOCKS BY 2-D BLOCKS
            addA=idA1; delC=idC1;
        else                  % 4) 1-D BLOCKS BY 1-D BLOCKS
            if declared_outer
                addA=idA1+1; addB=idB1;
            else
                addA=idA1; addB=idB1+1; delC=idC2;
            end
        end    
    else % Engine using TIMES (also used if SCALARA || SCALARB)
        
        sxtimesOK = true;
        if NidA==2 && NidB==2 % 1) 2-D BLOCKS BY 2-D BLOCKS
            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                checktimes=true;
            elseif scalarsinA || scalarsinB || ... % scal-by-mat
                (q==1 && r==1)  % vec-by-vec ("outer")
            elseif p==1 && s==1 % vec-by-vec ("inner")
                swapA=idA1; sumOK=idC1; checktimes=true;
            elseif s==1 % mat-by-vec
                swapB=idB1; sumOK=idC2;
            elseif p==1 % vec-by-mat
                swapA=idA1; sumOK=idC1;
            else % mat-by-mat
                addA=idA2+1; addB=idB1; sumOK=idC2; delC=idC2;
            end
        elseif NidA==2 % 2) 2-D BLOCKS BY 1-D BLOCKS
            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                addB=idB1; checktimes=true;
            elseif scalarsinA % scal-by-vec
                delA=idA1;
            elseif scalarsinB % mat-by-scal
                addB=idB1;
            elseif p==1 % vec-by-vec ("inner")
                delA=idA1; sumOK=idC1; checktimes=true;
            else % mat-by-vec
                addB=idB1; sumOK=idC2; delC=idC2;
            end
        elseif NidB==2 % 3) 1-D BLOCKS BY 2-D BLOCKS
            if scalarA || scalarB
                timesOK=true;                
            elseif scalarsinA && scalarsinB % scal-by-scal
                addA=idA1+1; checktimes=true;
            elseif scalarsinB % vec-by-scal
                delB=idB2;
            elseif scalarsinA % scal-by-mat
                addA=idA1+1;
            elseif s==1 % vec-by-vec ("inner")
                delB=idB2; sumOK=idC1; checktimes=true;
            else % vec-by-mat
                addA=idA1+1; sumOK=idC1; delC=idC1;
            end
        else % 4) 1-D BLOCKS BY 1-D BLOCKS
            if scalarA || scalarB
                timesOK=true;                
            elseif declared_outer % vec-by-vec ("outer")
                addA=idA1+1; addB=idB1;
            elseif scalarsinA && scalarsinB % scal-by-scal
                checktimes=true;
            elseif scalarsinA || scalarsinB % vec-by-scal
            else % vec-by-vec
                sumOK=idC1; checktimes=true;
            end
        end
    end
% STEP 3/3 - Adjusting the size of A and B. The size of C is adjusted
%            later, because it is not known yet.
    [sizeA, sizeisnew(1)] = adjustsize(sizeA0, shiftA, addA, delA, swapA);
    [sizeB, sizeisnew(2)] = adjustsize(sizeB0, shiftB, addB, delB, swapB);
    if checktimes % Faster than calling BBXFUN
        diff = length(sizeB) - length(sizeA);
        if isequal([sizeA ones(1,diff)], [sizeB ones(1,-diff)])
            timesOK = true;
        end
    end
end
function [sizeA, sizeisnew] = adjustsize(sizeA0, shiftA, addA, delA, swapA)
% ADJUSTSIZE  Adjusting size of a block array.
    % Dimension shifting (by adding or deleting trailing singleton dim.)
    if     shiftA>0, [sizeA,newA1] = addsing(sizeA0, 1, shiftA);
    elseif shiftA<0, [sizeA,newA1] = delsing(sizeA0, 1,-shiftA); 
    else , sizeA = sizeA0;  newA1  = false;
    end
    % Modifying block size (by adding, deleting, or moving singleton dim.)
    if      addA, [sizeA,newA2] = addsing(sizeA, addA+shiftA, 1); % 1D-->2D 
    elseif  delA, [sizeA,newA2] = delsing(sizeA, delA+shiftA, 1); % 2D-->1D
    elseif swapA, [sizeA,newA2] = swapdim(sizeA,swapA+shiftA); % ID Swapping
    else ,               newA2  = false;
    end
    sizeisnew = newA1 || newA2;
end
function [newsize, flag] = addsing(size0, dim, ns)
%ADDSING   Adding NS singleton dimensions to the size of an array.
%   Warning: NS is assumed to be a positive integer.
%   Example: If the size of A is ..... SIZE0 = [5 9 3]
%            NEWSIZE = ADDSING(SIZE0, 3, 2) is [5 9 1 1 3]
    if dim > length(size0)
        newsize = size0;
        flag = false;
    else 
        newsize = [size0(1:dim-1), ones(1,ns), size0(dim:end)];
        flag = true;
    end
end  
function [newsize, flag] = delsing(size0, dim, ns)
%DELSING   Removing NS singleton dimensions from the size of an array.
%   Warning: Trailing singletons are not removed
%   Example: If the size of A is SIZE0 = [1 1 1 5 9 3]
%            NEWSIZE = DELSING(SIZE, 1, 3) is  [5 9 3]
    if dim > length(size0)-ns % Trailing singletons are not removed
        newsize = size0;
        flag = false;
    else % Trailing singl. added, so NEWSIZE is guaranteed to be 2D or more
        newsize = size0([1:dim-1, dim+ns:end, dim]);
        flag = true;
    end
end
function [newsize, flag] = swapdim(size0, dim)
%SWAPDIM   Swapping two adjacent dimensions of an array (DIM and DIM+1).
%   Used only when both A and B are multi-block arrays with 2-D blocks.
%   Example: If the size of A is .......... 5�(6�3)
%            NEWSIZE = SWAPIDS(SIZE0, 2) is 5�(3�6)
    newsize = [size0 1]; % Guarantees that dimension DIM+1 exists.
    newsize = newsize([1:dim-1, dim+1, dim, dim+2:end]);
    flag = true;
end
function b = multitransp(a, dim)
%MULTITRANSP  Transposing arrays of matrices.
%    B = MULTITRANSP(A) is equivalent to B = MULTITRANSP(A, DIM), where
%    DIM = 1.
%
%    B = MULTITRANSP(A, DIM) is equivalent to
%    B = PERMUTE(A, [1:DIM-1, DIM+1, DIM, DIM+2:NDIMS(A)]), where A is an
%    array containing N P-by-Q matrices along its dimensions DIM and DIM+1,
%    and B is an array containing the Q-by-P transpose (.') of those N
%    matrices along the same dimensions. N = NUMEL(A) / (P*Q), i.e. N is
%    equal to the number of elements in A divided by the number of elements
%    in each matrix.
%
%    MULTITRANSP, PERMUTE and IPERMUTE are a generalization of TRANSPOSE
%    (.') for N-D arrays.
%
%    Example:
%       A 5-by-9-by-3-by-2 array may be considered to be a block array
%       containing ten 9-by-3 matrices along dimensions 2 and 3. In this
%       case, its size is so indicated:  5-by-(9-by-3)-by-2 or 5x(9x3)x2.
%       If A is ................ a 5x(9x3)x2 array of 9x3 matrices,
%       C = MULTITRANSP(A, 2) is a 5x(3x9)x2 array of 3x9 matrices.
%
%    See also PERMUTE, IPERMUTE, MULTIPROD.
% $ Version: 1.0 $
% CODE      by:                 Paolo de Leva (IUSM, Rome, IT) 2005 Sep 9
% COMMENTS  by:                 Code author                    2006 Nov 21
% OUTPUT    tested by:          Code author                    2005 Sep 13
% -------------------------------------------------------------------------
% Setting DIM if not supplied.
if nargin == 1, dim = 1; end
% Transposing
order = [1:dim-1, dim+1, dim, dim+2:ndims(a)];
b = permute(a, order);
end