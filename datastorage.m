classdef datastorage
  %static
  properties(Constant)
    %default value of some internal parameters
    parameter_list={...
      'start',          time.zero_date,@(i) isdatetime(i);...
      'stop',           time.inf_date, @(i) isdatetime(i);...
      'debug',          false,         @(i) islogical(i);...
      'start_timestamp_only',true,     @(i) islogical(i);...
      'inclusive',      false,         @(i) islogical(i);...
    };
  end
  %public
  properties(GetAccess=public,SetAccess=public)
    debug
    inclusive
  end
  %read only
  properties(SetAccess=private)
    data
  end
  %private (visible only to this object)
  properties(GetAccess=private)
    starti
    stopi
  end
  %calculated only when asked for
  properties(Dependent)
    start
    stop
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(datastorage.parameter_list); end
      out=v.picker(varargin{:});
    end

%     function out=id(in,data_type,varargin)
%       switch lower(data_type)
%       case {'dataname','datanames','dn'}
%         assert(numel(varargin)==1,['If datastorage.id is called with ''data_type'' as ''',data_type,...
%           ''', then can only handle one optional input argument, not ',num2str(numel(varargin)),'.'])
%         out=datanames(in,varargin{:});
%       case {'product','products','dp'}
%         out=dataproduct(in,varargin{:});
%       otherwise
%         assert(numel(varargin)==0,['If datastorage.id is called with ''data_type'' as ''',data_type,...
%           ''', then can not handle any optional input argument (while ',num2str(numel(varargin)),' were given).'])
%         out=datanames(in);
%         out=out.(data_type);
%       end
%     end
  end
  methods
    %% constructor
    function obj=datastorage(varargin)
      obj.log('@','in','varargin',varargin)
      % reset data type list
      obj=obj.data_init;
      %create argument object, declare and parse parameters, save them to obj
      [~,~,obj]=varargs.wrap('sinks',{obj},'sources',{datastorage.parameters('obj')},varargin{:});
      % data is added to this initialized object with the 'init' method (see below)
      obj.log('@','out')
    end
    function out=isempty(obj)
      out=isempty(obj.data);
    end
    %% debug logging
    function log(obj,varargin)
      %don't do anything unless debug mode is on
      if ~obj.debug
        return
      end
      %sanity on the nr of input arguments
      assert(mod(numel(varargin),2)==0,['the number of input arguments must be even, not ',num2str(numel(varargin)),'.'])
      %make room for message
      msg=cell((nargin-1)/2,1);
      %loop over all arguments
      for i=1:numel(msg)
        idx_name=2*i-1;
        idx_value=2*i;
        %don't show empty arguments
        if isempty(varargin{idx_value})
          continue;
        end
        msg{i}=[varargin{idx_name},':',str.show(varargin{idx_value})];
      end
      %clean up empty cells
      msg=cells.rm_empty(msg);
      %show the debug message
      str.say('stack_delta',1,strjoin(msg,', '))
    end
    %% dataname operations (done at the root level of the obj.data structure)
    function obj=data_init(obj)
      obj.data=struct([]);
    end
    function out=isdata_empty(obj,dn)
      %reduce dataname to common object
      dn=datanames(dn);
      %check is empty
      out=~isfield(obj.data,dn.name_clean) || isempty(obj.data.(dn.name_clean));
      %in case field_path is non-empty, check that too
      if ~out && ~isempty(dn.field_path)
        out=isempty(structs.get_value(obj.data.(dn.name_clean),dn.field_path));
      end
    end
    function out=dataname_list(obj,non_empty)
      out=datanames.array(fieldnames(obj.data));
      if exist('non_empty','var') && non_empty
        empty_idx=cellfun(@(i) isempty(obj.data.(i.name)),out);
        out=out(~empty_idx);
      end
    end
    function peek(obj,dn_list,varargin)
      if obj.isempty
        return
      end
      %look at all datanames by default
      if ~exist('dn','var') || isempty(dn_list)
        dn_list='all';
      end
      %parse optional arguments
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('info_methods',{'size','nr_gaps','start','stop'}, @(i) ischar(i) || iscellstr(i));
      p.addParameter('tab',[32,20], @(i) isnumeric(i));
      p.parse(varargin{:});
      %show header
      disp(str.tablify(p.Results.tab,'product',p.Results.info_methods))
      %retrieve global field path list
      obj_list=obj.data_get(dn_list);
      dn_list=obj.data_list(dn_list);
      %loop over all retrieved objects
      for i=1:numel(obj_list)
        msg=cell(size(p.Results.info_methods));
        if isempty(obj_list{i})
          msg(:)={'-'};
        else
          for m=1:numel(p.Results.info_methods)
            %any(isprop(...)) catches when obj_list{i} is a string 
            if ismethod(obj_list{i},p.Results.info_methods{m}) || any(isprop(obj_list{i},p.Results.info_methods{m}))
              msg{m}=obj_list{i}.(p.Results.info_methods{m});
            else
              msg{m}='N/A';
            end
          end
        end
        disp(str.tablify(p.Results.tab,dn_list{i}.str,msg))
      end
    end
    function da(obj,dn_list,varargin)
      if obj.isempty
        return
      end
      %look at all datanames by default
      if ~exist('dn','var') || isempty(dn_list)
        dn_list='all';
      end
      %parse optional arguments
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('info_method','nr_gaps', @(i) ischar(i));
      p.addParameter('tab',[32,11], @(i) isnumeric(i));
      p.addParameter('period',days(1), @(i) isduration(i));
      p.addParameter('group',10, @(i) isscalar(i) && isnumeric(i));
      p.parse(varargin{:});
      %gather list of periods
      [startlist,stoplist]=time.list(obj.start,obj.stop,p.Results.period);
      %split list of periods into groups of 10 (otherwise it' difficult to see anything in the screen)
      n=ceil(numel(startlist)/p.Results.group);
      list0=cell(1,n);
      list1=cell(1,n);
      for i=1:n
        start_idx=(i-1)*p.Results.group+1;
        stop_idx=min([i*p.Results.group,numel(startlist)]);
        list0{i}=startlist(start_idx:stop_idx);
        list1{i}=stoplist( start_idx:stop_idx);
      end
      %show header
      disp(['Showing ',p.Results.info_method,' for periods of ',str.show(p.Results.period),'.'])
      for l=1:n
        disp(str.tablify(p.Results.tab,'product',list0{l}))
        %retrieve global field path list
        obj_list=obj.data_get(dn_list);
        dn_list=obj.data_list(dn_list);
        %loop over all retrieved objects
        for i=1:numel(obj_list)
          msg=cell(size(list0{l}));
          if isempty(obj_list{i})
            msg(:)={'-'};
          else
            for m=1:numel(list0{l})
              %try to get the requested information
              try
                %get the requested info for the current data period
                 msg{m}=obj_list{i}.trim(list0{l}(m),list1{l}(m)).(p.Results.info_method);
              catch
                %patch with N/A
                msg{m}='N/A';
              end
            end
          end
          disp(str.tablify(p.Results.tab,dn_list{i}.str,msg))
        end
      end
    end
    %% value operations (basically wraps some methods in the structs object)
    function out=value_get(obj,dn)
      %reduce dataname to common object
      dn=datanames(dn);
      %retrieve value ('true' allows for structs.get_value to loop along dn.field_path
      %                until it refers to field available in obj.data.(dn.name_clean))
      out=structs.get_value(obj.data.(dn.name_clean),dn.field_path,true);
    end
    function obj=value_set(obj,dn,value)
      %reduce dataname to common object
      dn=datanames(dn);
      %make sure this dataname has been initiated
      if ~isfield(obj.data,dn.name_clean)
        obj.data(1).(dn.name_clean)={};
      end
      obj.data.(dn.name_clean)=structs.set_value(...
        obj.data.(dn.name_clean),...
        dn.field_path,...
        value...
      );
    end
    %% data operations
    function out=isdata_leaf(obj,dn,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('only_non_empty',false,@(i) islogical(i) && isscalar(i));
      p.parse(varargin{:});
      %reduce dataname to common object
      dn=datanames(dn);
      %determine leafness (structs.isleaf returns empty structures as leafs by default, only_non_empty changes that)
      out=~isfield(obj.data,dn.name_clean) || structs.isleaf(obj.data.(dn.name_clean),dn.field_path,p.Results.only_non_empty);
    end
    function dn_list=data_list(obj,dn,varargin)
      if numel(dn)==0
        %trivial call
        dn_list={};
        return
      elseif ischar(dn) && strcmp(dn,'all')
        %handle named inputs
        dn_list=obj.data_list(obj.dataname_list(true));
        return
      elseif numel(dn)>1 && ~ischar(dn)
        %vector mode
        assert(iscell(dn),['when dealing with non-scalar input ''dn'', can only deal with cells, not with ',class(dn),'.'])
        dn_list=cells.flatten(cellfun(@obj.data_list,dn,'UniformOutput',false));
        return
      elseif isa(dn,'datanames') && any(strcmp(dn.field_path,'*'))
        %get field_path list that matches input globbed field_path
        field_path_globbed=structs.field_list_glob(obj.data.(dn.name_clean),dn.field_path);
        %build list of datanames for this globbed field_path
        dn_list=cell(size(field_path_globbed));
        for i=1:numel(field_path_globbed)
          %patch globbed field_path
          dn_list{i}=datanames(dn.name,field_path_globbed{i});
        end
        return
      else
        %reduce dataname to common object
        dn=datanames(dn);
      end
      %scalar mode
      %check if this dataname is a leaf product
      if obj.isdata_leaf(dn,varargin{:})
        field_path_list={dn.field_path};
      else
        %get all field names under dn.field_path
        field_path_list=structs.field_list(structs.get_value(obj.data.(dn.name_clean),dn.field_path));
        %prepend current dn.field_path
        for i=1:numel(field_path_list)
          field_path_list{i}=[dn.field_path,field_path_list{i}];
        end
      end
      %handle empty field_path_list
      if isempty(field_path_list)
        dn_list={};
        return
      end
      %be sure this is a cell of cells
      assert(cells.iscellofcells(field_path_list,1),'variable ''field_path_list'' must be a cell of cells.')
      %propagate field path list to list of dataname objects
      dn_list=cellfun(@(i) dn.set_field_path(i),field_path_list,'UniformOutput',false);
    end
    function out=data_get(obj,dn)
      %reduce dataname to common object
      dn_list=obj.data_list(dn);
      %retrieve
      out=cell(size(dn_list));
      for i=1:numel(dn_list)
        if isfield(obj.data,dn_list{i}.name_clean)
          out{i}=obj.value_get(dn_list{i});
        end
      end
    end
    %scalar mode is useful when calling this function as part of an expression
    function out=data_get_scalar(obj,dn)
      %call mother routine
      out=data_get(obj,dn);
      %some sanity
      assert(numel(out)==1,'This routine cannot handle non-scalar data entries.')
      %reduce to scalar
      out=out{1};
    end
    function obj=data_set(obj,dn,values)
      %reduce dataname to common object
      dn_list=obj.data_list(dn);
      values=cells.scalar(values,'set');
      
      %TODO: I question the usefulness of this, need to be sure it does not break things
%       %handle scalar values, unwrap them to all field paths
%       if numel(values)==1 && ~iscell(values)
%         tmp=values;
%         values=cell(size(dn_list));
%         values{:}=tmp;
%       end
      
      %check  if fiels paths and values have the same length
      assert(numel(dn_list)==numel(values),'Cannot handle inputs ''in'' and ''values''.')
      %propagate 
      for i=1:numel(dn_list)
        %try to match field names in structures with field paths
        if isstruct(values{i})
          tmp=structs.get_value(values{i},dn_list{1}.field_path);
          if ~isempty(tmp)
            obj=obj.value_set(dn_list{i},tmp);
          else
            assert(isempty(dn_list{1}.field_path),[...
              'Cannot find field path ''',dn_list{1}.field_path_str,''' in product ''',dn_list{1}.filename,'''.'])
            obj=obj.value_set(dn_list{i},values{i});
          end
        else
          obj=obj.value_set(dn_list{i},values{i});
        end
      end
    end
    %% vector data operations, done over multiple datanames and all corresponding 'field_path's
    % Applies the function f to all data entries and returns the result,
    % so it transforms (tr) a list of objects into something else.
    function out=vector_func_tr(obj,dn,f,varargin)
      %collapse requested data into cell array and operate on the cell array
      out=cellfun(@(i) f(i,varargin{:}),obj.data_get(dn));
    end
    % Applies the function f to all data entries under 'dn' and
    % saves the resulting objects back to the same set, using data_set(...)
    function obj=vector_func_op(obj,dn,f,varargin)
      %operate on the requested set
      values=obj.vector_tr(dn,f,varargin{:});
      %sanity
      if ~iscell(values)
        error([mfilename,': function ',fun2str(f),' must return a cell array, not a ',class(values),'.'])
      end
      %propagate cell array back to object
      obj=obj.data_set(dn,values);
    end
    % Retrieves the values of particular field or a method from the
    % object stored at the leaves of all data entries (used to be called vector_sts)
    function [out,dn_list]=vector_method_tr(obj,dn,method,varargin)
      %get data entries under 'dn'
      dn_list=obj.data_list(dn);
      %retrieve list of objects, i.e. the values given by obj.data_get
      obj_list=obj.data_get(dn_list);
      %make room for outputs
      out=cell(size(obj_list));
      %loop over all objects
      for i=1:numel(obj_list)
        if isempty(obj_list{i})
          out{i}=[];
        else
          %check if this object responds to this
          if structs.respondto(obj_list{i},method)
            if numel(varargin)==0
              out{i}=obj_list{i}.(method);
            else
              out{i}=obj_list{i}.(method)(varargin{:});
            end
          end
        end
      end
    end
    % Applies the given method to all data entries under 'dn' and
    % saves the resulting objects back to the same set, using data_set(...)
    function obj=vector_method_op(obj,dn,method,varargin)
      %operate on the requested set
      [values,dn_list]=vector_method_tr(obj,dn,method,varargin{:});
      %propagate cell array back to object
      obj=obj.data_set(dn_list,values);
    end
    % Sets the values of particular field or zero-input method of the
    % object stored at the leaves of all data entries
    function obj=vector_method_set(obj,dn,method,value)
      %get data entries under 'dn'
      dn_list=obj.data_list(dn);
      %retrieve list of objects, i.e. the values given by obj.data_get
      obj_list=obj.data_get(dn_list);
      %loop over all objects
      for i=1:numel(obj_list)
        if ~isempty(obj_list{i})
          %check if this object responds to this
          if structs.respondto(obj_list{i},method)
            obj_list{i}.(method)=value;
          end
        end
      end
      %propagate cell array back to object
      obj=obj.data_set(dn_list,obj_list);
    end
    % Applies a 2-argument method to the data given by obj and obj2
    function obj=vector_obj_op2(obj,dn,obj2,method,varargin)
      %operate on the requested set
      values1=obj.data_get(dn);
      values2=obj2.data_get(dn);
      %sanity
      if numel(values1) ~= numel(values2)
        error([mfilename,': the given dataname list does not correspond to the same number of data entries in both input obj.'])
      end
      if any(numel(cells.rm_empty(values1))~=numel(cells.rm_empty(values2)))
        error([mfilename,': the given dataname list does not have data in both objects. Debug needed.'])
      end
      %outputs
      result=cell(size(values1));
      %operate
      for i=1:numel(values1)
        result{i}=values1{i}.(method)(values2{i},varargin{:});
      end
      %propagate cell array with results back to object
      obj=obj.data_set(dn,result);
    end
    % Applies a 2-argument method to the data given by dn1 and dn2 and
    % saves it in dn_out
    function obj=vector_obj_op3(obj,dn1,dn2,method,dn_out,varargin)
      %operate on the requested set
      values1=obj.data_get(dn1);
      values2=obj.data_get(dn2);
      %sanity
      if numel(values1) ~= numel(values2)
        error([mfilename,': the given dataname list does not correspond to the same number of data entries in both input obj.'])
      end
      %outputs
      result=cell(size(values1));
      %operate
      for i=1:numel(values1)
        %enforce common time domain both ways
        values2{i}=values2{i}.interp2(values1{i},varargin{:},'check_width',false,'skip_par_check',{'lmax','labels','y_units'});
        values1{i}=values1{i}.interp2(values2{i},varargin{:},'check_width',false,'skip_par_check',{'lmax','labels','y_units'});
        %apply method
        result{i}=values1{i}.(method)(values2{i},varargin{:});
      end
      %propagate cell array with results back to object
      obj=obj.data_set(dn_out,result);
    end
    %% start/stop operations
    function io=start_criteria(obj,io)
      assert(~isempty(io),'Cannot handle empty ''io''.')
      assert(iscell(io) && all(cellfun(@isdatetime,io)),'Can only handle cell arrays of datetime')
      %need to filter out invalids (zero and inf dates)
      invalid_idx=~time.isvalid(io);
      if any(invalid_idx); io(invalid_idx)={time.zero_date}; end
      %hanle inclusive datasets (span over all boundaries)
      if obj.inclusive
        if isempty(io)
          io=time.zero_date;
        else
          io=min([io{:}]);
        end
      else
        io=max([io{:}]);
      end
    end
    function io=stop_criteria(obj,io)
      assert(~isempty(io),'Cannot handle empty ''values''.')
      assert(iscell(io) && all(cellfun(@isdatetime,io)),'Can only handle cell arrays of datetime')
      %need to filter out invalids (zero and inf dates)
      invalid_idx=~time.isvalid(io);
      if any(invalid_idx); io(invalid_idx)={time.inf_date}; end
      %hanle inclusive datasets (span over all boundaries)
      if obj.inclusive
        if isempty(io)
          io=time.inf_date;
        else
          io=max([io{:}]);
        end
      else
        io=min([io{:}]);
      end
    end
    function out=start_retrieve(obj,dn)
      out=obj.start_criteria(obj.vector_method_tr(dn,'start'));
    end
    function out=stop_retrieve(obj,dn)
      out=obj.stop_criteria(obj.vector_method_tr(dn,'stop'));
    end
    function out=get.start(obj)
      out=obj.starti;
    end
    function out=get.stop(obj)
      out=obj.stopi;
    end
    function obj=set.start(obj,start)
      %save current global start value
      old_start=obj.starti;
      %check if updating is needed
      if isempty(obj.starti) || start==obj.start_criteria({obj.starti,start})
        obj.log('@','start update','from',obj.starti,'to',start)
        %update internal record
        obj.starti=start;
      end
      %check if something changed (also trims newly added data that starts before old_start)
      if isempty(old_start) || start~=old_start
        %trim all data entries
        obj=obj.vector_method_set('all','start',obj.starti);
      end
    end
    function obj=set.stop(obj,stop)
      %save current global stop value
      old_stop=obj.stopi;
      %update if needed
      if isempty(obj.stopi) || stop==obj.stop_criteria({obj.stopi,stop})
        obj.log('@','stop update','from',obj.stopi,'to',stop)
        %update internal record
        obj.stopi=stop;
      end
      %check if something changed (also trims newly added data that end after old_stop)
      if isempty(old_stop) || stop~=old_stop
        %trim all data entries
        obj=obj.vector_method_set('all','stop',obj.stopi);
      end
    end
    function obj_out=trim(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('start',  obj.start,@(i) isdatetime(i) && isscalar(i));
      p.addParameter('stop',   obj.stop, @(i) isdatetime(i) && isscalar(i));
      p.addParameter('dn_list',obj.data_list('all'), @(i) iscell(i));
      p.parse(varargin{:});
      obj_out=obj.data_init.data_set(p.Results.dn_list,obj.data_get(p.Results.dn_list));
      obj_out.start=p.Results.start;
      obj_out.stop=p.Results.stop;
    end
    function obj=startstop_update(obj,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(varargin{:});
      % check if start/stop were part of varargin
      for i={'start','stop'}
        if ~any(strcmp(p.UsingDefaults,i{1}))
          obj.(i{1})=p.Results.(i{1});
        end
      end
    end
    function obj=startstop_retrieve_update(obj,product,inclusive_now)
      %save current value of object-wide inclusive flag
      inclusive_old=obj.inclusive;
      %handle inclusive flag
      if exist('inclusive_now','var') && ~isempty(inclusive_now)
        obj.inclusive=inclusive_now;
      elseif product.ismdfield('inclusive')
        obj.inclusive=product.mdget('inclusive');       
      end
      %retrieve start/stop
      start_now=obj.start_retrieve(product);
      stop_now=obj.stop_retrieve(product);
      %update start/stop
      obj=obj.startstop_update('start',start_now,'stop',stop_now);
      %inform
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
      %recover old value of inclusive flag
      obj.inclusive=inclusive_old;
    end
    %% time operations
    function out=isteq(obj,dn)
      %retrieve global field path list
      global_field_path_list=obj.data_list(dn);
      %default value
      out=true;
      %use the first object as reference
      obj1=obj.data_get_scalar(global_field_path_list{1});
      %loop over all remaining objects
      for i=2:numel(global_field_path_list)
        %check of the time domains agree
        if ~obj1.isteq(obj.data_get_scalar(global_field_path_list{i}))
          out=false;
          return
        end
      end
    end
    function obj=interp(obj,dn)
      %trivial call
      if obj.isteq(dn)
        return
      end
      %retrieve global field path list
      global_field_path_list=obj.data_list(dn);
      %use the first object as reference
      obj_ref=obj.data_get_scalar(global_field_path_list{1});
      %loop over all remaining objects (do it twice to be sure all time domains are the same)
      for j=0:1
        %the last object is for sure already merged on the second run, skip that
        for i=2:numel(global_field_path_list)-j
          %interpolate (with merging the time domains)
          obj=obj.data_set(global_field_path_list{i},...
            obj.data_get_scalar(global_field_path_list{i}).interp2(obj_ref)...
          );
        end
      end
      %sanity
      assert(obj.isteq(dn),'Could not merge all time domains. Debug needed.')
    end
    %% length operations
    function out=length(obj,dn)
      if ~exist('dn','var') || isempty(dn); dn='all'; end
      out=cell2mat(obj.vector_method_tr(dn,'length'));
    end
    %% debug utils
    function print(obj,dn,varargin)
      %show everything by default
      if ~exist('dn','var') || isempty(dn)
        dn='all';
      end
      %retrieve global field path list
      obj_list=obj.data_get(dn);
      dn_list=obj.data_list(dn);
      %loop over all retrieved objects
      for i=1:numel(obj_list)
        disp(dn_list{i}.str)
        if isempty(obj_list{i})
          disp('empty')
        else
          obj_list{i}.print(varargin{:});
        end
      end
    end
    %% product interface
    function [product,dn]=product_get(~,id,varargin)
      product=dataproduct(id,varargin{:});
      if nargout>1; dn=product.dataname; end
    end
    %% load/save operations
    function out=data_edges(obj,product,mode)
      switch mode
      case 'start'; f=@max;
      case 'stop';  f=@min;
      otherwise
        error(['Unknown mode ''',mode,'''.'])
      end
      if product.ismdfield(mode)
        out=f([obj.(mode),product.(mode)]);
      else
        out=obj.(mode);
      end
    end
    function [obj,success]=load(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      if product.mdget('plot_product','default',false)
        product_type='plot';
      else
        product_type='data';
      end
      % get the file list
      file_list=product.file(product_type,...
        'start',obj.data_edges(product,'start'),... %this (and stop) only reduces the number of files to be loaded, it does not garantee that the data is trimmed correctly
        'stop',obj.data_edges(product,'stop'),...
        'start_timestamp_only',datastorage.parameters('start_timestamp_only'),... %this needs to be in agreement with what is used in datastorage.save
        'discover',false,... %TODO: it may make sense to set this to true, dunno why it's not
        'ensure_dir',false,...
        varargin{:}...
      );
      % get file existence
      file_exists=product.isfile(product_type,...
        'start',obj.data_edges(product,'start'),...
        'stop',obj.data_edges(product,'stop'),...
        'start_timestamp_only',datastorage.parameters('start_timestamp_only'),... %this needs to be in agreement with what is used in datastorage.save
      varargin{:});
      % check if any file is missing
      if any(~file_exists)
        % debug output
        success=false;
        if isempty(file_list)
          obj.log('@','out','no files to load for',product)
        else
          obj.log('@','out','file(s) missing for',[product.str,newline,strjoin(file_list(~file_exists),newline)])
        end
        return
      end
      obj.log('nr of files to load',numel(file_list),'from dir',fileparts(file_list{1}))
      %if this is a plot product, then skip loading data
      if ~product.mdget('plot_product','default',false)
        %loop over all files
        for f=1:numel(file_list)
          %load data in this file
          load(file_list{f},'s');
          if f==1
            %move to cumulative variable
            s_out=s;
          else
            %%append to cumulative data
            if isstruct(s_out)
              %do that for all fields in this structure
              s_out=structs.objmethod2('append',s_out,s);
            else
              %simple call
              s_out=s_out.append(s);
            end
          end
          % debug output
          start_now=obj.start_criteria(cells.scalar(cells.rm_empty(structs.get_value_all(structs.objmethod('start',  s_out))),'set'));
          stop_now =obj.stop_criteria( cells.scalar(cells.rm_empty(structs.get_value_all(structs.objmethod('stop',   s_out))),'set'));
          nr_gaps  =max(cells.c2m(cells.rm_empty(structs.get_value_all(structs.objmethod('nr_gaps',s_out)))));
          [~,file_now,e]=fileparts(file_list{f});
          obj.log(['loaded file ',num2str(f),' '],[' ',file_now,e],'cum start',start_now,'cum stop',stop_now,'gaps',nr_gaps)
        end
        % enforce start/stop in metadata
        %NOTICE: this handles the case when data is saved first and the start/stop metadata entries are increased/decreased (trimming the saved data)
        %NOTICE: this does not handle decreasing/increased the start/stop metadata such that no new data file(s) is created (e.g. when using yearly or global storage_period)
        s_out=structs.objmethod('trim',s_out,obj.data_edges(product,'start'),obj.data_edges(product,'stop'));
        %TEMPORARY: enforce proper gravity object labels
        gravitylabelfix=false;
        if isstruct(s_out)
          fl=structs.field_list(s_out);
          for i=1:numel(fl)
            fv=structs.get_value(s_out,fl{i});
            if isa(fv,'gravity')
              [fv,r]=fv.setlabels;
              if r; s_out=structs.set_value(s_out,fl{i},fv);gravitylabelfix=true; end
            end
          end
        else
          if isa(s_out,'gravity')
            [s_out,gravitylabelfix]=s_out.setlabels;
          end
        end
        % save the data in the object
        obj=obj.data_set(product,s_out);
        %TEMPORARY: enforce proper gravity object labels
        if gravitylabelfix; obj.save(product); end
        % update start/stop
        obj=obj.startstop_retrieve_update(product);
      end
      % debug output
      success=true;
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function save(obj,product,varargin)
      %ignore plot products, file saving is done inside the init method
      if product.mdget('plot_product','default',false)
        return
      end
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('force',false,@(i) islogical(i) && isscalar(i)); %forces saving the data even if the data is already saved 
      % parse it
      p.parse(varargin{:});
      % get the data to be saved
      % NOTICE: it is *always* everything under this product.name, never resolved at the level of the field_path
      % this is to be in agreement with how the filenames are built in datanames
      s_out=obj.value_get(product.name);
      if isstruct(s_out)
        % serialize the data
        s_cells=structs.get_value_all(s_out);
      else
        s_cells={s_out};
      end
      % catch unvalid s_out
      if any(~cells.allrespondto(s_cells,{'start','stop'}))
        obj.log('@','out','no data to save for',product)
        return
      end
      % get list of start/stop dates from the data
      startlist=cells.c2m(cellfun(@(i) i.start,s_cells,'UniformOutput',false));
      stoplist =cells.c2m(cellfun(@(i) i.stop, s_cells,'UniformOutput',false));
      % get the file list
      [file_list,startlist,stoplist]=product.file('data',...
        'start',min(startlist),...
        'stop',max(stoplist),...
        'start_timestamp_only',datastorage.parameters('start_timestamp_only'),... %this needs to be in agreement with what is used in datastorage.load
        'discover',false,...
        varargin{:}...
      );
      %check if nothing is to be save
      if isempty(file_list)
        obj.log('@','out','no files to save for',product)
        return
      end
      obj.log('nr of files to save',numel(file_list))
      %loop over all files
      for f=1:numel(file_list)
        % global storage periods not have start/stop lists
        if isempty(startlist) || isempty(stoplist)
          %paranoid sanity
          assert(numel(file_list)==1,['expecting there to be a single file name to save, not ',num2str(numel(file_list)),'.'])
          s=s_out;
        else
          %trim the data to the current start/stop epochs
          if isstruct(s_out)
            %do that for all fields in this structure
            s=structs.objmethod('trim',s_out,startlist(f),stoplist(f));
          else
            %simple call
            s=s_out.trim(startlist(f),stoplist(f));
          end
        end
        %skip if this entry is empty
        if structs.isempty(s)
          obj.log(['skip saving file ',num2str(f)],['no data for ',datestr(startlist(f),29)])
          continue
        end
        %check if this file already exists
        if exist(file_list{f},'file')
          %save the data if force is true
          if p.Results.force
            replace_flag=true;
          else
            %save new variable
            s_new=s;
            %load this file
            load(file_list{f},'s');
            %re-align data if needed
            if ~strcmp(class(s),class(s_new))
              disp(['Data type in newly-computed object (',class(s_new),...
              ') differs from what is saved in file ',file_list{f},' (',class(s),')',...
              '; replacing data saved in file'])
              %discard saved data
              replace_flag=true;
              s=s_new;
            elseif isstruct(s)
              replace_flag=false;
              [good_fnl,fnl,fnl_new]=structs.iseq_field_list(s,s_new);
              if good_fnl
                %skip saving if start/stop dates in s_new are not wider than those already in file
                for i=1:numel(fnl)
                  s_new_now=structs.get_value(s_new,fnl{i});
                  s_old    =structs.get_value(s,    fnl{i});
                  if isempty(s_new_now)
                    replace_flag=false;
                  elseif isempty(s_old)
                    replace_flag=true;
                  elseif s_new_now.start>=s_old.start && s_new_now.stop<=s_old.stop
                    %do nothing
                  elseif isempty(s_new_now) && isempty(s_old)
                    %also do nothing
                  elseif s_new_now.start<s_old.start || s_new_now.stop>s_old.stop
                    %augment saved data
                    s=structs.set_value(s,fnl{i},s_old.augment(s_new_now));
                    replace_flag=true;
                  else
                    replace_flag=true;
                  end
                end
              else
                disp(['structure fields in newly-computed object (',...
                  strjoin(cellfun(@(i) strjoin(i,'.'),fnl_new,'UniformOutput',false),', '),...
                  ') differ from what is saved in file ',file_list{f},' (',...
                  strjoin(cellfun(@(i) strjoin(i,'.'),fnl,'UniformOutput',false)),...
                  '); replacing data saved in file'])
                %discard saved data
                replace_flag=true;
                s=s_new;
              end
            else
              assert(...
                (ismethod(s,'start') || isprop(s,'start')) && ...
                (ismethod(s,'stop' ) || isprop(s,'stop' )),...
                ['Cannot handle data of type ',class(s),', no support for start/stop.']...
              );
              assert(...
                (ismethod(s_new,'start') || isprop(s_new,'start')) && ...
                (ismethod(s_new,'stop' ) || isprop(s_new,'stop' )),...
                ['Cannot handle data of type ',class(s_new),', no support for start/stop.']...
              );
              replace_flag=false;
              if s_new.start>=s.start && s_new.stop<=s.stop
                %do nothing
              elseif s_new.start<s.start || s_new.stop>s.stop
                %augment saved data
                s=s.augment(s_new);
                replace_flag=true;
              else
                replace_flag=true;
              end
            end
          end
          if replace_flag
            if p.Results.force
              obj.log(['replacing file ',num2str(f)],'force flag is true')
            else
              obj.log(['replacing file ',num2str(f)],'new data found')
            end
          else
            obj.log(['skip saving file ',num2str(f)],'all data already saved')
            s=[];
          end
        end
        %save this section of the data to the file (unless empty)
        if ~isempty(s) && numel(s)>0 && all(size(s)>0)
          % debug output
          obj.log(['saving file ',num2str(f)],file_list{f})
          save(file_list{f},'s');
        end
      end
      % debug output
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function export(obj,id,varargin)
      obj.log('@','in','id',id,'varargin',varargin,'start',obj.start,'stop',obj.stop)
      %retrieve product info
      product=obj.product_get(id,varargin{:});
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('filename',product.codename,@(i) ischar(i));
      % parse it
      p.parse(varargin{:});
      % get the data 
      s_out=obj.value_get(product);
      % get the data names
      s_names=structs.field_list(s_out);
      if isempty(s_names)
        s_out.export([product.codename,'.dat'],'ascii');
      else
        % loop over all data elements
        for i=1:numel(s_names)
          s_now=structs.get_value(s_out,s_names{i});
          %----------------------------------------------------------
          %this can be changed according to what needs to be exported
          s_now.lmax=4;
          %----------------------------------------------------------
          s_now.export([product.codename,'.',strjoin(s_names{i},'.'),'.dat'],'ascii');
        end
      end
      % debug output
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    %% datatype initialization
    function obj=init(obj,id,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('reload', false, @(i) islogical(i) && isscalar(i)); %re-loads the data if already loaded
      p.addParameter('force',  false, @(i) islogical(i) && isscalar(i)); %calls the init method even if the data is already saved
      p.addParameter('subcat',    '', @(i) ischar(i));
      % parse it
      p.parse(varargin{:});
      %update start/stop if given in inputs
      obj=obj.startstop_update(varargin{:});
      obj.log('@','in','id',id,'varargin',varargin,'start',obj.start,'stop',obj.stop)
%       % load the metadata for this product
%       obj=obj.product_set(id,varargin{:});
      %retrieve product info
      product=obj.product_get(id,varargin{:});
      obj.log('@','in','product',product,'product start',product.start,'product stop',product.stop)
      % make sure all data sources are loaded
      obj=obj.init_sources(product,varargin{:});
      % get init method
      ih=str2func(product.mdget('method'));
      %resolve leafs: either from level-wrapping, from sources or from existing leafs
      if product.is_wrapped        
        %unwarp products in this list and feed output to input, to unwrap multitple wrapped parts
        product_list=dataproduct.unwrap_product({product});
        %maybe need to prepend some source fields
        if product.ismdfield('source_fields_from')
          %get field path to prepend (from specified source product)
          prepend_product=product_from_source_leafs(obj,{product});
          %make room for outputs
          product_list_out=cell(numel(product_list),numel(prepend_product));
          %prepend the field path of the prepend_product to all elements of the product list
          for i=1:numel(product_list)
            for j=1:numel(prepend_product)
              product_list_out{i,j}=product_list{i};
              product_list_out{i,j}.dataname=product_list{i}.dataname.prepend_field_root(...
                prepend_product{j}.dataname.field_path...
              );
            end
          end
          %assign outputs
          product_list=product_list_out(:);
        end
      elseif product.mdget('explicit_fields','default',false) %|| product.mdget('plot_product','default',false) (this breaks expanding plot product fields)
        %the handling of the data flow between source products and this product is handled explicitly in the product method
        product_list={product};
      elseif product.nr_sources>0
        %expand to source leafs (avoids having to declare 'levelX_name/vals' in all downstream products)
        product_list=obj.product_from_source_leafs({product});
      else
        %enforce sub-category, if requested
        if ~isempty(p.Results.subcat)
          product.dataname=product.dataname.append_field_leaf(p.Results.subcat);
        end
        %expand existing leafs in this product (should generally return itself, since this product has not yet been initialized)
        product_list=product.field_path_expand(obj.data_list(product));
      end
      %sort product_list (needed to resolve save_idx below)
      [~,sort_idx]=sort(cellfun(@(i) i.name,product_list,'UniformOutput',false));
      product_list=product_list(sort_idx);
      %compute indexes of the products that are to be saved (this is needed so that a product with multiple field_paths gets saved in one file)
      save_idx=true(size(product_list));
      for i=1:numel(save_idx)-1
        if strcmp(product_list{i}.name,product_list{i+1}.name)
          save_idx(i)=false;
        end
      end
      %wrapped products need to be loaded differently because they are all stored in the same file,
      %i.e. can't iterate over product_list to load them (but need to iterate to init them)
      if product.is_wrapped && ~p.Results.force
        %try loading this product
        [obj,success]=obj.load(product,varargin{:});
        %if that works, then empty product_list to skip the loading/init loop
        if success
          obj.log('@','iter','loaded wrapped product',product)
          product_list={};
        end
      end
      %loop over all products
      for i=1:numel(product_list)
        %check if data is already loaded
        if obj.isdata_empty(product_list{i}) || p.Results.reload
          %check if init method is to be called even if data was already saved as mat file
          if p.Results.force
            success=false;
          else
            %try to load saved data
            [obj,success]=obj.load(product_list{i},varargin{:});
          end
          %check if data was not loaded
          if success
            obj.log('@','iter',['loaded product_list{',num2str(i),'}'],product_list{i})
          else
            % invoke init method, for all unwrapped leaf products (if any)
%             try
             obj=ih(obj,product_list{i},varargin{:});
%             catch ME
%               if strcmp( ME.identifier,'MATLAB:UndefinedFunction') && ...
%                  str.contains(ME.message,func2str(ih)) && ...
%                  strcmp(cells.first(split(func2str(ih),'.')),'datastorage')
%                 com=['obj.',cells.last(split(func2str(ih),'.')),'(product_list{i},varargin{:})'];
%                 disp(com)
%                 obj=eval(com);
%                 disp('done!')
%               else
%                 error(ME.message)
%               end    
%             end
            obj.log(...
              '@','iter',...
              'applied method',product.mdget('method'),...
              ['to load product_list{',num2str(i),'}'],product_list{i}...
            )
            % save data if this product has been completed (plot_product check done inside)
            if save_idx(i)
              obj.save(product_list{i},varargin{:});
            end
            % update start/stop
            % NOTE: this has to come after saving so that data that is saved in long chunks (i.e. monthly) can be effectively
            %       saved.
            switch product.mdget('method')
            case 'csr.estimate_temp_corr'
              %TODO: maybe product.inclusive can help with this
              obj.log('@','save: skipped saving because this product truncates data internally.')
            otherwise
              if ~product.mdget('plot_product','default',false)
                obj=obj.startstop_retrieve_update(product_list{i});
              end
            end
          end
        else
          obj.log('@','iter','already loaded',['product_list{',num2str(i),'}=',product_list{i}])
        end
      end
      %user feedback
      if obj.debug
        obj.peek
      else
        str.say('initialized',product.str)
      end
      %check if this product is too old for the requested start/stop dates
      %NOTICE: the == test is needed to handle static models (which have the same start/stop)
      assert(obj.start<=obj.stop,['Requested start/stop dates are incompatible with product ',product.str,...
        '; consider updating the data in that product.'])
      obj.log('@','out','id',id,'start',obj.start,'stop',obj.stop)
    end
    function source_list=source_leafs(obj,id_list)
      if iscell(id_list)
        id_list=cells.flatten(cellfun(@obj.source_leafs,id_list,'UniformOutput',false));
        return
      end
      obj.log('@','in','id_list',id_list)
      %type conversion
      product_list=obj.product_get(id_list);
      %easier names
      field_path_now=product_list.dataname.field_path;
      n=numel(field_path_now);
      %results containers
      source_inventory=cell(1,product_list.nr_sources);
      %go over all sources to collect their leafs
      for i=1:product_list.nr_sources
        %search for existing data: start by checking the complete field_path (k=n-1,j=1:1,field_path_now(1:n))
        for k=n-1:-1:0
          for j=1:n-k
            source_inventory{i}=obj.data_list(...
              product_list.sources(i).set_field_path(field_path_now(j:j+k)),...
            'only_non_empty',true);
            %update bail flag
            bail_flag=~isempty(source_inventory{i});
            %are we there yet?
            if bail_flag,break,end
          end
          if bail_flag,break,end
        end
%         %if nothing was found, save all leafs of this source
%         if isempty(source_inventory{i})
%           source_inventory{i}=obj.data_list(product_list.sources(i));
%         end
      end
      %convert to product
      source_list=cellfun(@(i) dataproduct(i), cells.flatten(cells.rm_empty(source_inventory)),'UniformOutput',false);
      obj.log('@','out','source_list',source_list)
    end
    function product_list=product_from_source_leafs(obj,product_list)
      if iscell(product_list)
        product_list=cells.flatten(cellfun(@obj.product_from_source_leafs,product_list,'UniformOutput',false));
        return
      end
      obj.log('@','in','product_list',product_list)
      %ensure existence of sources
      assert(product_list.nr_sources>0,'need positive number of sources')
      %go over all sources and collect their leafs, store each list of datanames in a cell array
      source_inventory=cell(1,product_list.nr_sources);
      for i=1:product_list.nr_sources
        source_inventory{i}=obj.data_list(...
          product_list.sources(i).set_field_path(product_list.dataname.field_path)...
        );
      end
      obj.log('@','1','source_inventory',source_inventory)
      %some products simply copy the field paths from sources
      if product_list.ismdfield('source_fields_from')
        source_idx=find(strcmp(...
          cellfun(@(i) i.name, product_list.source_list,'UniformOutput',false),...
          product_list.mdget('source_fields_from')...
        ));
        %make sure source_fields_from points to a valid source
        assert(~isempty(source_idx),['The value of the metadata field ''source_fields_from'' (',...
          product_list.mdget('source_fields_from'),') is not part of the ''sources''.'])
        % limit the depth of the products
        if product_list.ismdfield('source_fields_max_depth')
          source_depth=product_list.mdget('source_fields_max_depth');
        else
          source_depth=-1;
        end
      %other products have explicit field paths
      elseif product_list.ismdfield('source_field_path')
        %TODO: not sure how to handle this...
        keyboard
      %most products operate on a one-to-one basis in the field paths
      else
        %get strings with field paths for all sources
        field_path_str=...
          cellfun(@(i) ...
            cellfun(@(j) ...
              j.field_path_str,...
            i,'UniformOutput',false),...
          source_inventory,'UniformOutput',false);
        %make sure all sources have the same leafs (otherwise can't really do anything)
        for i=2:numel(field_path_str)
          if ~cells.isequal(field_path_str{1},field_path_str{i})
            disp(['WARNING: ',...
              'The data names under ',product_list.sources(1).str,...
              ' (i.e. ',strjoin(field_path_str{1},', '),') ',...
              ' are in not agreement with those under ',product_list.sources(i).str,...
              ' (i.e. ',strjoin(field_path_str{i},', '),').',...
            ])
            %get common field_path
            field_path_str{1}=intersect(field_path_str{1},field_path_str{i});
            field_path_str{i}=field_path_str{1};
          end
        end
        %check if field paths have vanished
        assert(all(~cells.isempty(field_path_str)),'Found empty field paths, possibly these sources are not compatible.')
        %any source will do and do not limit the depth of products
        source_idx=1;
        source_depth=-1;
      end
      obj.log('@','2','source_idx',source_idx,'source_depth',source_depth)
      %reduce source inventory
      source_inventory=source_inventory{source_idx};
      %retrieve only the required depth (this is untested for source_depth>1
      if source_depth>0
        source_inventory_out={source_inventory{1}.set_field_path({})};
        for i=1:source_depth
          %get unique field names for this depth
          source_depth_now=unique(cellfun(@(dn) dn.field_path{i},source_inventory,'UniformOutput',false));
          %save current source inventory
          source_inventory_now=source_inventory_out;
          %make room for updated source inventory
          source_inventory_out=cell(1,numel(source_depth_now)*numel(source_inventory_now));
          %loop over all sources at current depth and number of current sources in inventory
          for j=1:numel(source_depth_now)
            for k=1:numel(source_inventory_now)
              source_inventory_out{k+numel(source_inventory_now)*(j-1)}=...
                source_inventory_now{k}.append_field_leaf(source_depth_now{j});
            end
          end
        end
        %propagate data
        source_inventory=source_inventory_out;
      end
      %convert to dataproduct
      product_list=product_list.field_path_expand(source_inventory);
      obj.log('@','out','product_list',product_list)
    end
    function obj=init_sources(obj,product,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',        @(i) isa(i,'dataproduct'));
      p.addParameter('reload', false, @(i) islogical(i) && isscalar(i));
      p.addParameter('force',  false, @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(product,varargin{:});
      %loop over all source data
      for i=1:product.nr_sources
        %load this source if it is empty (use obj.init explicitly to re-load or reset data)
        if obj.isdata_empty(product.sources(i)) || p.Results.reload || p.Results.force
          obj.log('@','iter','product',product,'loading source',product.sources(i))
          obj=obj.init(product.sources(i),varargin{:});
        end
      end
    end
    function obj=init_nrtdm(obj,product,varargin)
      obj.log('@','in','product',product,'varargin',varargin,'start',obj.start,'stop',obj.stop)
      % sanity
      assert(time.isvalid(obj.start) && time.isvalid(obj.stop),'Need valid obj.start and obj.stop.')
      %retrieve relevant parameters
      sats         =product.mdget('sats');
      indir        =product.mdget('nrtdm_data_dir');
      nrtdm_sats   =product.mdget('nrtdm_sats');
      nrtdm_product=product.mdget('nrtdm_product');
      %sanity
      str.sizetrap(sats,nrtdm_sats)
      %loop over the satellites
      for s=1:numel(sats)
        %a little tweaking
        if dateshift(obj.stop,'start','day')==obj.stop
          stop_now=dateshift(obj.stop,'end','day');
        else
          stop_now=obj.stop;
        end
        %load the data
        tmp=nrtdm([nrtdm_sats{s},'_',nrtdm_product],obj.start,stop_now,'data_dir',indir);
        %get and save the data
        obj=obj.data_set(product.dataname.set_field_path(sats(s)),tmp.ts);
      end
        obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    %% plot utils
    function out=plotdatafilename(obj,product)
      out=cells.scalar(product.mdset('storage_period','direct').file('plot',...
        'start',obj.start,'stop',obj.stop,...
        'ext','mat',...
        'sub_dirs','single',...
        'start_timestamp_only',false,...
        'suffix','plot-data'...
      ),'get');
    end
    function out=default_plot_columns(obj,product)
      %gather inputs
      v=varargs.wrap('sources',{...
        {...
          'plot_columns',[],@(i) isnumeric(i);...
        },...
        product.plot_args....
      });
      %need to columns to plot
      if isempty(v.plot_columns)
        out=~obj.data_get_scalar(product).iszero_cols;
        for i=1:product.nr_sources
          out=out & ~obj.data_get_scalar(product.sources(i)).iszero_cols;
        end
        out=find(out);
      else
        %NOTICE: this is literal, you need to know the order of the columns, no translation is done.
        out=v.plot_columns;
      end
    end
    function out=default_plot_column_names(obj,product)
      %gather inputs
      v=varargs.wrap('sources',{...
        {...
          'plot_column_names',[],@(i) isnumeric(i);...
        },...
        product.plot_args....
      });
      %need the labels of the columns
      if isempty(v.plot_column_names)
        if product.nr_sources>0
          out=obj.data_get_scalar(product.sources(1)).labels;
        else
          out=obj.data_get_scalar(product).labels;
        end
      else
        %NOTICE: this is usually only good for a quickfix, you first plot the data, then decide to change their name and
        %        define this in the metadata, in the same order the data is plotted. The replacement is done blindly.
        out=v.plot_column_names;
      end
    end
    %% plot functions
    %plots a single data entry (with field paths resolved elsewhere), honours plot_columns (defaults to all non-zero)
    %does not save or annotate the plot in any way, just plots the data using the object plot method.
    function out=justplot(obj,dn,varargin)
      obj.log('@','in','dn',dn,'start',obj.start,'stop',obj.stop)
      %type convertion
      [product,dn]=obj.product_get(dn,varargin{:});
      %gather inputs
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'dn_reference',       dn, @(i) isa(i,'datanames'); ...
          'plot_psd',        false, @(i) islogical(i);...
          'plot_title',   dn.title, @(i) ischar(i); ...
          'plot_columns', obj.default_plot_columns(product), @(i) isnumeric(i);...
        },...
        product.plot_args....
      },varargin{:});
      %if dn_reference is not dn, parse plot arguments from that product
      if ~strcmp(v.dn_reference.str,dn.str)
        v=v.join(obj.product_get(v.dn_reference,varargin{:}).plot_args).join(varargin);
      end
      %retrieve the requested data
      d=obj.data_get_scalar(dn);
      %transmute to frequency series object if requested
      if v.plot_psd; d=simplefreqseries.transmute(d); end
      %get plot arguments: those that start with plot_* and translate them to remote that bit from their name
      plot_args=dataproduct.parse_commands(structs.fieldname_strip(structs.filter(v,'plot_'),'plot_'));
      %further strip the class of the data to be plotted: allows for things like plot_gravity_method
      plot_args=structs.fieldname_strip(plot_args,[class(d),'_']);
      %inform user of the columns to be plotted
      obj.log('@','iter','plot_columns',plot_args.columns)
      %checking data class
      switch class(d)
      case {'simpletimeseries','gravity'}
        out=d.plot(plot_args);
      case 'simplefreqseries'
        out=d.plot_psd(plot_args);
      otherwise
        if isempty(d)
          str.say('Skip plotting empty data',dn.name)
          out=[];
        else
          error([mfilename,': cannot plot data of class ',class(d),'; implementation needed!'])
        end
      end
      obj.log('@','iter','plot_legend',out.legend)
      obj.log('@','out','dn',dn,'start',obj.start,'stop',obj.stop)
    end
    %plots a single data entry (i.e. field paths separately), honours plot_columns (defaults to all non-zero)
    function out=plot_single(obj,dn,varargin)
      obj.log('@','in','dn',dn,'start',obj.start,'stop',obj.stop)
      %parse mandatory args
      dn_list=obj.data_list(dn);
      % recursive call multiple plots
      if numel(dn_list)>1
        out=cellfun(@(i) obj.plot(i,varargin{:}),dn_list,'UniformOutput',false);
        return
      end
      %save the single entry in dataname, convert to product (i.e. load metadata)
      product=obj.product_get(dn,varargin{:});
      %sanity
      assert(~product.metadata.plot_product,['Product ''',product.str,...
        ''' is a plot product, so the ''init'' command is needed, not the ''plot'' command.'])
      %gather inputs
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_file_suffix',      '',@(i) ischar(i); ...
          'plot_file_prefix',      '',@(i) ischar(i); ...
        },...
        product.plot_args...
      },varargin{:});
      %plot filename arguments
      filename_args=[product.file_args('plot'),{...
        'add_field_path',true,...
        'start',obj.start,...
        'stop',obj.stop,...
        'timestamp',true,...
        'prefix',v.plot_file_prefix,...
        'suffix',v.plot_file_suffix...
      }];
      %plot filename
      filename=dn.file(filename_args{:});
      % check if plot is already there
      if isempty(dir(filename))
        plotting.figure(v.varargin{:});
        %retrive data names to be plotted here
        out=obj.justplot(dn,v.varargin{:});
        %annotate plot
        out=structs.copy(plotting.enforce(v.varargin{:}),out);
        out.filename=filename;
        %save this plot
        saveas(gcf,filename)
        str.say('Created plot',filename)
      else
        out=[];
        str.say('Skipped plot',filename)
      end
      obj.log('@','in','dn',dn,'start',obj.start,'stop',obj.stop)
    end
    % plot columns separately for all products given in dn_list, honours plot_columns (defaults to all non-zero)
    function out=plot_separate(obj,dn_list,varargin)
      obj.log('@','in','dn_list',dn_list,'start',obj.start,'stop',obj.stop)
      %parse mandatory args
      dn_list=obj.data_list(dn_list);
      %get produce where the plot preferences are defined
      v=varargs.wrap('sources',{...
        {...
          'dn_reference',dn_list{1},@(i) isa(i,'datanames') || isa(cells.scalar(obj.data_list(i),'get'),'datanames'); ...
        },...
      },varargin{:});
      %plot arguments come from dn_reference (defaults to fist product in dn_list)
      product=obj.product_get(v.dn_reference,varargin{:});
      %sanity
      assert(~product.metadata.plot_product,['Product ''',product.str,...
        ''' is a plot product, so the ''init'' command is needed, not the ''plot'' command.'])
      %gather more inputs (cannot discard cv.dn_reference);
      %NOTICE: varargs.append ignores common entries, so preference is given to what comes first (unlike varargs.wrap)
      v=v.append(product.plot_args).append(....
        varargs({...
          'plot_order',   1:numel(dn_list), @(i) isnumeric(i) && numel(i)==numel(dn_list);...
          'plot_columns',      obj.default_plot_columns(     product), @(i) isnumeric(i);...
          'plot_column_names', obj.default_plot_column_names(product), @(i) iscellstr(i);...
          'plot_file_prefix',      '', @(i) ischar(i); ...
          'plot_file_suffix',      '', @(i) ischar(i); ...
          'plot_title_suffix',     '', @(i) ischar(i); ...
        })).append(plotting.default);
      %make room for handles
      out=cell(1,numel(v.plot_columns));
      %loop over all data columns to plot
      for i=1:numel(v.plot_columns)
        %plot filename arguments
        filename_args=[product.file_args('plot'),{...
          'add_field_path',false,...
          'start',obj.start,...
          'stop',obj.stop,...
          'timestamp',true,...
          'prefix',v.plot_file_prefix,...
          'suffix',{strrep(v.plot_column_names{v.plot_columns(i)},' ','_'),v.plot_file_suffix}... %datanames.file accepts cellstr
        }];
        %plot filename
        filename=product.dataname.file(filename_args{:});
        % check if plot is already there
        if isempty(dir(filename))
          plotting.figure(v.varargin{:});
          out{i}.filename=filename;
          %inites
          out{i}.lines=cell(1,numel(v.plot_order));
          %loop over all data
          for j=1:numel(v.plot_order)
            %retrive data names to be plotted here
            out{i}.lines{j}=obj.justplot(dn_list{v.plot_order(j)},v.varargin{:},...
              'plot_columns',v.plot_columns(i),...
              'plot_title_suffix',[v.plot_column_names{v.plot_columns(i)},' ',v.plot_title_suffix]...
            );
          end
          %resolve strings
          title_str=strjoin(datanames.common(dn_list),' ');
          legend_str=cellfun(@(i) strjoin(i,' '),datanames.unique(dn_list),'UniformOutput',false);
          %patch empty legend entries (this expects there to be only one empty legend entry) 
          if any(cells.isempty(legend_str))
            legend_str(cells.isempty(legend_str))={title_str};
          end
          %annotate plot
          out{i}=structs.copy(...
            plotting.enforce(v.varargin{:},...
              'plot_legend',legend_str,...
              'plot_title',title_str...
            ),...
          out{i});
          %save this plot
          saveas(gcf,filename)
          str.say('Created plot',filename)
        else
          out{i}=[];
          str.say('Skipped plot',filename)
        end          
      end
      %done
      obj.log('@','out','dn_list',dn_list,'start',obj.start,'stop',obj.stop)
    end
    %% legacy plotting (maybe useful to resurect some of this)
    function v=plot_legend(~,h,dn_list,v)
      %get particular defaults from 
      plot_default=varargs(plotting.default).pluck({'plot_zeromean','plot_scale_legend_str'}).cell;
      %v receives the new entries in obj_new. Common entries are ignored.
      v=v.append([{...
        'plot_legend',           {}, @(i) iscellstr(i);...
        'plot_legend_prefixes',  {}, @(i) iscellstr(i);...
        'plot_legend_fontname','Cambria',@(i) ischar(i);...
        'plot_normalize',       false,@(i) islogical(i);...TODO: this needs to be moved to simpledata.plots
        };plot_default]);
      %if there's only one data source, then leave legend as it is
      if numel(dn_list)==1
        legend_str=get(legend,'String');
      else
        %get number of plotted lines
        data_width=sum(cellfun(@(i) numel(i.line_handle),h));
        %add legend only if there are multiple lines
        if data_width>1
          %add the legend given as input, if there
          if ~isempty(v.plot_legend)
            %some sanity
            assert(data_width==numel(v.plot_legend),[...
              'The number of legend entries (',num2str(numel(v.plot_legend)),...
              ') is not in agreement with the number of plotted lines (',num2str(data_width),').'...
            ])
            %propagate
            legend_str=v.plot_legend;
          else
            if ~isempty(v.plot_legend_prefixes)
              prefixes=v.plot_legend_prefixes;
            else
              %get unique parts in datanames
              prefixes=datanames.unique(dn_list);
              %join the dataname parts
              prefixes=cellfun(@(i) strjoin(i,' '),prefixes,'UniformOutput',false);
            end
            %pad with blanks
            prefixes=str.cellpad(prefixes);
            % paranoid sanity
            str.sizetrap(h,prefixes)
            %get how many lines have been plotted
            n=0;
            for j=1:numel(h)
              if ~isempty(h{j}) && isfield(h{j}, 'handle')
                n=n+numel(h{j}.line_handle);
              end
            end
            %make room for legend strings
            legend_str=cell(1,n);
            %loop over all legend entries to add y_mean; TODO: this needs to be moved to simpledata.plots
            c=0;
            for j=1:numel(h)
              if ~isempty(h{j})
                for k=1:numel(h{j}.line_handle)
                  %add y_mean if lines have zero mean
                  if isfield(h{j}, 'y_mean') && v.plot_zeromean %TODO: this needs to be moved to simpledata.plots
                    y_mean=num2str(h{j}.y_mean{k},'%+.3g');
                    v.plot_legend_fontname='FixedWidth';
                  else
                    y_mean='';
                  end
                  %build legend strings
                  c=c+1;
                  legend_str{c}=strtrim([prefixes{j},' ',y_mean]);
                end
              end
            end
            c=0;
            %loop over all legend entries to add y_scale; TODO: this needs to be moved to simpledata.plots
            for j=1:numel(h)
              if ~isempty(h{j})
                for k=1:numel(h{j}.line_handle)
                  %add scale if lines have been normalized
                  if isfield(h{j}, 'y_scale') && v.plot_normalize
                    y_scale=[v.plot_scale_legend_str,num2str(h{j}.y_scale{k},3)];
                  else
                    y_scale='';
                  end
                  %build legend strings
                  c=c+1;
                  legend_str{c}=strtrim([legend_str{c},' ',y_scale]);
                  v.plot_legend_fontname='FixedWidth';
                end
              end
            end
          end
        end
        %save the legend 
        v.plot_legend=legend_str;
      end
    end
    function v=plot_title( ~,~,dn_list,v)
      v=v.append(varargs(plotting.default).pluck({...
        'plot_title',...
      }));
      %handle title keywords
      switch lower(v.plot_title)
      case {'auto'}
        %get title parts
        v.plot_title=strjoin(datanames.common(dn_list),' ');
      end
    end
    function v=plot_ylabel(~,h,~,      v)
      v=v.append(varargs(plotting.default).pluck({...
        'plot_ylabel'...
      }));
      %add the y-label given as input, if there
      if ~isempty(v.plot_ylabel)
        ylabel_str=v.plot_ylabel;
      elseif all(cellfun(@(i) ~isfield(i,'y_units'),h))
        %do nothing
        ylabel_str='';
      else
        %get non-empty indexes
        good_idx=find(cellfun(@(i) ~isempty(i) && ~isempty(i.y_units),h));
        %check if the labels of all lines are compatible
        for i=2:numel(good_idx)
          if ~strcmp(h{good_idx(1)}.y_units,h{good_idx(i)}.y_units)
            error([mfilename,':BUG TRAP: y-units are not consistent in all plotted lines.'])
          end
        end
        %fix y-axis label
        if numel(good_idx)==1
          ylabel_str=h{good_idx(1)}.ylabel;
        else
          ylabel_str=h{1}.y_units;
        end
      end
      v.plot_ylabel=ylabel_str;
    end
    function out=plot_annotate(obj,h,dn,dn_list,varargin)
      %parse mandatory args
      assert(iscell(h),['can only handle input ''h'' as cell, not of class ',class(h),'.'])
      % paranoid sanity
      str.sizetrap(h,dn_list)
      %get product
      product=obj.product_get(dn,varargin{:});
      %gather inputs
      v=varargs.wrap('sources',{product.plot_args,varargin});
      %call annotation methods
      v=obj.plot_legend(h,dn_list,v);
      v=obj.plot_title( h,dn_list,v);
      v=obj.plot_ylabel(h,dn_list,v);
      %enforce plot preferences (using the metadata of the current dataname)
      out=product.enforce_plot(v.varargin{:});
    end
    %produces a structure with fields:
    % - sources: cell of cellstr with product leafs
    % - source_names: cells.flatten(datanames.unique(product.sources))
    % - plot_columns: indices of the columns to plot
    % - plot_column_names: labels of the columns:
    %    - in product, if there's no sources, or
    %    - in the first source, if there are sources
    % - startlist: datetime vectors with periodic starts
    % - soptlist : datetime vectors with periodic stops
    function e=plot_elements(obj,product,varargin)
      %gather inputs
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_source_idx',   [],@(i) isnumeric(i); ...
          'plot_source_names', {},@(i) iscellstr(i); ...
          'plot_columns',      [],@(i) isnumeric(i);...
          'plot_column_names', {},@(i) iscellstr(i); ...
        },...
        product.plot_args....
      },varargin{:});
      %need the source list
      e.sources=obj.source_leafs(product);
      if ~isempty(v.plot_source_idx); e.sources=e.sources(v.plot_source_idx);end
      %need the source list names
      if isempty(v.plot_source_names)
        e.source_names=cells.flatten(datanames.unique(product.source_list));
      else
        e.source_names=v.plot_source_names;
      end
      %enforce plot_source_idx
      if ~isempty(v.plot_source_idx); e.source_names=e.source_names(v.plot_source_idx);end
      %need to columns to plot
      if isempty(v.plot_columns)
        e.plot_columns=~obj.data_get_scalar(product).iszero_cols;
        for i=1:product.nr_sources
          e.plot_columns=e.plot_columns & ~obj.data_get_scalar(product.sources(i)).iszero_cols;
        end
        e.plot_columns=find(e.plot_columns);
      else
        %NOTICE: this is literal, you need to know the order of the columns, no translation is done.
        e.plot_columns=v.plot_columns;
      end
      %need the labels of the columns
      if isempty(v.plot_column_names)
        if product.nr_sources>0
          e.plot_column_names=obj.data_get_scalar(product.sources(1)).labels;
        else
          e.plot_column_names=obj.data_get_scalar(product).labels;
        end
      else
        %NOTICE: this is usually only good for a quickfix, you first plot the data, then decide to change their name and
        %        define this in the metadata, in the same order the data is plotted. The replacement is done blindly.
        e.plot_column_names=v.plot_column_names;
      end
      %gather list of daily/monthly/yearly/... data files
      [~,e.startlist,e.stoplist]=product.file('data',v.varargin{:},'start',obj.start,'stop',obj.stop);
    end
    %plots multiple data entries (more limited than datastorage.plot but flexible for multiple data sources)
    % - does not save plots in any way, that has to be done externally;
    % - does not use the 'plot_elements' method
    function [obj,out]=plot_mult(obj,id,dn_list,plot_columns,varargin)
      obj.log('@','in','dn',id,'dn_list',dn_list,'plot_columns',plot_columns,'start',obj.start,'stop',obj.stop)
      %type conversion
      [product,dn]=obj.product_get(id,varargin{:});
      %gather inputs
      dn_list=datanames.array(dn_list,dn.field_path);
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_order', [1:numel(dn_list),0], @(i) isnumeric(i) && numel(i)==numel(dn_list)+1;... (index 0 represents dn, non-zero represent dn_list)
          'plot_normalize',       false,@(i) islogical(i);...TODO: this needs to be moved to plotting.enforce
        },...
        product.plot_args...
      },varargin{:});      
      %easier names
      plot_scale_default=ones(...
        max(cellfun(@(i) i.width,obj.data_get(dn_list))),...
        numel(v.plot_order)...
      );
      %derived inputs: v receives the new entries. Common entries are ignored.
      v=v.append(varargs({...
        'plot_scale',      plot_scale_default,       @(i) all(size(i)==size(plot_scale_default));...
        'plot_smooth_span',zeros(size(v.plot_order)),@(i) all(size(i)==size(v.plot_order));...
      }));
      %update plot_columns in v
      if v.isparameter('plot_columns'); v.plot_columns=plot_columns; end
      %assign ordered product list
      dn_plot_list=cell(size(v.plot_order));
      zero_idx=v.plot_order==0;
      assert(sum(zero_idx)<=1,['Expecting ''plot_order'' to have only one zero entry, not ',num2str(sum(zero_idx)),'.'])
      dn_plot_list(~zero_idx)=dn_list(v.plot_order(~zero_idx));
      dn_plot_list( zero_idx)={dn};
      %get rid of zero entry if this is a plot product and it is part of the plot_order
      if v.plot_product
        dn_plot_list( zero_idx)=[];
      end
%       for i=1:numel(v.plot_order)
%         if v.plot_order(i)==0
%           dn_plot_list{i}=dn;
%         else
%           dn_plot_list{i}=dn_list{v.plot_order(i)};
%         end
%       end
      %expand scalar plot_columns into cell array (usually plot_columns is not a scalar)
      plot_columns=cells.deal(plot_columns,size(dn_plot_list));
%       err=false;
%       switch numel(plot_columns) %this is not v.plot_columns!
%         case 0
%           plot_columns=num2cell(ones(1,numel(dn_plot_list)));
%         case 1
%           if isnumeric(v.plot_columns)
%             plot_columns=num2cell(plot_columns*ones(1,numel(dn_plot_list)));
%           elseif iscell(v.plot_columns)
%             plot_columns=num2cell(plot_columns{1}*ones(1,numel(dn_plot_list)));
%           else
%             err=true;
%           end
%         case numel(dn_plot_list)
%           if isnumeric(plot_columns)
%             plot_columns=num2cell(plot_columns);
%           elseif iscell(plot_columns)
%             %do nothing
%           else
%             err=true;
%           end
%         otherwise
%           tmp=cell(size(dn_plot_list));
%           tmp(:)={plot_columns};
%           plot_columns=tmp;
%       end
%       %error handling
%       assert(~err,['Cannot handle input ''plot_columns'' of class ''',class(v.plot_columns),...
%         ''' and/or its length (',num2str(numel(v.plot_columns)),') is in conflict with length ',...
%         'of input ''dn_list'' (',num2str(numel(dn_plot_list)),').'])
      %expand smooth_span (only done if needed and sizes are checked)
      v.plot_smooth_span=cells.c2m(cells.deal(v.plot_smooth_span,size(v.plot_order)));
      %make room for outputs
      h=cell(size(dn_plot_list));
      %loop over all datanames
      for i=1:numel(dn_plot_list)
        h{i}=obj.justplot(dn_plot_list{i},...
          v.varargin{:},...
          'dn_reference',dn,... %otherwise the justplot method picks the plot_* parameters define in dn_list{i}
          'plot_scale',v.plot_scale(plot_columns{i},i),... %this is not v.plot_columns!
          'plot_smooth_span',v.plot_smooth_span(i),... 
          'plot_columns',plot_columns{i}...
        );
      end
      %sanity
      assert(all(~cells.isempty(h)),['Could not plot data for ',...
        strjoin(cellfun(@(i) i.str,dn_plot_list(cells.isempty(h)),'UniformOutput',false),', ')...
      ])
      %annotate plot
      out=obj.plot_annotate(h,dn,dn_plot_list,v.varargin{:});
      %propagate plot handles
      out.plot_handles=h;
      %reset ylabel to use intersetion of all labels involved (if plot is not normalized)
      if ~v.plot_normalize
        label_str=cells.deal_struct(h,'ylabel');
        label_str=cellfun(@(i) strsplit(i,'\[|\]','DelimiterType','RegularExpression'),label_str,'UniformOutput',false);
        label_str=unique(cellfun(@(i) i{2},label_str,'UniformOutput',false));
        ylabel(strjoin(label_str,' or '))
      end
      %done
      obj.log('@','out','dn',dn,'start',obj.start,'stop',obj.stop)
    end
    %parses the sources of a product and calls plot_mult (useful as value of the 'method' metadata field)
    function [obj,out]=plot_auto(obj,id,varargin)
      obj.log('@','in','id',id,'start',obj.start,'stop',obj.stop)
      %type conversion
      [product,dn]=obj.product_get(id,varargin{:});
      %sanity
      assert(obj.isdata_leaf(dn),'Can only handle leaf products.')
      %gather inputs
      v=varargs.wrap('sources',{...
        plotting.default,...
        {...
          'plot_columns',          -1,@(i) isnumeric(i);...
          'plot_together',       {''},@(i) iscellstr(i);...
        },...
        product.plot_args...
      },varargin{:});
      %retrieve plot elements (repetitive processing of parameters)
      e=obj.plot_elements(product,v.varargin{:});
      %inform user if nothing to plot
      if isempty(e.startlist)
        str.say('Nothing to plot for',product,'start =',obj.start,'stop =',obj.stop,v.varargin{:})
      end
      %make room for outputs
      out=cell(size(v.plot_columns));
      %loop over all columns to plot
      for c=1:numel(v.plot_columns)
        obj.log('@','iter','column',c,'start',obj.start,'stop',obj.stop)
        %plot filename arguments
        filename_args=[product.file_args('plot'),{...
          'start',obj.start,...
          'stop', obj.stop,...
          'timestamp',true,...
          'remove_part','',...
          'prefix',v.plot_file_prefix...
          'suffix',strjoin([product.dataname.field_path,e.plot_column_names(v.plot_columns(c)),{v.plot_title_suffix}],'.')...
        }];
        %plot filename
        filename=dn.file(filename_args{:});
        if isempty(dir(filename))
          %make sure there is data
          if any(cell2mat(obj.vector_method_tr('all','nr_valid'))>1)
            %plot it
            plotting.figure(v.varargin{:});
            [~,out{c}]=obj.plot_mult(dn,product.source_list,v.plot_columns(c),...
              v.varargin{:},...
              'plot_title',strjoin([product.dataname.field_path,e.plot_column_names(v.plot_columns(c)),{v.plot_title_suffix}],' ')...
            );
            out{c}.filename=filename;
            out{c}.fig_handle=gcf;
            %check if any data was plotted
            if all(isempty(out{c}));   str.say('Skipped plot',filename,'(no data plotted)'); close(gfc) %nothing plotted
            else saveas(gcf,filename); str.say('Created plot',filename)                                 %save this plot
            end
          else str.say('Skipped plot',filename,['(no data in ',product.name,')']); out{c}={};
          end
        else   str.say('Skipped plot',filename,'(file already exists)');           out{c}={};
        end
      end
      obj.log('@','out','id',id,'start',obj.start,'stop',obj.stop)
    end
    %% generalized operations
    function obj=stats(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      %sanity
      assert(product.nr_sources==1,['Can only handle one source model, not ',num2str(product.nr_sources),'.'])
      assert(obj.isdata_leaf(product.sources(1)),'Can only handle leaf products.')
      %retrieve data
      d=obj.data_get_scalar(product.sources(1));
      assert(structs.respondto(d,'stats'),['Product ',product.sources(1).str,' does not handle the ''stats'' method. Debug needed'])
      %compute stats
      stats_data=d.stats(...
        'period',product.mdget('stats_period'),...
        'overlap',product.mdget('stats_overlap'),...
        'outlier',product.mdget('stats_outlier'),...
        'struct_fields',product.mdget('stats')...
      );
      %save data
      obj=obj.data_set(product,stats_data);
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=arithmetic(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      %gather inputs
      v=varargs.wrap('sources',{...
        {...
          'operation',       'plus',               @(i) ischar(i) || iscellstr(i);...
          'operation_order', 1:product.nr_sources, @(i) isnumeric(i) && numel(i)<=product.nr_sources;...
        },...
        product.metadata,...
      },varargin{:});   
      obj.log('@','in','operation',v.operation)
      %expand the operations, if scalar
      v.operation=cells.deal(cells.scalar(v.operation,'set'),size(v.operation_order));
      %operate on all sources
      for i=1:numel(v.operation_order)
        idx=v.operation_order(i);
        if i==1
          %get first source
          out=obj.data_get_scalar(product.sources(idx));
          msg=[product.str,' = ',product.sources(idx).str];
        else
          %get the current source
          in=obj.data_get_scalar(product.sources(idx));
          msg=[msg,' ',v.operation,' ',product.sources(idx).str]; %#ok<AGROW>
          %enforce common time domain, assume it is defined by the first argument
          in=in.interp(out.t);
          %operate
          out=out.(v.operation{i})(in);
        end
        obj.log('@','iter','operation',msg)
      end
      %propagate result
      obj=obj.data_set(product,out);
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    %TODO: merge this with the above
    function obj=operation(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      %gather inputs
      v=varargs.wrap('sources',{...
        {...
          'operation',       'error',               @(i) ischar(i) || iscellstr(i);...
          'operation_order', 1:product.nr_sources, @(i) isnumeric(i) && numel(i)<=product.nr_sources;...
        },...
        product.metadata,...
      },varargin{:});   
      obj.log('@','in','operation',v.operation)
      %expand the operations, if scalar
      v.operation=cells.deal(cells.scalar(v.operation,'set'),size(v.operation_order));
      %operate on all sources
      for i=1:numel(v.operation_order)
        idx=v.operation_order(i);
        if i==1
          %propagate first source
          obj=obj.data_set(product,obj.data_get(product.sources(idx)));
          msg=[product.str,' = ',product.sources(idx).str];
        else
          %operate on this source
          obj=obj.vector_obj_op3(product,product.sources(idx),v.operation{i},product);
          msg=[msg,' ',v.operation,' ',product.sources(idx).str]; %#ok<AGROW>
        end
        obj.log('@','iter','operation',msg)
      end
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=filter(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      assert(product.nr_sources==1,['The filter method cannot operated on product ''',product.str,...
        ''' because it only accept one source, not ',num2str(product.nr_sources),'.'])
      %retrieve required operation
      operation=product.mdget('operation');
      args=product.mdget('arguments');
      %get source
      in=obj.data_get_scalar(product.sources(1));
      %operate
      out=in.(operation)(args{:});
      %propagate result
      obj=obj.data_set(product,out);
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    %% special operations
    function obj=component_split(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      assert(product.nr_sources==1,['The component_split method cannot operated on product ''',product.str,...
        ''' because it only accept one source, not ',num2str(product.nr_sources),'.'])
      %gather inputs
      v=varargs.wrap('sources',{product.metadata},varargin{:});
      %get source
      in=obj.data_get_scalar(product.sources(1));
      %TEMPORARY: patch descriptor
      in.descriptor=product.sources(1).codename;
      %split into components
      out=in.component_split(v.rename_silent('comp_plot_dir','plot_dir').rename_silent('comp_data_dir','data_dir').varargin{:});
      %append 'signal' to field path, so that downstream procedures can identify this as valid data
      if ~strcmp(product.dataname.field_path{end},'signal'); product.dataname=product.dataname.append_field_leaf('signal'); end
      %propagate result
      obj=obj.data_set(product,out);
      obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)
    end
    function obj=corr(obj,product,varargin)
      error('implementation needed')
      % parse mandatory arguments
      p=inputParser;
      p.addRequired('product', @(i) isa(i,'dataproduct'));
      p.parse(product);
      %sanity
      if product.nr_sources~=1
        error([mfilename,': number of sources in product ',product.dataname,...
          ' is expected to be 1, not ',num2str(product.nr_sources),'.'])
      end
      %retrieve source and part lists
      sourcep=obj.product_get(product.sources(1),varargin{:});
      [~,levels,fields,sats]=sourcep.partslist;
      %check if data is already there
      if ~product.isfile('data')
        %loop over all data
        for i=1:numel(levels)
          for j=1:numel(fields)
            %retrive data for both satellites
            d=obj.field_get(sourcep.dataname.type,levels{i},fields{j});
            %can only deal with two sats at the moment
            if numel(fieldnames(d))~=2
              error([mfilename,': can only deal with two sats at the moment. Implementation needed!'])
            end
            %compute stats
            stats_data=simpletimeseries.stats2(...
              d.(sats{1}),...
              d.(sats{2}),...
              'period',product.mdget('period'),...
              'overlap',product.mdget('overlap'),...
              'outlier',product.mdget('outlier'),...
              'struct_fields',product.mdget('stats')...
            );
            %loop over all requested stats
            stats_fields=fieldnames(stats_data);
            for si=1:numel(stats_fields)
              %save stats in 'sat' names given by the requested stats name
              obj=obj.sat_set(...
                product.dataname.type,levels{i},fields{j},stats_fields{si},...
                stats_data.(stats_fields{si})...
              );
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
        for i=1:numel(levels)
          obj=obj.level_set(product.dataname.type,levels{i},s.(levels{i})); %#ok<NODEF>
        end
      end
    end
%     %% usefull stuff
%     function [sourcep,df]=dataflow(obj,product)
%       %retrive source metadata
%       sourcep=cell(1,product.nr_sources);
%       for i=1:product.nr_sources
%         sourcep{i}=obj.product_get(product.sources(i),varargin{:});
%       end
%       %retrieve dataflow structure
%       df=product.mdget('dataflow');
%       %loop over the types/levels/fields metadata fields
%       for i=1:numel(obj.parts)
%         %easier names
%         partname=[obj.parts{i},'s'];
%         %check if this part is there
%         if isfield(df,partname)
%           %retrive metadata value
%           partvalues=df.(partname);
%           %multiple rows means this product spans multiple types/levels/fields
%           if iscell(partvalues{1})
%             %need to unwrap nested cells case there are multiple rows
%             df.(partname)=cell(numel(partvalues),numel(sourcep)+1);
%             for j=1:numel(partvalues)
%               %sanity on types/levels/fields: number of columns must be equal to number of sources
%               assert(numel(partvalues{j}) == numel(sourcep)+1,[mfilename,': ',...
%                 'ilegal ''',partname,''' entry of the ''dataflow'' metadata field, ',...
%                 'it must have the same number of columns as the number of product sources +1 (',num2str(numel(sourcep)+1),'), ',...
%                 'not ',num2str(numel(partvalues{j})),'.'])
%               df.(partname)(j,:)=partvalues{j};
%             end
%           else
%             %variable partvalues already contains types/levels/fields
%             assert(numel(partvalues) == numel(sourcep)+1,[mfilename,': ',...
%               'ilegal ''',obj.parts{i},'s'' metadata field, ',...
%               'it must have the same number of columns as the number of product sources (',num2str(numel(sourcep)+1),'), ',...
%               'not ',num2str(numel(partvalues)),'.'])
%           end
%         else
%           %if the dataflow structure is missing a part, then patch from this product or sources
%           if product.ismdfield(partname)
%             df.(partname)=product.mdget(partname,'always_cell_array',true);
%             %this avoids having to define repetitive (e.g.) sats values in the 'dataflow' metadata field,
%             %the scalar 'sats' metadata fiels is enough (doesn't work for 'types')
%             if i>1 && numel(df.(partname))==1 && numel(df.([obj.parts{1},'s']))>1
%               partvalues=cell(size(df.([obj.parts{1},'s'])));
%               partvalues(:)=df.(partname);
%               df.(partname)=partvalues;
%             end
%           else
%             found_part=false;
%             %search in all sources for an explicit list of this partname
%             for j=1:product.nr_sources
%               if sourcep{j}.ismdfield(partname)
%                 partvalues=sourcep{j}.mdget(partname);
%                 df.(partname)=cell(numel(partvalues),numel(sourcep)+1);
%                 for k=1:numel(sourcep)+1
%                   df.(partname)(:,k)=partvalues(:);
%                 end
%                 found_part=true;
%               end
%             end
%             if ~found_part
%               %if nothing found, search for values of this partname in the data of the sources
%               partvalues=obj.data_list(sourcep{1}.name);
%               %check if all sources share the same values of this part
%               for j=2:product.nr_sources
%                 %cell array value comparisson:
%                 %http://stackoverflow.com/questions/3231580/matlab-comparison-of-cell-arrays-of-string
%                 assert(cells.isequal(obj.data_list(sourcep{j}.name),partvalues),...
%                   ['The ',partname,' names differ in products ',sourcep{1}.name,' and ',...
%                   sourcep{j}.name,'.'])
%               end
%               %propagate this part values
%               df.(partname)=cell(numel(partvalues),numel(sourcep)+1);
%               for j=1:numel(partvalues)
%                 df.(partname)(j,:)=partvalues(j);
%               end
%             end
%           end
%         end
%       end
%       %return scalar source product if only one
%       if numel(sourcep)==1
%         sourcep=sourcep{1};
%       end
%     end
  end
end
