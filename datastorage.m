classdef datastorage
  %static
  properties(Constant)
    %default value of some internal parameters
    parameter_list=struct(...
      'start',          struct('default',time.zero_date,'validation',@(i) isdatetime(i)),...
      'stop',           struct('default',time.inf_date, 'validation',@(i) isdatetime(i)),...
      'debug',          struct('default',false,         'validation',@(i) islogical(i)),...
      'metadata_dir',   struct('default',dataproduct.default_list.metadata_dir,'validation',@(i) ischar(i))...
    );
  end
  %public
  properties(GetAccess=public,SetAccess=public)
    debug
  end
  %read only
  properties(SetAccess=private)
    par
    data
    metadata_dir
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
    function out=parameters
      out=fieldnames(datastorage.parameter_list);
    end
  end
  methods
    %% constructor
    function obj=datastorage(varargin)
      obj.log('@','in','varargin',varargin)
      % parameter names
      pn=datastorage.parameters;
      % input parsing
      p=inputParser;
      p.KeepUnmatched=true;
      % declare parameters
      for i=1:numel(pn)
        p.addParameter(pn{i},datastorage.parameter_list.(pn{i}).default,datastorage.parameter_list.(pn{i}).validation)
      end
      % parse it
      p.parse(varargin{:});
      % reset data type list
      obj=obj.data_init;
      % save parameters
      for i=1:numel(pn)
        obj.(pn{i})=transpose(p.Results.(pn{i})(:));
        obj.log('@','iter',pn{i},obj.(pn{i}),['p.Results.(',pn{i},')'],p.Results.(pn{i}))
      end
      obj.log('@','out')
      % data is added to this initialized object with the 'init' method (see below)
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
    function dn_list=datanames(obj,non_empty)
      dn_list=datanames.array(fieldnames(obj.data));
      if exist('non_empty','var') && non_empty
        empty_idx=cellfun(@(i) isempty(obj.data.(i.name)),dn_list);
        dn_list=dn_list(~empty_idx);
      end
    end
    function peek(obj,dn)
      if obj.isempty
        return
      end
      if ~exist('dn','var') || isempty(dn)
        dn='all';
      end
      dn=obj.data_list(dn);
      cellfun(@(i) disp(str.tablify([32,6],i.codename,obj.data_get_scalar(i).size)),dn)
    end
    %% value operations (basically wraps some methods in the structs object)
    function out=value_get(obj,dn)
      %reduce dataname to common object
      dn=datanames(dn);
      %retrieve value
      out=structs.get_value(obj.data.(dn.name_clean),dn.field_path);
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
    function out=isdata_leaf(obj,dn)
      %reduce dataname to common object
      dn=datanames(dn);
      %determine leafness
      out=~isfield(obj.data,dn.name_clean) || structs.isleaf(obj.data.(dn.name_clean),dn.field_path);
    end
    function dn_list=data_list(obj,dn)
      if numel(dn)==0
        %trivial call
        dn_list={};
        return
      elseif ischar(dn) && strcmp(dn,'all')
        %handle named inputs
        dn_list=obj.data_list(obj.datanames(true));
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
      if obj.isdata_leaf(dn)
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
      %handle scalar values, unwrap them to all field paths
      if numel(values)==1 && ~iscell(values)
        tmp=values;
        values=cell(size(dn_list));
        values{:}=tmp;
      end
      %make sure fiels paths and values have the same length
      assert(numel(dn_list)==numel(values),'variables ''dn_list'' and ''values'' must have the same length.')
      %propagate
      for i=1:numel(dn_list)
        obj=obj.value_set(dn_list{i},values{i});
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
    % Applies the 2-argument function f to the data given by dataname_list
    function obj=vector_obj_op2(obj,dn,obj2,f,varargin)
      %operate on the requested set
      values1=obj.data_get(dn);
      validv1=cells.rm_empty(values1);
      values2=obj2.data_get(dn);
      validv2=cells.rm_empty(values2);
      %sanity
      if numel(values1) ~= numel(values2)
        error([mfilename,': the given dataname list does not correspond to the same number of data entries in both input obj.'])
      end
      if any(validv1~=validv2)
        error([mfilename,': the given dataname list does not have data in both objects. Debug needed.'])
      end
      %outputs
      result=cell(size(values1));
      %operate
      for i=1:numel(values1)
        result{i}=values1{i}.(f)(values2{i},varargin{:});
      end
      %propagate cell array with results back to object
      obj=obj.data_set(dn,result);
    end
    %% start/stop operations
    function out=start_retrieve(obj,dn)
      out=obj.vector_method_tr(dn,'start');
      %need to filter out zero_dates
      out(~time.isvalid(out))={time.zero_date};
      %take the maximum of all starts
      out=max([out{:}]);
    end
    function out=stop_retrieve(obj,dn)
      out=obj.vector_method_tr(dn,'stop');
      %need to filter out zero_dates
      out(~time.isvalid(out))={time.inf_date};
      %take the minimum of all stops
      out=min([out{:}]);
    end
    function out=get.start(obj)
      out=obj.starti;
    end
    function out=get.stop(obj)
      out=obj.stopi;
    end
    function obj=set.start(obj,start)
      %check if updating is needed
      if isempty(obj.stopi) || start>obj.starti
        %update internal record
        obj.starti=start;
        %trim all data entries
        obj=obj.vector_method_set('all','start',start);
      end
    end
    function obj=set.stop(obj,stop)
      %check if updating is needed
      if isempty(obj.stopi) || stop<obj.stopi
        %update internal record
        obj.stopi=stop;
        %trim all data entries
        obj=obj.vector_method_set('all','stop',stop);
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
    function obj=startstop_retrieve_update(obj,dn)
      obj=obj.startstop_update('start',obj.start_retrieve(dn),'stop',obj.stop_retrieve(dn));
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
    %% length operations
    function out=length(obj,dn)
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
      %loop over all retrieved objects
      for i=1:numel(obj_list)
        disp(obj_list{i}.descriptor)
        obj_list{i}.print(varargin{:});
      end
    end
    %% metadata interface
    function obj=product_set(obj,dn,varargin)
      %reduce dataname to common object
      dn=datanames(dn);
      %set the requested metadata
      obj.par.(dn.name_clean)=dataproduct(dn,...
        'metadata_dir',obj.metadata_dir,... %add more metadata-relevant parameters here
        varargin{:}...
      );
    end
    function md=product_get(obj,dn)
      %reduce dataname to common object
      dn=datanames(dn);
      %some sanity
      assert(isfield(obj.par,dn.name_clean),['cannot find product ',dn.name,'. ',...
          'Have you called obj.init(''',dn.name,''')? ',...
          'Maybe the requested product is not sourced in any of the products already loaded.'])
      %retrieve the requested metadata
      md=obj.par.(dn.name_clean);
    end
    %% load/save operations
    function [obj,success]=load(obj,product,varargin)
      %type sanity
      assert(isa(product,'dataproduct'),['input ''product'' has to be of class ''dataproduct'', not ''',class(product),'''.'])
%       %operate on leaves
%       if ~obj.isdata_leaf(product)
%         product_leafs=product.field_path_expand(obj.data_list(product));
%         success=true(size(product_leafs));
%         for i=1:numel(product_leafs)
%           [obj,success(i)]=obj.load(product_leafs{i},varargin{:});
%         end
%         success=all(success);
%         return
%       end
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % get the file list
      file_list=product.file('data','discover',false,'start',obj.start,'stop',obj.stop,varargin{:});
      % check if any file is missing
      if any(~product.isfile('data','start',obj.start,'stop',obj.stop,varargin{:}))
        % debug output
        success=false;
        if isempty(file_list)
          obj.log('@','out','no files to load for',product)
        else
          obj.log('@','out','file(s) missing for',[product.str,char(10),strjoin(file_list,char(10))])
        end
        return
      end
      obj.log('nr of files to load',numel(file_list))
      %loop over all files
      for f=1:numel(file_list)
        % debug output
        obj.log(['loading file ',num2str(f),' '],[' ',file_list{f}])
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
      end
      % save the data in the object
      obj=obj.data_set(product,s_out);
      % update start/stop
      obj=obj.startstop_retrieve_update(product);
      % debug output
      success=true;
      obj.log('@','out','product',product)
    end
    function save(obj,product,varargin)
      %type sanity
      assert(isa(product,'dataproduct'),['input ''product'' has to be of class ''dataproduct'', not ''',class(product),'''.'])
%       %operate on leaves
%       if ~obj.isdata_leaf(product)
%         product_leafs=product.field_path_expand(obj.data_list(product));
%         for i=1:numel(product_leafs)
%           obj.save(product_leafs{i},varargin{:});
%         end
%         return
%       end
      %ignore plot products, file saving is done inside the init method
      if product.ismd_field('plot_product') && product.mdget('plot_product')
        return
      end
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      % get the file list
      [file_list,startlist,stoplist]=product.file('data',...
        'discover',false,...
        'start',obj.start,...
        'stop',obj.stop,...
        varargin{:}...
      );
      %check if nothing is to be save
      if isempty(file_list)
        obj.log('@','out','no files to save for',product)
        return
      end
      obj.log('nr of files to save',numel(file_list))
      %get the data
      s_out=obj.value_get(product);
      %loop over all files
      for f=1:numel(file_list)
        % debug output
        obj.log(['saving file ',num2str(f)],file_list{f})
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
        %save this section of the data to the file
        save(file_list{f},'s');
      end
      % debug output
      obj.log('@','out','product',product)
    end
    %% datatype initialization
    function obj=init(obj,dn,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('force',    false,@(i) islogical(i) && isscalar(i));
      p.addParameter('recompute',false,@(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(varargin{:});
      %reduce dataname to common object
      dn=datanames(dn);
      %update start/stop if given in inputs
      obj=obj.startstop_update(varargin{:});
      obj.log('@','in','dataname',dn,'varargin',varargin,'start',obj.start,'stop',obj.stop)
      % load the metadata for this product
      obj=obj.product_set(dn,varargin{:});
      %retrieve product info
      product=obj.product_get(dn);
      % make sure all data sources are loaded
      obj=obj.init_sources(product,varargin{:});
      % get init method
      ih=str2func(product.mdget('method'));
      %resolve leafs: either from level-wrapping, from sources or from existing leafs
      if product.is_wrapped
        %unwarp products in this list and feed output to input, to unwrap multitple wrapped parts
        product_list=dataproduct.unwrap_product({product});
      elseif product.nr_sources>0
        %expand to source leafs (avoids having to declare 'levelX_name/vals' in all downstream products)
        product_list=obj.source_leafs({product});
      else
        %expand existing leafs in this product (should generally return itself, since this product has not yet been initialized)
        product_list=obj.product_leafs({product});
      end
      %loop over all products
      for i=1:numel(product_list)
        obj.log(...
          '@','iter',...
          'method',product.mdget('method'),...
          ['product_list{',num2str(i),'}'],product_list{i}...
        )
        %check if data is already loaded
        if obj.isdata_empty(product_list{i}) || p.Results.force
          %check if init method is to be called even if data was already saved as mat file
          if p.Results.recompute
            success=false;
          else
            %try to load saved data
            [obj,success]=obj.load(product_list{i},varargin{:});
          end
          %check if data was not loaded
          if ~success
            % invoke init method, for all unwrapped leaf products (if any)
            obj=ih(obj,product_list{i},varargin{:});
            % update start/stop
            obj=obj.startstop_retrieve_update(product_list{i});
            % save data
            obj.save(product_list{i},varargin{:});
          end
        else
          obj.log('@','iter: data already loaded')
        end
      end
      obj.log('@','out','dataname',dn,'start',obj.start,'stop',obj.stop)
      %user feedback
      if ~obj.debug
        str.say('initialized',product.str)
      end
      %check if this product is too old for the requested start/stop dates
      assert(obj.start<obj.stop,['Requested start/stop dates are incompatible with product ',dn.str,...
        '; consider updating the data in that product.'])
    end
    function out=product_leafs(obj,product_list)
      if iscell(product_list)
        out=cells.flatten(cellfun(@obj.product_leafs,product_list,'UniformOutput',false));
      else
        out=product_list.field_path_expand(obj.data_list(product_list));
      end
    end
    function product_list=source_leafs(obj,product_list)
      if iscell(product_list)
        product_list=cells.flatten(cellfun(@obj.source_leafs,product_list,'UniformOutput',false));
      else
        %branch on the existence of sources
        if product_list.nr_sources>0
          %go over all sources and collect their leafs, store each list of datanames in a cell array
          source_inventory=cell(1,product_list.nr_sources);
          for i=1:product_list.nr_sources
            source_inventory{i}=obj.data_list(...
              product_list.sources(i).set_field_path(product_list.dataname.field_path)...
            );
          end
          %some products simply copy the field paths from sources
          if product_list.ismd_field('source_fields_from')
            source_idx=find(strcmp(...
              cellfun(@(i) i.name, product_list.source_list,'UniformOutput',false),...
              product_list.mdget('source_fields_from')...
            ));
            % limit the depth of the products
            if product_list.ismd_field('source_fields_max_depth')
              source_depth=product_list.mdget('source_fields_max_depth');
            else
              source_depth=-1;
            end
          %most products operate on a one-to-one basis in the field paths
          else
            %get strings with field paths for all source
            field_path_str=cellfun(...
              @(i) cellfun(@(j) j.field_path_str,i,'UniformOutput',false),...
            source_inventory,'UniformOutput',false);
            %make sure all sources have the same leafs (otherwise can't really do anything)
            for i=2:numel(source_inventory)
              assert(cells.isequal(field_path_str{1},field_path_str{i}),[...
                'The data names under ',product_list.sources(1).str,...
                ' (i.e. ',strjoin(field_path_str{1},', '),') ',...
                ' are in not agreement with those under ',product_list.sources(i).str,...
                ' (i.e. ',strjoin(field_path_str{i},', '),').',...
              ])
            end
            source_idx=1;
            source_depth=-1;
          end
          %reduce source inventory
          source_inventory=source_inventory{source_idx};
          %retrieve only the required depth (this is untested for source_depths>1
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
          %expand this product list
          product_list=product_list.field_path_expand(source_inventory);
        else
          %need implementation
          error('implementation needed')
        end
      end
    end
    function obj=init_sources(obj,product,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',@(i) isa(i,'dataproduct'));
      % parse it
      p.parse(product,varargin{:});
      %loop over all source data
      for i=1:product.nr_sources
        %load this source if it is empty (use obj.init explicitly to re-load or reset data)
        if obj.isdata_empty(product.sources(i))
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
        %load the data
        tmp=nrtdm([nrtdm_sats{s},'_',nrtdm_product],obj.start,obj.stop,'data_dir',indir);
        %get and save the data
        obj=obj.data_set(product.dataname.set_field_path(sats(s)),tmp.ts);
      end
        obj.log('@','out','product',product,'start',obj.start,'stop',obj.stop)  
    end
    %% plot utils
    function h=justplot(obj,dn,varargin)
      obj.log('@','in','dn',dn,'start',obj.start,'stop',obj.stop)
      %sanity on type
      assert(isa(dn,'datanames')&&numel(dn)==1,...
        ['Can only handle input ''in'' as scalar of class dataname, not a ',...
        class(dn),' with ',num2str(numel(dn)),' entries.'])
      %parse optional parameters as defined in the metadata
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('dn_reference',dn, @(i) isa(i,'datanames'));
      p=obj.product_get(dn).plot_args(p,varargin{:});
      %if dn_reference is not the default, parse plot arguments from that product
      if ~strcmp(p.UsingDefaults,'dn_reference')
        p=obj.product_get(p.Results.dn_reference).plot_args(p,varargin{:});
      end
      %retrieve the requested data
      d=obj.data_get_scalar(dn);
      %checking data class
      switch class(d)
      case 'gravity'
        h=d.plot(...
          'method',  p.Results.plot_method,...
          'columns', p.Results.plot_columns,...
          'outlier', p.Results.plot_outlier,...
          'zeromean',p.Results.plot_zeromean...
        );
      case 'simpletimeseries'
        h=d.plot(...
          'columns', p.Results.plot_columns,...
          'outlier', p.Results.plot_outlier,...
          'zeromean',p.Results.plot_zeromean...
        );
      otherwise
        if isempty(d)
          str.say('Skip plotting empty data',dn.name)
          h=[];
        else
          error([mfilename,': cannot plot data of class ',class(d),'; implementation needed!'])
        end
      end
      %check units
      if ~isempty(h) && all(p.Results.plot_columns>0)
        h.y_units=d.y_units{p.Results.plot_columns(1)};
        for i=2:numel(p.Results.plot_columns)
          if ~isempty(h.y_units) && ~strcmp(h.y_units,d.y_units(p.Results.plot_columns(i)))
            error([mfilename,':BUG TRAP: y-units are not consistent in all plotted lines.'])
          end
        end
      end
      obj.log('@','out','dn',dn,'start',obj.start,'stop',obj.stop)
    end
    function plot_legend(~,p,h,dn_list,varargin)
      %add legend only if there are multiple lines
      if numel(h)>1
        %add the legend given as input, if there
        if ~isempty(p.Results.plot_legend)
          %some sanity
          str.sizetrap(h,p.Results.plot_legend)
          %propagate
          legend_str=p.Results.plot_legend;
        else
          %get unique parts in datanames 
          prefixes=datanames.unique(dn_list);
          % paranoid sanity
          str.sizetrap(h,prefixes)
          %get how many lines have been plotted
          n=0;
          for j=1:numel(h)
            if ~isempty(h{j}) && isfield(h{j}, 'handle')
              n=n+numel(h{j}.handle);
            end
          end
          %make room for legend strings
          legend_str=cell(1,n);
          %loop over all legend entries
          c=0;
          for j=1:numel(h)
            if ~isempty(h{j})
              for k=1:numel(h{j}.handle)
                %define sufix, add y_mean if lines have zero mean
                if isfield(h{j}, 'y_mean') && p.Results.plot_zeromean
                  suffix=num2str(h{j}.y_mean{k});
                else
                  suffix='';
                end
                %build legend strings
                c=c+1;
                legend_str{c}=str.clean(...
                  strjoin([prefixes{j};{suffix}],' '),...
                  {'title','succ_blanks'}...
                );
              end
            end
          end
        end
        %remove strings as needed
        legend_str=cellfun(@(i) str.clean(i,p.Results.plot_legend_suppress),legend_str,'UniformOutput',false);
        %put the legend in the current plot
        set(legend(legend_str),'fontname','courier')
      end
    end
    function plot_title( ~,p,~,dn_list,varargin)
      %add the title given as input, if there
      if ~isempty(p.Results.plot_title)
        title_str=p.Results.plot_title;
      else
        %get title parts
        title_str=datanames.common(dn_list);
        %suppress some parts, if requested
        title_str=setdiff(title_str,p.Results.plot_title_suppress,'stable');
        %clean up the tile and put it there
        title_str=strjoin([{p.Results.plot_title_prefix};title_str;{p.Results.plot_title_suffix}],' ');
      end
      %put the title in the current plot
      title(str.clean(title_str,'_'))
      %handle keywords
      if strcmp(title_str,'clear')
        title('')
      end
    end
    function plot_ylabel(~,p,h,~,      varargin)
      %add the y-label given as input, if there
      if ~isempty(p.Results.plot_ylabel)
        ylabel_str=p.Results.plot_ylabel;
      elseif ~isfield(h,'y_units')
        %do nothing
        ylabel_str='';
      else
        %get non-empty indexes
        good_idx=find(cellfun(@(i) ~isempty(i) && ~isempty(i.y_units),h));
        %check the labels of all lines are compatible
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
      %put the ylabel in the current plot, if not empty
      if ~isempty(ylabel_str)
        ylabel(ylabel_str)
      end
    end
    function plot_annotate(obj,h,dn,dn_list,varargin)
      %parse mandatory args
      assert(iscell(h),['can only handle input ''h'' as cell, not of class ',class(dn),'.'])
      % paranoid sanity
      str.sizetrap(h,dn_list)
      %parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=obj.product_get(dn).plot_args(p,varargin{:});
      %call annotation methods
      obj.plot_legend(p,h,dn_list,varargin{:})
      obj.plot_title( p,h,dn_list,varargin{:})
      obj.plot_ylabel(p,h,dn_list,varargin{:})
      %enforce plot preferences (using the metadata of the current dataname)
      obj.product_get(dn).enforce_plot(varargin{:})
      %grid, if requested
      if p.Results.plot_grid;grid on;end
    end
    function e=plot_elements(obj,p,product,varargin)
      %build filename sufix
      if isempty(p.Results.plot_file_suffix)
        e.suffix='';
      else
        e.suffix=['.',p.Results.plot_file_suffix];
      end
      %need the source list
      e.sources=product.source_list;
      %need the labels of the columns
      e.col_names=p.Results.plot_column_names;
      if isempty(e.col_names)
        if product.nr_sources>0
          col_names_dn=e.sources{1};
        else
          col_names_dn=product.dataname;
        end
        cnd=obj.data_get_scalar(col_names_dn);
        e.col_names=cnd.labels;
      end
      %gather list of daily data files
      [~,e.startlist,e.stoplist]=product.file('data',varargin{:},'start',obj.start,'stop',obj.stop);
    end
    %% generalized plotting
    %plot a single data entry
    function h=plot(obj,dn,varargin)
      obj.log('@','in','dn',dn,'start',obj.start,'stop',obj.stop)
      %parse mandatory args
      dn_list=obj.data_list(dn);
      % recursive call multiple plots
      if numel(dn_list)>1
        h=cellfun(@(i) obj.plot(i,varargin{:}),dn_list,'UniformOutput',false);
        return
      end
      %save the single entry in dataname
      dn_now=dn_list{1};
      product=obj.product_get(dn_now);
      %parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=product.plot_args(p,varargin{:});
      %if columns are not to be plotted together, need to expand the calls to obj.plot to include each column
      if ~p.Results.plot_columns_together
        %retrieve column names (a.o.)
        e=obj.plot_elements(p,product,varargin{:});
        %make room for handles
        h=cell(1,numel(p.Results.plot_columns));
        %loop over all data columns to plot
        for i=1:numel(p.Results.plot_columns)
          %buid suffixes for the plots of this data column
%           if ~isempty(p.Results.plot_file_suffix)
%             file_suffix=[e.col_names(plot_columns(i)),{p.Results.plot_file_suffix}];
%           else
%             file_suffix=e.col_names(plot_columns(i));
%           end
%           if ~isempty(p.Results.plot_title_suffix)
%             title_suffix=[e.col_names(plot_columns(i)),{p.Results.plot_title_suffix}];
%           else
%             title_suffix=e.col_names(plot_columns(i));
%           end
          h{i}=obj.plot(dn_now,...
            varargin{:},...
            'plot_file_suffix', strjoin([e.col_names(p.Results.plot_columns(i)),{p.Results.plot_file_suffix}], '.'),...
            'plot_title_suffix',strjoin([e.col_names(p.Results.plot_columns(i)),{p.Results.plot_title_suffix}],' '),...
            'plot_columns_together',true,...
            'plot_columns',p.Results.plot_columns(i)...
          );
        end
      else
        %plot filename arguments
        filename_args=[product.file_args('plot'),{...
          'start',obj.start,...
          'stop',obj.stop,...
          'timestamp',true,...
          'remove_part',p.Results.plot_together,...
          'prefix',p.Results.plot_file_prefix,...
          'suffix',p.Results.plot_file_suffix...
        }];
        %plot filename
        filename=dn_now.file(filename_args{:});
        % check if plot is already there
        if isempty(dir(filename))
          figure('visible',p.Results.plot_visible);
          %retrive data names to be plotted here
          %NOTICE: this will cause plots to appear multiple times if they are not saved
          dn_list_to_plot=obj.data_list(dn_now.edit_field_part(p.Results.plot_together,'*'));
          h=cell(size(dn_list_to_plot));
          %loop over all data to plot
          for i=1:numel(dn_list_to_plot)
            h{i}=obj.justplot(dn_list_to_plot{i},varargin{:});
          end
          %annotate plot
          obj.plot_annotate(h,dn_now,dn_list_to_plot,varargin{:})
          %save this plot
          saveas(gcf,filename)
          str.say('Created plot',filename)
        else
          str.say('Skipped plot',filename)
          h=filename;
        end
      end
      obj.log('@','in','dn',dn,'start',obj.start,'stop',obj.stop)
    end
    %plot multiple data entries (more limited than datastorage.plot)
    function h=plot_mult(obj,dn,dn_list,plot_columns,varargin)
      obj.log('@','in','dn',dn,'dn_list',dn_list,'plot_columns',plot_columns,'start',obj.start,'stop',obj.stop)
      %parse mandatory args
      dn =datanames(dn);
      dn_list=datanames.array(dn_list);
      if ~isnumeric(plot_columns)
        error([mfilename,': can only handle input ''plot_columns'' that is numeric, not of class ''',...
          class(plot_columns),'''.'])
      end
      %expand scalar plot_columns into cell array (usually plot_columns is not a scalar)
      if isscalar(plot_columns)
        plot_columns=num2cell(plot_columns*ones(1,numel(dn_list)));
      elseif numel(plot_columns) ~= numel(dn_list)
        tmp=cell(size(dn_list));
        for i=1:numel(tmp)
          tmp{i}=plot_columns;
        end
        plot_columns=tmp;
      end
      %get current product
      product=obj.product_get(dn);
      %parse optional
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the first metadata
      p=product.plot_args(p,varargin{:});
      %create list with plot parameters for dn_now
      dn_now_plot_args_list=product.plot_args_list;
      %make room for outputs
      h=cell(size(dn_list));
      %loop over all datanames
      for i=1:numel(dn_list)
        h{i}=obj.justplot(dn_list{i},...
          varargin{:},...
          dn_now_plot_args_list{:},... %otherwise the justplot method picks the plot_* parameters define in dn_list{i}
          'plot_columns_together',true,...
          'plot_columns',plot_columns{i}... 
        );
      end
      %annotate plot
      obj.plot_annotate(h,dn,dn_list,varargin{:})
      obj.log('@','iout','dn',dn,'start',obj.start,'stop',obj.stop)
    end
    %parses the sources of a product and calls plot_mult (useful as value of the 'method' metadata field)
    function obj=plot_auto(obj,dn_or_prod,varargin)
      obj.log('@','in','dn_or_prod',dn_or_prod,'start',obj.start,'stop',obj.stop)
      %handle both datanames and dataproduct objects
      switch class(dn_or_prod)
      case 'datanames'
        product=obj.product_get(dn_or_prod);
      case 'dataproduct'
        product=dn_or_prod;
      end
      dn=product.dataname;
      %sanity
      assert(obj.isdata_leaf(dn),'Can only handle leaf products.')
      %parse inputs
      p=inputParser;
      p.KeepUnmatched=true;
      %parse optional parameters as defined in the metadata
      p=product.plot_args(p,varargin{:});
      %retrieve plot elements (repetitive processing of parameters)
      e=obj.plot_elements(p,product);
      %loop over all data
      for t=1:numel(e.startlist)
        %loop over all columns to plot
        for c=1:numel(p.Results.plot_columns)
          obj.log('@','iter','column',c,'start',e.startlist(t),'stop',e.stoplist(t))
          %plot filename arguments
          filename_args=[obj.product_get(dn).file_args('plot'),{...
            'start',e.startlist(t),...
            'stop', e.stoplist(t),...
            'timestamp',true,...
            'remove_part','',...
            'prefix',p.Results.plot_file_prefix...
            'suffix',[e.col_names{p.Results.plot_columns(c)},e.suffix]...
          }];
          %plot filename
          filename=dn.file(filename_args{:});
          if isempty(dir(filename))
            %get the data for the current segment
            obj_curr=obj.trim('start',e.startlist(t),'stop',e.stoplist(t),'dn_list',e.sources,'filename',filename);
            %make sure there is data 
            if any(cell2mat(obj_curr.vector_method_tr('all','nr_valid'))>1)
              %plot it
              figure('visible',p.Results.plot_visible);
              h=obj_curr.plot_mult(dn,e.sources,p.Results.plot_columns(c),...
                varargin{:},...
                'plot_title_suffix',strjoin([e.col_names(p.Results.plot_columns(c)),{p.Results.plot_title_suffix}],' ')...
              );
              %check if any data was plotted
              if all(isempty(h))
                %nothing to save
                str.say('Skipped plot',filename,'(no data plotted)')
                close(gfc)
              else
                %save this plot
                saveas(gcf,filename)
                str.say('Created plot',filename)
              end
            else
              str.say('Skipped plot',filename,['(no data in ',dn.name,')'])
            end
          else
            str.say('Skipped plot',filename,'(file already exists)')
          end
        end
      end
      obj.log('@','in','dn_or_prod',dn_or_prod,'start',obj.start,'stop',obj.stop)
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
      %retrieve required operation
      operation=product.mdget('operation');
      obj.log('@','in','operation',operation)
      %operate on all sources
      for i=1:product.nr_sources
        if i==1
          %get first source
          out=obj.data_get_scalar(product.sources(1));
          msg=[product.str,' = ',product.sources(i).str];
        else
          %get the current source
          in=obj.data_get_scalar(product.sources(i));
          msg=[msg,' ',operation,' ',product.sources(i).str]; %#ok<AGROW>
          %enforce common time domain, assume it is defined by the first argument
          in=in.interp(out.t);
          %operate
          out=out.(operation)(in);
        end
        obj.log('@','iter','operation',msg)
      end
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
      sourcep=obj.product_get(product.sources(1));
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
%         sourcep{i}=obj.product_get(product.sources(i));
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
%           if product.ismd_field(partname)
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
%               if sourcep{j}.ismd_field(partname)
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