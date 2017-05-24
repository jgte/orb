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
  %read only
  properties
    debug
  end
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
      % parameter names
      pn=datastorage.parameters;
      p=inputParser;
      p.KeepUnmatched=true;
      % declare parameters
      for i=1:numel(pn)
        %declare parameters
        p.addParameter(pn{i},datastorage.parameter_list.(pn{i}).default,datastorage.parameter_list.(pn{i}).validation)
      end
      % parse it
      p.parse(varargin{:});
      % reset data type list
      obj=obj.data_init;
      % save parameters with defaults first, they may be needed below
      for i=1:numel(pn)
        obj.(pn{i})=collapse(p.Results.(pn{i}));
        obj.log('@','iter',pn{i},obj.(pn{i}),['p.Results.(',pn{i},')'],p.Results.(pn{i}))
      end
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
      %get the stack
      stack=dbstack(1);
      %make room for message
      msg=cell((nargin-1)/2+1,1);
      %first item is the calling function
      msg{1}=stack.name;
      %loop over all arguments
      for i=2:numel(msg)
        idx0=2*(i-1)-1;
        idx1=2*(i-1);
        %don't show empty arguments
        if isempty(varargin{idx1})
          continue;
        end
        msg{i}=[varargin{idx0},':',str.show(varargin{idx1})];
      end
      %clean up empty cells
      msg=msg(cellfun(@(i)(~isempty(i)),msg));
      %show the debug message
      disp(strjoin(msg,', '))
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
    function out=tab(obj)
      dn_list=obj.data_list('all');
      out=max(cellfun(@(i) max(cellfun(@numel,i.global_field_path)),dn_list));
    end
    function peek(obj,dn)
      if ~exist('dn','var') || isempty(dn)
        dn='all';
      end
      tab=obj.tab;
      cellfun(@(i) disp(str.tablify(tab,i.global_field_path)),obj.data_list(dn))
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
            out{i}=obj_list{i}.(method)(varargin{:});
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
      validv1=cellfun(@(i)(~isempty(i)),values1);
      values2=obj2.data_get(dn);
      validv2=cellfun(@(i)(~isempty(i)),values2);
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
      %retrieve the requested metadata
      md=obj.par.(dn.name_clean);
      %some sanity
      assert(~isempty(md),['cannot find product ',dn.name,'. ',...
          'Have you called obj.init(''',dn.name,''')? ',...
          'Maybe the requested product is not sourced in any of the products already loaded.'])
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
      file_list=product.file('data','discover',true,'start',obj.start,'stop',obj.stop,varargin{:});
      % check if any file is missing
      if any(~product.isfile('data','start',obj.start,'stop',obj.stop,varargin{:}))
        % debug output
        success=false;
        obj.log('@','out','no files to load for',product)
        return
      end
      obj.log('nr of files to load',numel(file_list))
      %loop over all files
      for f=1:numel(file_list)
        % debug output
        obj.log(['loading file ',num2str(f)],file_list{f})
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
      obj.log('@','in','dataname',dn,'start',obj.start,'stop',obj.stop)
      % load the metadata for this product
      obj=obj.product_set(dn,varargin{:});
      %retrieve product info
      product=obj.product_get(dn);
      % make sure all data sources are loaded
      obj=obj.init_sources(product,varargin{:});
%       %save start/stop (can be changed in the init method, because some data is stored globally)
%       start_here=obj.start;
%       stop_here =obj.stop;
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
%       % enforce start/stop
%       if isempty(start_here)
%         obj.starti=obj.start;
%       else
%         obj.start=start_here;
%       end
%       if isempty(stop_here)
%         obj.stopi=obj.stop;
%       else
%         obj.stop=stop_here;
%       end
      obj.log('@','out','dataname',dn,'start',obj.start,'stop',obj.stop)
      %user feedback
      if ~obj.debug
        disp(['initialized ',product.str])
      end
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
          %expand this product list
          product_list=product_list.field_path_expand(cells.flatten(source_inventory{1}));
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
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired('product',           @(i) isa(i,'dataproduct'));
      p.addParameter('start', obj.start, @(i) isdatetime(i)  &&  isscalar(i));
      p.addParameter('stop',  obj.stop,  @(i) isdatetime(i)  &&  isscalar(i));
      % parse it
      p.parse(product,varargin{:});
      % sanity
      if isempty(p.Results.start) || isempty(p.Results.stop)
        error([mfilename,': need ''start'' and ''stop'' parameters (or non-empty obj.start and obj.stop).'])
      end
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
        tmp=nrtdm([nrtdm_sats{s},'_',nrtdm_product],p.Results.start,p.Results.stop,'data_dir',indir);
        %get the data
        obj=obj.sat_set(product.dataname.type,product.dataname.level,product.dataname.field,sats{s},...
          tmp.ts...
        );
      end
      %make sure start/stop options are honoured (if non-empty)
      if ~isempty(obj)
        obj.start=p.Results.start;
        obj.stop= p.Results.stop;
      end
    end
    %% plot utils
    function h=justplot(obj,dn,varargin)
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
      d=obj.data_get(dn);
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
          disp(['Skip plotting empty data ',dn.name])
          h=[];
        else
          error([mfilename,': cannot plot data of class ',class(d),'; implementation needed!'])
        end
      end
      %check units
      if ~isempty(h) && all(p.Results.plot_columns>0)
        h.y_units=d.y_units{p.Results.plot_columns(1)};
        for i=2:numel(p.Results.plot_columns)
          if ~strcmp(h.y_units,d.y_units(p.Results.plot_columns(i)))
            error([mfilename,':BUG TRAP: y-units are not consistent in all plotted lines.'])
          end
        end
      end
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
        title_str=strjoin([title_str;{p.Results.plot_title_suffix}],' ');
      end
      %put the title in the current plot
      title(str.clean(title_str,'_'))
      %handle keywords
      if strcmp(title_str,'clear')
        title('')
      end
    end
    function plot_ylabel(~,p,h,~,              varargin)
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
        %retrieve plot columns indexes, these are naturally numeric arrays when read from the metadata files. keep it that way
        plot_columns=p.Results.plot_columns;
        %retrieve column names (a.o.)
        e=obj.plot_elements(p,product,varargin{:});
        %make room for handles
        h=cell(1,numel(plot_columns));
        %loop over all data columns to plot
        for i=1:numel(plot_columns)
          %buid suffixes for the plots of this data column
          if ~isempty(p.Results.plot_file_suffix)
            file_suffix=[e.col_names(plot_columns(i)),{p.Results.plot_file_suffix}];
          else
            file_suffix=e.col_names(plot_columns(i));
          end
          if ~isempty(p.Results.plot_title_suffix)
            title_suffix=[e.col_names(plot_columns(i)),{p.Results.plot_title_suffix}];
          else
            title_suffix=e.col_names(plot_columns(i));
          end
          h{i}=obj.plot(dn_now,...
            varargin{:},...
            'plot_file_suffix', strjoin(file_suffix, '.'),...
            'plot_title_suffix',strjoin(title_suffix,' '),...
            'plot_columns_together',true,...
            'plot_columns',plot_columns(i)...
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
          'suffix',p.Results.plot_file_suffix,...
          'full_path',p.Results.plot_file_full_path,...
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
        else
          h=filename;
        end
      end
    end
    %plot multiple data entries (more limited than datastorage.plot)
    function h=plot_mult(obj,dn,dn_list,plot_columns,varargin)
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
    end
    %parses the sources of a product and calls plot_mult (useful as value of the 'method' metadata field)
    function obj=plot_auto(obj,dn_or_prod,varargin)
      obj.log('@','in','dn_or_prod',dn_or_prod)
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
          obj.log('@','iter','column',c,'start',e.startlist(t),'stop',e.stoptlist(t))
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
            obj_curr=obj.trim('start',e.startlist(t),'stop',e.stoplist(t),'dn_list',e.sources);
            %make sure there is data 
            if any(cell2mat(obj_curr.vector_method_tr('all','nr_valid'))>1)
              %plot it
              figure('visible',p.Results.plot_visible);
              h=obj_curr.plot_mult(dn,e.sources,p.Results.plot_columns(c),...
                varargin{:},...
                'plot_title_suffix',e.col_names{p.Results.plot_columns(c)});
              %check if any data was plotted
              if all(isempty(h))
                %nothing to save
                obj.log('@','iter','skipped',dn.name,'file not created',filename)
                close(gfc)
              else
                %save this plot
                saveas(gcf,filename)
                obj.log('@','iter','plot saved',filename)
              end
              % user feedback
              if strcmp(p.Results.plot_visible,'off')
                disp(['plot_auto: plotted ',dn.name,' to file ',filename])
              end
            else
              disp(['plot_auto: not enough data to plot ',dn.name,' to file ',filename,' (skipped)'])
            end
          else
            obj.log('@','iter','skipped',dn.name,'file already exists',filename)
          end
        end
      end
      obj.log('@','out')
    end
    %% usefull stuff
    function [sourcep,df]=dataflow(obj,product)
      %retrive source metadata
      sourcep=cell(1,product.nr_sources);
      for i=1:product.nr_sources
        sourcep{i}=obj.product_get(product.sources(i));
      end
      %retrieve dataflow structure
      df=product.mdget('dataflow');
      %loop over the types/levels/fields metadata fields
      for i=1:numel(obj.parts)
        %easier names
        partname=[obj.parts{i},'s'];
        %check if this part is there
        if isfield(df,partname)
          %retrive metadata value
          partvalues=df.(partname);
          %multiple rows means this product spans multiple types/levels/fields
          if iscell(partvalues{1})
            %need to unwrap nested cells case there are multiple rows
            df.(partname)=cell(numel(partvalues),numel(sourcep)+1);
            for j=1:numel(partvalues)
              %sanity on types/levels/fields: number of columns must be equal to number of sources
              assert(numel(partvalues{j}) == numel(sourcep)+1,[mfilename,': ',...
                'ilegal ''',partname,''' entry of the ''dataflow'' metadata field, ',...
                'it must have the same number of columns as the number of product sources +1 (',num2str(numel(sourcep)+1),'), ',...
                'not ',num2str(numel(partvalues{j})),'.'])
              df.(partname)(j,:)=partvalues{j};
            end
          else
            %variable partvalues already contains types/levels/fields
            assert(numel(partvalues) == numel(sourcep)+1,[mfilename,': ',...
              'ilegal ''',obj.parts{i},'s'' metadata field, ',...
              'it must have the same number of columns as the number of product sources (',num2str(numel(sourcep)+1),'), ',...
              'not ',num2str(numel(partvalues)),'.'])
          end
        else
          %if the dataflow structure is missing a part, then patch from this product or sources
          if product.ismd_field(partname)
            df.(partname)=product.mdget(partname,'always_cell_array',true);
            %this avoids having to define repetitive (e.g.) sats values in the 'dataflow' metadata field,
            %the scalar 'sats' metadata fiels is enough (doesn't work for 'types')
            if i>1 && numel(df.(partname))==1 && numel(df.([obj.parts{1},'s']))>1
              partvalues=cell(size(df.([obj.parts{1},'s'])));
              partvalues(:)=df.(partname);
              df.(partname)=partvalues;
            end
          else
            found_part=false;
            %search in all sources for an explicit list of this partname
            for j=1:product.nr_sources
              if sourcep{j}.ismd_field(partname)
                partvalues=sourcep{j}.mdget(partname);
                df.(partname)=cell(numel(partvalues),numel(sourcep)+1);
                for k=1:numel(sourcep)+1
                  df.(partname)(:,k)=partvalues(:);
                end
                found_part=true;
              end
            end
            if ~found_part
              %if nothing found, search for values of this partname in the data of the sources
              partvalues=obj.data_list(sourcep{1}.name);
              %check if all sources share the same values of this part
              for j=2:product.nr_sources
                %cell array value comparisson:
                %http://stackoverflow.com/questions/3231580/matlab-comparison-of-cell-arrays-of-string
                assert(cells.isequal(obj.data_list(sourcep{j}.name),partvalues),...
                  ['The ',partname,' names differ in products ',sourcep{1}.name,' and ',...
                  sourcep{j}.name,'.'])
              end
              %propagate this part values
              df.(partname)=cell(numel(partvalues),numel(sourcep)+1);
              for j=1:numel(partvalues)
                df.(partname)(j,:)=partvalues(j);
              end
            end
          end
        end
      end
      %return scalar source product if only one
      if numel(sourcep)==1
        sourcep=sourcep{1};
      end
    end
    %% generalized operations
    function obj=stats(obj,product,varargin)
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
    end
    function obj=corr(obj,product,varargin)
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
    function obj=arithmetic(obj,product,varargin)
      obj.log('@','in','product',product,'start',obj.start,'stop',obj.stop)
      %parse mandatory arguments
      p=inputParser;
      p.addRequired('product', @(i) isa(i,'dataproduct'));
      p.parse(product);
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
    %% costumized operations
    function out=data_op(obj,opname,datatype,level,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addParameter('outlier',  0, @(i) isfinite(i));
      % parse it
      p.parse(varargin{:});
      switch lower(datatype)

      case 'acc'
        switch lower(opname)
        case 'load-stats'
          %compute residuals relative to the available models
          models=obj.field_list('acc','mod');
          for m=1:numel(models)
            filenames=simpletimeseries.unwrap_datafiles(...
              obj.id('data','accres_csr',level,models{m},'DATE_PLACE_HOLDER','stats'),...
              'start',obj.start,...
              'stop' ,obj.stop,...
              'only_existing',false...
            );
            for i=1:numel(filenames)
              if isempty(dir(filenames{i}))
                for s=1:numel(grace.sats)
                  %get the l1b acc data
                  l1bcsr=obj.sat_get('acc','l1b','csr',grace.sats{s});
                  %get the calibration model
                  calmod=obj.sat_get('acc','calmod',level,grace.sats{s});
                  %get the modeled acc
                  accmod=obj.sat_get('acc','mod',models{m},grace.sats{s});
                  if any([...
                    ~isa(l1bcsr,'simpletimeseries'),...
                    ~isa(calmod,'simpletimeseries'),...
                    ~isa(accmod,'simpletimeseries')...
                  ])
                    %patch nan residuals
                    accres=simpletimeseries(...
                      [obj.start;obj.stop],...
                      nan(2,numel(obj.par.acc.data_col_name))...
                    );
                  else
                    %calibrate the accelerometer
                    acccal=l1bcsr+calmod;
                    %compute residual
                    accres=acccal-accmod;
                  end
                  %save descriptor
                  accres.descriptor=str.clean(filenames{i},{'file','.'});
                  %compute stats per period
                  s.(grace.sats{s})=accres.stats(...
                    'period',days(1),...
                    'overlap',seconds(0),...
                    'outlier',p.Results.outlier,...
                    'struct_fields',obj.par.modes.stats...
                  );
                end
                save(filenames{i},'s')
              else
                load(filenames{i})
              end
            end
            %save this data to output var
            out.(models{m})=s;
          end
        case 'load-stats2'
          error([mfilename,': implementation needed.'])
        otherwise
          error([mfilename,': cannot handle operation ''',opname,''' for datatype ''',datatype,'''.'])
        end
      otherwise
        error([mfilename,': cannot handle datatype ''',datatype,'''.'])
      end
    end
    %% costumized plotting routine
    function oldplot(obj,plotname,datatype,level,varargin)
      switch lower(datatype)
      case 'calpar_csr'
        switch lower(plotname)
        case ''
          % loop over CSR cal pars
          fields=obj.field_list(datatype,level);
          for j=1:numel(fields)
            filename=obj.id('plot',datatype,level,fields{j},'');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %retrive data
              sat=struct(...
                'A',obj.sat_get(datatype,level,fields{j},'A'),...
                'B',obj.sat_get(datatype,level,fields{j},'B')...
              );
              %plot data
              sat.A.plot('columns',obj.par.calpar_csr.data_col_idx)
              sat.B.plot('columns',obj.par.calpar_csr.data_col_idx)
              %legend
              legend({'GRACE-A','GRACE-B'});
              obj.enforce_plot_par('autoscale',true,'automean',true,'outlier',3)
              grid on
              ylabel(sat.A.y_units{1}) %could also be sat.B.y_units{1}
              title([fields{j},' ',level])
              saveas(gcf,filename)
            end
          end
        case 'stats'
          s=obj.data_op('load-stats',datatype,level,varargin{:});
          fields=fieldnames(obj.par.calpar_csr.fields);
          %plot number of data
          obj.par.plot.prefix='length.bars';
          for i=1:numel(fields)
            filename=obj.id('plot',datatype,level,fields{i},'');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              b=bar(...
                mean([...
                  datenum(s.(fields{i}).A.length.t),...
                  datenum(s.(fields{i}).B.length.t)...
                ],2),[...
                  s.(fields{i}).A.length.scale(0.5).y(:,1),...
                  s.(fields{i}).B.length.scale(0.5).y(:,1)...
                ]...
              );
              b(1).EdgeColor = 'red';
              b(1).FaceColor = 'red';
              b(2).EdgeColor = 'black';
              b(2).FaceColor = 'black';
              %plot tweaking
              obj.enforce_plot_par
              legend('GRACE-A','GRACE-B')
              datetick('x','yyyy-mm')
              xlabel('time')
              grid on
              title(['Monthly count of ',fields{i},' ',level])
              saveas(gcf,filename)
            end
          end
          %plot more stats
          for mode=obj.par.modes.stats;
            %plot the requested statistics of the data
            obj.par.plot.prefix=mode{1};
            for i=1:numel(fields)
              filename=obj.id('plot',datatype,level,fields{i},'');
              if isempty(dir(filename))
                switch mode{1}
                case 'length'
                  continue
                  mask=true(s.(fields{i}).A.(mode{1}).length,1); %#ok<UNRCH>
                  scl=0.5;
                otherwise
                  mask=s.(fields{i}).A.length.y(:,1)>10;
                  scl=1;
                end
                sat=struct(...
                  'A',s.(fields{i}).A.(mode{1}).mask_and(mask),...
                  'B',s.(fields{i}).B.(mode{1}).mask_and(mask)...
                );
                figure('visible',obj.par.plot.visible);
                sat.A.scale(scl).plot('columns',obj.par.calpar_csr.data_col_idx);
                sat.B.scale(scl).plot('columns',obj.par.calpar_csr.data_col_idx);
                %plot tweaking
                obj.enforce_plot_par('autoscale',true,'automean',true,'outlier',1)
                legend('GRACE-A','GRACE-B')
                ylabel(['[',sat.A.y_units{obj.par.calpar_csr.data_col_idx},']'])
                grid on
                title(['monthly ',mode{1},' of ',fields{i},' ',level])
                saveas(gcf,filename)
              end
            end
          end
        case 'corr'
          s=obj.data_op('load-stats2',datatype,level,varargin{:});
          fields=fieldnames(obj.par.calpar_csr.fields);
          % loop over all required modes
          for mode=obj.par.modes.stats2;
            %plot the requested statistics of the data
            obj.par.plot.prefix=mode{1};
            % loop over CSR cal pars
            for i=1:numel(fields)
              filename=obj.id('plot',datatype,level,fields{i},'');
              if isempty(dir(filename))
                switch mode{1}
                case 'length'
                  continue
                  mask=true(s.(fields{i}).(mode{1}).length,1); %#ok<UNRCH>
                  ylimits=[-Inf,Inf];
                  scl=0.5;
                otherwise
                  mask=s.(fields{i}).length.y(:,1)>=10;
                  ylimits=[-1,1];
                  scl=1;
                end    
                figure('visible',obj.par.plot.visible);
                s.(fields{i}).(mode{1}).mask_and(mask).scale(scl).plot(...
                  'columns',obj.par.calpar_csr.data_col_idx...
                );
                obj.enforce_plot_par('ylimits',ylimits)
                datetick('x','yyyy-mm')
                xlabel('time')
                ylabel('[ ]')
                grid on
                title(['monthly ',mode{1},' of ',fields{i},' ',level])
                saveas(gcf,filename)
              end
            end
          end
        case 'overview'
          fields={...
            'AC0X',...
            'AC0Y1',...
            'AC0Z',...
            'AC0XD',...
            'AC0YD1',...
            'AC0ZD',...
            'AC0XQ',...
            'AC0YQ1',...
            'AC0ZQ'...
          };
          %single-sat stats
          s=obj.data_op('load-stats',datatype,level,varargin{:});
          stats_modes=obj.par.modes.stats;
          for j=1:numel(stats_modes);
            %skip length, not so interesting in this mode
            if strcmp(stats_modes{j},'length')
              continue
            end
            obj.par.plot.prefix=stats_modes{j};
            filename=obj.id('plot',datatype,level,'','');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %gathering data
              y=nan(numel(fields),2);
              for i=1:numel(fields)
                if ~isfield(s,fields{i})
                  continue
                end
                y(i,1)=s.(fields{i}).A.(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.calpar_csr.data_col_idx);
                y(i,2)=s.(fields{i}).B.(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.calpar_csr.data_col_idx);
              end
              pos_idx=(y>0);
              %plot negative parameters downwards
              y_now=abs(y); y_now(pos_idx)=nan;
              b{1}=bar(log10(y_now)); hold on
              %plot positive parameters upwards
              y_now=abs(y); y_now(~pos_idx)=nan;
              b{2}=bar(-log10(y_now)); hold on
              %annotating
              obj.enforce_plot_par
              for k=1:2
                b{k}(1).EdgeColor = 'red';
                b{k}(1).FaceColor = 'red';
                b{k}(2).EdgeColor = 'black';
                b{k}(2).FaceColor = 'black';
              end
              if all(pos_idx(~isnan(y)))
                axis([0 numel(fields)+1 0 10])
              else
                axis([0 numel(fields)+1 -10 10])
              end
              set(gca,'XTick',1:numel(fields))
              set(gca,'XTickLabel',fields)
              set(gca,'XTickLabelRotation',45)
              yticks=get(gca,'YTickLabel');
              for k=1:numel(yticks)
                ytickval=str2double(yticks{k});
                if ytickval==0
                  yticks{k}='1';
                elseif ytickval>0
                  yticks{k}=['y.10^{-',yticks{k},'}'];
                elseif ytickval<0
                  yticks{k}=['-y.10^{',yticks{k},'}'];
                end
              end
              set(gca,'YTickLabel',yticks);
              legend('GRACE-A','GRACE-B')
              grid on
              title([...
                stats_modes{j},' of ',level,' for ',...
                datestr(obj.par.plot.xlimits(1)),' to ',...
                datestr(obj.par.plot.xlimits(2))...
              ]);
              saveas(gcf,filename)
            end
          end
          %all-sat stats
          s=obj.data_op('load-stats2',datatype,level,varargin{:});
          stats_modes=obj.par.modes.stats2;
          for j=1:numel(stats_modes);
            %skip length, not so interesting in this mode
            if strcmp(stats_modes{j},'length')
              continue
            end
            obj.par.plot.prefix=stats_modes{j};
            filename=obj.id('plot',datatype,level,'','');
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              %gathering data
              y=nan(numel(fields),1);
              for i=1:numel(fields)
                if ~isfield(s,fields{i})
                  continue
                end
                y(i,1)=s.(fields{i}).(stats_modes{j}).trim(...
                  obj.par.plot.xlimits(1),...
                  obj.par.plot.xlimits(2)...
                ).stats(...
                  'period',days(inf),...
                  'outlier',1,...
                  'nsigma',3 ...
                ).mean(obj.par.(datatype).data_col_idx);
              end
              bar(y,'k')
              %annotating
              obj.enforce_plot_par
              set(gca,'XTick',1:numel(fields))
              set(gca,'XTickLabel',fields)
              set(gca,'XTickLabelRotation',45)
              axis([0 numel(fields)+1 -1 1])
              grid on
              title([...
                stats_modes{j},' of ',level,' for ',...
                datestr(obj.par.plot.xlimits(1)),' to ',...
                datestr(obj.par.plot.xlimits(2))...
              ]);
              saveas(gcf,filename)
            end
          end  
%           case 'calpar_csr-tables'
%             fields=fieldnames(obj.par.calpar_csr.fields);
%             tab=12;
%             %single-sat stats
%             s=obj.data_op('load-stats',datatype,level,varargin{:});
%             stats_modes=obj.par.modes.stats;
%             for s=1:numel(grace.sats)
%               disp(['GRACE-',upper(grace.sats{s})])
%               for i=0:numel(fields)
%                 out=cell(1,numel(stats_modes)+1);
%                 if i==0
%                   calpar_clean='cal par';
%                 else
%                   calpar_clean=str.clean(fields{i},'_');
%                 end
%                 out{1}=str.tabbed(calpar_clean,tab);
%                 for j=1:numel(stats_modes);
%                   if i==0
%                     out{j+1}=str.tabbed(stats_modes{j},tab,true);
%                   else
%                     d=s{i}.(grace.sats{s}).(stats_modes{j}).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
%                     out{j+1}=str.tabbed(num2str(d(obj.par.calpar_csr.data_col_idx),'% .3g'),tab,true);
%                   end
%                 end
%                 disp(strjoin(out,' '))
%               end
%             end
%             %all-sat stats
%             s=obj.data_op('load-stats2',datatype,level,varargin{:});
%             stats_modes=obj.par.modes.stats2;
%             for i=0:numel(fields)
%               out=cell(1,numel(stats_modes)+1);
%               if i==0
%                 calpar_clean='cal par';
%               else
%                 calpar_clean=str.clean(fields{i},'_');
%               end
%               out{1}=str.tabbed(calpar_clean,tab);
%               for j=1:numel(stats_modes);
%                 if i==0
%                   out{j+1}=str.tabbed(stats_modes{j},tab,true);
%                 else
%                   d=s{i}.(stats_modes{j})(:,1).stats('period',days(inf),'outlier',true,'nsigma',3).mean;
%                   out{j+1}=str.tabbed(num2str(mean(d(obj.par.calpar_csr.data_col_idx)),'% .3f'),tab,true);
%                 end
%               end
%               disp(strjoin(out,' '))
%             end
        otherwise
          error([mfilename,': cannot handle plotname ''',plotname,''' for datatype ''',datatype,'''.'])
        end
      case 'acc'
        switch lower(plotname)
        case {'x','y','z'}
          for s=1:numel(grace.sats)
            data_col_idx=find(~cellfun(@isempty,strfind(obj.par.acc.data_col_name,upper(plotname))));
            data_col_name=obj.par.acc.data_col_name{data_col_idx}; 
            filename=obj.id('plot',datatype,level,...
              [data_col_name,'.',datestr(obj.start,'yyyymmddTHHMMSS'),'-',datestr(obj.stop,'yyyymmddTHHMMSS')],...
            grace.sats{s});
            if isempty(dir(filename))
              figure('visible',obj.par.plot.visible);
              ts={...
                obj.data.(datatype).l1b.csr.(       grace.sats{s}),...
                obj.data.(datatype).calmod.(level).(grace.sats{s}),...
                obj.data.(datatype).mod.csr.(       grace.sats{s}),...
                obj.data.(datatype).mod.nrtdm.(     grace.sats{s})...
              };
              idx=0;
              legend_str=cell(0);
              idx=idx+1;
              if isa(ts{1},'simpletimeseries')
                h=ts{1}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['L1B         ',num2str(h.y_mean{1})];
                y_units=ts{1}.y_units{data_col_idx};
              end
              idx=idx+1;
              if isa(ts{3},'simpletimeseries')
                calib=ts{1}+ts{2};
                h=calib.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['Calibrated  ',num2str(h.y_mean{1})];
              end
              idx=idx+1;
              if isa(ts{2},'simpletimeseries')
                h=ts{2}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['Calib. Mod. ',num2str(h.y_mean{1})];
                y_units=ts{2}.y_units{data_col_idx};
              end
              idx=idx+1;
              if isa(ts{3},'simpletimeseries')
                h=ts{3}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['CSR Mod.    ',num2str(h.y_mean{1})];
                y_units=ts{3}.y_units{data_col_idx};
              end
              idx=idx+1;
              if isa(ts{4},'simpletimeseries')
                h=ts{4}.plot('columns',data_col_idx,'zeromean',true);
                legend_str{idx}=['NRLMSISE-00 ',num2str(h.y_mean{1})];
                y_units=ts{4}.y_units{data_col_idx};
              end
              if all(cellfun(@isempty,legend_str))
                disp(['No data to plot ',filename])
                return
              end
              %plot tweaking
              obj.enforce_plot_par('autoscale',true,'automean',true)
              %add legend
              legend_idx=cellfun(@isempty,legend_str);
              lh=legend(legend_str{~legend_idx});
              set(lh,'fontname','courier')
              ylabel(['[',y_units,']'])
              grid on
              title([level,' ACC ',lower(plotname),'-componend GRACE-',grace.sats{s}])
              saveas(gcf,filename)
            end            
          end
        case 'stats'
          stats=obj.data_op('load-stats',datatype,level,varargin{:});
          models=fieldnames(stats.(level));
          stat_modes=obj.par.modes.stats;
          for s=1:numel(grace.sats)
            for sm=1:numel(stat_modes)
              for i=1:numel(obj.par.(datatype).data_col_idx)
                data_col_idx=obj.par.(datatype).data_col_idx(i);
                data_col_name=obj.par.acc.data_col_names{data_col_idx};
                filename=obj.id('plot',datatype,level,[stat_modes{sm},'-',data_col_name],grace.sats{s});
                if isempty(dir(filename))
                  figure('visible',obj.par.plot.visible);
                  legend_str=cell(1,numel(models));
                  for m=1:numel(models)
                    %get current stats timeseries
                    ts=stats.(level).(models{m}).(grace.sats{s}).(stat_modes{sm});
                    %only show stats when there is enough data
                    switch stat_modes{1}
                    case 'length'
                      mask=true(ts.length,1);
                      scl=0.5;
                    otherwise
                      mask=stats.(level).(models{m}).(grace.sats{s}).length.y(:,1)>10;
                      scl=1;
                    end
                    %plot it
                    h=ts.mask_and(mask).scale(scl).plot('columns',data_col_idx);
                    legend_str{m}=[models{m},' \mu=',h.y_mean{data_col_idx}];
                  end
                  %plot tweaking
                  obj.enforce_plot_par('autoscale',true,'automean',true)
                  legend(legend_str)
                  ylabel(['[',ts.y_units{data_col_idx},']'])
                  grid on
                  title(['daily ',stat_modes{1},' of ',level,'-calib. acc. residuals along ',data_col_name,...
                    'for GRACE-',grace.sats{s}])
                  saveas(gcf,filename)
                end
              end
            end
          end
        otherwise
          error([mfilename,': cannot handle plotname ''',plotname,''' for datatype ''',datatype,'''.'])
        end
      otherwise
        error([mfilename,': unknown data of type ''',datatype,'''.'])
      end
    end
  end
end

function out=collapse(in)
  if ~isscalar(in) && ~iscell(in)
    %vectors are always lines (easier to handle strings)
    out=transpose(in(:));
  else
    out=in;
  end
end
