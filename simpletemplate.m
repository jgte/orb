classdef simpletemplate < simpletimeseries
  %static
  properties(Constant,GetAccess=private)
    %NOTE: edit this if you add a new parameter
    parameter_list={...
      'char',    'char',@ischar;...
      'scalar',   1    ,@num.isscalar;...
      'flag',     true,@(i) islogical(i) && isscalar(i);...
    };
    %These parameter are considered when checking if two data sets are
    %compatible (and only these).
    %NOTE: edit this if you add a new parameter (if relevant)
    compatible_parameter_list={'char'};
  end
  %These parameters should not modify the data in any way; they should
  %only describe the data or the input/output format of it.
  %NOTE: edit this if you add a new parameter (if read/write)
  properties(GetAccess=public,SetAccess=public)
    char % does nothing
    scalar %scales y at init
    flag %sets y to zero at init if true
  end
  %NOTICE: the properties below don't work with varargs' save method because it needs to
  %        access the properties (which is not possible unless both get/set access is public)
  %NOTE: edit this if you add a new parameter (if read only)
  properties(SetAccess=private)
  end
  %NOTE: edit this if you add a new parameter (if private)
  properties(GetAccess=private)
  end
  %calculated only when asked for
  properties(Dependent)
    %add if relevent
  end
  methods(Static)
    function out=parameters(varargin)
      persistent v
      if isempty(v); v=varargs(simpletemplate.parameter_list); end
      out=v.picker(varargin{:});
    end
    %NOTICE: This method (mainly) makes sense when this class (A) is derived from another
    %        class (B) and you want to be able to convert an object of class B into A. 
    %NOTICE: As example, this method illustates how to transmute simpledata and friends.
    function out=transmute(in)
      if isa(in,'simpletemplate')
        %trivial call
        out=in;
      else
        %transmute into this object
        if isprop(in,'t')
          out=simpletemplate(in.t,in.y,in.varargin{:});
        elseif isprop(in,'x')
          out=simpletemplate(in.x,in.y,in.varargin{:});
        else
          error('Cannot find ''t'' or ''x''. Cannot continue.')
        end
      end
    end
    %% general test for the current object
    function test(method,l,w)
      if ~exist('method','var') || isempty(method)
        method='all';
      end
      if ~exist('l','var') || isempty(l)
        l=10;
      end
      if ~exist('w','var') || isempty(w)
        w=3;
      end

      %get common parameters
      args=simpledata.test_parameters('args',l,w);
      now=juliandate(datetime('now'),'modifiedjuliandate');
      t=datetime(now,           'convertfrom','modifiedjuliandate'):...
        datetime(now+round(l)-1,'convertfrom','modifiedjuliandate');
      %init object
      a=simpletemplate(...
          t,...
          simpledata.test_parameters('y_all',l,w),...
          'mask',simpledata.test_parameters('mask',l,w),...
          'flag',false,...
          args{:}...
        );

      switch method
        case 'all'
          for i={'print'}
            simpletemplate.test(i{1},l);
          end
        case 'print'
          a.print
      end
    end
  end
  methods
    %% constructor
    function obj=simpletemplate(t,y,varargin)
      % input parsing
      p=machinery.inputParser;
      p.addRequired( 't' ); %this can be char, double or datetime
      p.addRequired( 'y', @(i) simpledata.valid_y(i));
      %create argument object, declare and parse parameters, save them to obj
      %NOTICE: the parameter (i.e., char, scalar and flag) are parsed by varargs and saved to the variable v)
      [v,p]=varargs.wrap('parser',p,'sources',{simpletemplate.parameters('obj')},'mandatory',{t,y},varargin{:});
      %do some operations on the basis of the parameters
      if p.Results.flag
        y(:)=0;
      else
        y=y*p.Results.scalar;
      end
      % call superclass
      obj=obj@simpletimeseries(t,y,varargin{:});
      % save the arguments v into this object
      % NOTICE: it is here that the parameters (i.e., char, scalar and flag) are saved to obj
      obj=v.save(obj,{'t','y'});
    end
    function obj=assign(obj,y,varargin)
      %pass it upstream
      obj=assign@simpletimeseries(obj,y,varargin{:});
      %update internal 
      if any(y(:)~=0)
        obj.flag=false;
        obj.scalar=1;
      else
        obj.flag=true;
      end
    end
    function obj=copy_metadata(obj,obj_in,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      obj=copy_metadata@simpletimeseries(obj,obj_in,[simpletemplate.parameters('list');more_parameters(:)]);
    end
    function out=metadata(obj,more_parameters)
      if ~exist('more_parameters','var')
        more_parameters={};
      end
      %call superclass
      out=metadata@simpletimeseries(obj,[simpletemplate.parameters('list');more_parameters(:)]);
    end
    %the varargin method can be called directly
    %% info methods
    function print(obj,tab)
     if ~exist('tab','var') || isempty(tab)
        tab=12;
      end
      %parameters
      relevant_parameters={'char','scalar','flag'};
      for i=1:numel(relevant_parameters)
        obj.disp_field(relevant_parameters{i},tab);
      end
      %print superclass
      print@simpletimeseries(obj,tab)
    end
  end
end