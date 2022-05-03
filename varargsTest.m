classdef varargsTest < matlab.unittest.TestCase
  methods(Static)
    function go
       t=varargsTest;
       disp(table(t.run))
    end
  end
  properties(Constant)
    par=[struct(...
      'name' ,'char',...
      'value','some string',...
      'valid',@ischar...
    ),struct(...
      'name' ,'numeric',...
      'value',[1,2],...
      'valid',@isnumeric...
    ),struct(...
      'name' ,'cell_array',...
      'value',{{3,4}},...
      'valid',@iscell...
    ),struct(...
      'name' ,'cell_str',...
      'value',{{'a','b','c'}},...
      'valid',@iscellstr...
    )];
    cell_matrix=permute(struct2cell(varargsTest.par),[3,1,2]);
  	cell_vector=reshape(transpose(varargsTest.cell_matrix(:,1:2)),1,numel(varargsTest.par)*2);
    par_names=arrayfun(@(i) i.name,varargsTest.par,'UniformOutput',false);
    S=arrayfun(@(i) varargs.initS(i.name,i.value,'validation',i.valid),varargsTest.par);
    obj = varargs(varargsTest.cell_matrix);
  end
  methods(Test)
    function isvararginTest(testCase)
      testCase.verifyTrue(varargs.isvarargs(varargsTest.cell_vector));
    end
    function iscellTest(testCase)
      testCase.verifyTrue(varargs.isvarargs(varargsTest.cell_matrix));
    end
    function isSTest(testCase)
      testCase.verifyTrue(varargs.isvarargs(varargsTest.S));
    end
    function isobjTest(testCase)
      testCase.verifyTrue(varargs.isvarargs(varargsTest.obj));
    end
    function constructorStructTest(testCase)
      testCase.verifyTrue(varargsTest.obj==varargs(varargsTest.S));
    end
    function constructorCellTest(testCase)
      testCase.verifyTrue(varargsTest.obj==varargs(varargsTest.cell_matrix));
    end
    function wrapTest(testCase)
      %create dummy required inputs
      required1=10;
      required2='string';
      optional1={'str',0};
      new_char='another string';
      new_numeric=[3 4];
      ignored1='another string';
      ignored2=NaN;
      %implement traditional parses
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'required1' ,         @isnumeric);
      p.addRequired( 'required2' ,         @ischar);
      p.addParameter('optional1' ,   {0} , @iscell);
      %create argument object, declare and parse parameters, save them to sink
      [v,p,sink]=varargs.wrap(...
        'sinks'    ,{struct('numeric',0,'cell_array',{0})},... %fake object, only receives parameters defined in 'sources'
        'parser'   , p,...
        'sources'  ,{varargsTest.cell_matrix},...
        'mandatory',{required1,required2},...
        'optional1', optional1,...
        'char'     , new_char,...
        'numeric'  , new_numeric,...
        'ignored1' , ignored1,...
        'ignored2' , ignored2...
      );
      testCase.verifyTrue(cells.isequal(v.rest,{'ignored1',ignored1,'ignored2',ignored2}),'v.ignored')
      testCase.verifyTrue(cells.isequal(v.required1,required1),'v.required1')
      testCase.verifyTrue(cells.isequal(v.required2,required2),'v.required2')
      testCase.verifyTrue(cells.isequal(v.char,new_char),      'v.new_char')
      testCase.verifyTrue(cells.isequal(v.numeric,new_numeric),'v.new_numeric')
      testCase.verifyTrue(cells.isequal(v.optional1,optional1),'v.optional1')
      testCase.verifyTrue(cells.isequal(v.cell_array,varargsTest.obj.cell_array),'v.cell_array') %this was not updated
      testCase.verifyTrue(cells.isequal(v.cell_str,  varargsTest.obj.cell_str  ),'v.cell_str')   %this was not updated
      testCase.verifyTrue(cells.isequal(sort(p.Parameters),...
        sort([varargsTest.par_names,{'optional1','required1','required2'}])),'p.Parameters')
      testCase.verifyTrue(cells.isequal(sort(p.UsingDefaults),...
        sort({'cell_array','cell_str'})),'p.UsingDefaults')
      testCase.verifyTrue(cells.isequal(p.Unmatched.ignored1,ignored1),'p.Unmatched.ignored1')
      testCase.verifyTrue(cells.isequal(p.Unmatched.ignored2,ignored2),'p.Unmatched.ignored2')
      testCase.verifyTrue(cells.isequal(p.Results.required1,required1),'p.required1')
      testCase.verifyTrue(cells.isequal(p.Results.required2,required2),'p.required2')
      testCase.verifyTrue(cells.isequal(p.Results.char,new_char),      'p.new_char')
      testCase.verifyTrue(cells.isequal(p.Results.numeric,new_numeric),'p.new_numeric')
      testCase.verifyTrue(cells.isequal(p.Results.optional1,optional1),'p.optional1')
      testCase.verifyTrue(cells.isequal(p.Results.cell_array,varargsTest.obj.cell_array),'p.cell_array') %this was not updated
      testCase.verifyTrue(cells.isequal(p.Results.cell_str,  varargsTest.obj.cell_str  ),'p.cell_str')   %this was not updated
      testCase.verifyTrue(cells.isequal(sink.numeric,    new_numeric),               'sink.numeric')
      testCase.verifyTrue(cells.isequal(sink.cell_array, varargsTest.obj.cell_array),'sink.cell_array') %this was not updated
    end
    function cellTest(testCase)
      testCase.verifyTrue(cells.isequal(varargsTest.obj.cell,varargsTest.cell_matrix))
    end
    function vararginTest(testCase)
      testCase.verifyTrue(cells.isequal(varargsTest.obj.varargin,varargsTest.cell_vector))
    end
    function deleteTest(testCase)
      obj1=varargsTest.obj.dup.delete(arrayfun(@(i) i.name,varargsTest.par(1:2),'UniformOutput',false));
      obj2=varargs(varargsTest.cell_matrix(3:end,:));
      testCase.verifyTrue(obj1==obj2);
    end
    function isolateTest(testCase)
      obj1=varargsTest.obj.dup.isolate(arrayfun(@(i) i.name,varargsTest.par(1:2),'UniformOutput',false));
      obj2=varargs(varargsTest.cell_matrix(1:2,:));
      testCase.verifyTrue(obj1==obj2);
    end
    function subsetTest(testCase)
      name_part='cell';
      obj1=varargsTest.obj.dup.subset(name_part);
      obj2=varargs({});
      for i=1:numel(varargsTest.S)
        if str.contains(varargsTest.S(i).name,name_part)
          obj2.set(varargs.initS(varargsTest.S(i).name,varargsTest.S(i).value,'validation',varargsTest.S(i).validation));
        end
      end
      testCase.verifyTrue(obj1==obj2);
    end
  end
end