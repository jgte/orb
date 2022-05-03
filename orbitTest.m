classdef orbitTest < matlab.unittest.TestCase
  properties(Constant)
    l=100;
    step=10;
    fmt={'tudelft','aiub','gswarm','swarm'};
    satn={'swarma-kin','swarmb-kin','swarmc-kin-ifg','swarma-l2'};
    data_dir=utilsTest.testdir_in;
  end
  methods(Static)
    function out=pos(l)
      out=randn(l,3);
    end
    function out=vel(l)
      out=randn(l,3)+1;
    end
    function out=acc(l)
      out=randn(l,3)+2;
    end
    function out=start
      out=datetime('2020-01-01');
    end
    function out=stop(l)
      out=orbitTest.start+seconds(orbitTest.step*l);
    end
    function out=time
      out=orbitTest.start:seconds(orbitTest.step):orbitTest.stop;
    end
    function go(testMethod)
      t=orbitTest;
      if ~exist('testMethod','var')
        disp(table(t.run))
      else
        disp(table(t.run(testMethod)))
      end
    end
  end
  methods(Test)
    function testCase=filenames(testCase)
      test='filenames';
      fmti=orbitTest.fmt;
      satni=orbitTest.satn;
      for i=1:numel(fmti)
        a.(fmti{i})=orbit.filename(...
          fmti{i},satni{i},orbitTest.start,...
          'data_dir',orbitTest.data_dir...
        );
      disp(a.(fmti{i}))
      end
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=formats(testCase)
      test='formats';
      fmti=orbitTest.fmt;
      satni=orbitTest.satn;
      for i=1:numel(fmti)
        a.(fmti{i})=orbit.import(...
          fmti{i},satni{i},orbitTest.start,orbitTest.stop(orbitTest.l),...
          'data_dir',orbitTest.data_dir...
        );
      end
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
  end
end