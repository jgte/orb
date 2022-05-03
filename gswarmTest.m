% NOTICE: call this with gswarmTest.go
classdef gswarmTest < matlab.unittest.TestCase
  properties(Constant)
     start_date='2015-01-01';
     stop_date ='2019-12-31';
     project_list={'test-precombval','test-validation'};
     debug=false;
     max_degree=20;
  end
  methods(Static)
    function go(testMethod)
      global PROJECT
      PROJECT.start_date=gswarmTest.start_date;
      PROJECT.stop_date =gswarmTest.stop_date;
      t=gswarmTest;
      if ~exist('testMethod','var')
        disp(table(t.run))
      else
        disp(table(t.run(testMethod)))
      end
    end
    function update_project_name(in)
      %set expected values in PROJECT
      global PROJECT
      PROJECT.name=in;
      %NOTICE: this is needed so that the metadata_dir property of dataproduct is 
      %        refreshed with the PROJECT.name just updated in the previous line
      clear functions %#ok<CLFUNC>
    end
  end
  methods(Test)
    function testCase=deepoceanmaskTest(testCase)
      gswarm.plotdeepoceanmask(4);
      testCase=utilsTest.check_single_plot(testCase,'deepoceanmask');
    end
    function testCase=c20modelTest(testCase)
      for i=1:numel(gswarmTest.project_list)
        %set expected values in PROJECT
        gswarmTest.update_project_name(gswarmTest.project_list{i});
        %delete it
        gswarm.c20model('clear',file.orbdir('plot'));
        %plot it
        gswarm.c20model('plot',file.orbdir('plot'));
        close(gcf);
        %create latex tables
        gswarm.c20model('latex',file.orbdir('plot'));
      end
      testCase=utilsTest.compare_files(testCase,'c20model');
    end
    function testCase=grace_modelTest(testCase)
      global PROJECT
      %NOTICE: the start/stop times are not taken from PROJECT.start/stop_date
      gswarm.grace_model(...
        'max_degree',gswarmTest.max_degree,...
        'force',true,...
        'start',datetime(PROJECT.start_date),...
        'stop', datetime(PROJECT.stop_date)...
      );
      close all
      testCase=utilsTest.compare_files(testCase,'grace_model');
    end
    %This needs to come after grace_modelTest
    function testCase=precombvalTest(testCase)
      %set expected values in PROJECT
      global PROJECT
      gswarmTest.update_project_name(gswarmTest.project_list{1});
      %call main routine, forcing re-computing everything
      [d,p]=gswarm.precombval(...
        'overwrite_common_t',true,...
        'get_input_data',false,... 
        'c20model'      ,false,...
        'grace_model'   ,false,...
        'force'         ,true,...
        'plot_force'    ,true...
      );
      close all
      %test start/stop dates
      testCase.verifyTrue(cells.isequal(d.start,datetime(PROJECT.start_date)),'start_date')
      testCase.verifyTrue(cells.isequal(d.stop ,datetime(PROJECT.stop_date )),'stop_date' )
      %test common time domain
      for h=1:numel(p)
        for i=1:numel(p{h}.metadata.sources)
          %skip the parametric model
          if str.contains(p{h}.metadata.sources{i}.name,'pd.ts')
            continue
          end
          %make sure this source is available (some may not be when loading data from saved data files)
          if ~isfield(d.data,p{h}.metadata.sources{i}.codename)
            str.say(['WARNING: skipped checking time domain of ',p{h}.metadata.sources{i}.name,...
              ' because cannot find it in output ''d''.'])
            continue
          end
          for j=i+1:numel(p{h}.metadata.sources)
            %skip the parametric model
            if str.contains(p{h}.metadata.sources{j}.name,'pd.ts')
              continue
            end
            %make sure this source is available (some may not be when loading data from saved data files)
            if ~isfield(d.data,p{h}.metadata.sources{j}.codename)
              str.say(['WARNING: skipped checking time domain of ',p{h}.metadata.sources{j}.name,...
                ' because cannot find it in output ''d''.'])
              continue
            end
            disp([' -- testCase.verifyTrue -- ','t: ',p{h}.metadata.sources{i}.name,' vs ',p{h}.metadata.sources{j}.name])
            testCase.verifyTrue(all(...
              d.data.(p{h}.metadata.sources{i}.codename).signal.t == ...
              d.data.(p{h}.metadata.sources{j}.codename).signal.t...
            ),['t: ',p{h}.metadata.sources{i}.name,' vs ',p{h}.metadata.sources{j}.name])
          end
        end
      end
      %test output figures
      for h=1:numel(p)
        testCase=utilsTest.compare_files(testCase,p{h}.name);
      end
    end
  end
end