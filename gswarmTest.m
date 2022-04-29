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
    function testCase=compare_plots(testCase,name)
      %define directory with test records, to have something to compare with
      dirtest =fullfile(file.orbdir('plot'),[name,'.test']);
      %define directory with recent results
      dircheck=fullfile(file.orbdir('plot'),name);
      %make sure dirs exist
      if ~file.exist(dirtest)
        str.say(['WARNING: cannot check plots from ',name,...
          ' because test records directory ',dirtest,' is missing.'])
        return
      end
      if ~file.exist(dircheck)
        str.say(['WARNING: cannot find plots from ',name,...
          ' because plot directory ',dircheck,' is missing.'])
        return
      end
      %get list of plots to check
      plot_list=dir(dircheck);
      %loop over all plots
      for i=1:numel(plot_list)
        if plot_list(i).isdir
          continue
        end
        [~,~,e]=fileparts(plot_list(i).name);
        if ~strcmp(e,'.png')
          if gswarmTest.debug
            str.say(['NOTICE: ignoring file ',plot_list(i).name,' inside ',dircheck,...
              ' because it is not a valid image file.'])
          end
          continue
        end
        same_image=file.im_count_diff_pixels(...
          fullfile(dirtest,plot_list(i).name),...
          fullfile(dircheck,plot_list(i).name)...
        )==0;
        %show the plots if they don't match
        if ~same_image
          %NOTICE: this only works on OSX
          file.system(['open ',...
            fullfile(dirtest ,plot_list(i).name),' ',...
            fullfile(dircheck,plot_list(i).name)...
          ]);
        end
        %save test results
        testCase.verifyTrue(same_image,plot_list(i).name)
      end
    end
    function testCase=compare_ascii(testCase,name,ext)
      %define directory with test records, to have something to compare with
      dirtest =fullfile(file.orbdir('plot'),[name,'.test']);
      %define directory with recent results
      dircheck=fullfile(file.orbdir('plot'),name);
      %make sure dirs exist
      if ~file.exist(dirtest)
        str.say(['WARNING: cannot check plots from ',name,...
          ' because test records directory ',dirtest,' is missing.'])
        return
      end
      if ~file.exist(dircheck)
        str.say(['WARNING: cannot find plots from ',name,...
          ' because plot directory ',dircheck,' is missing.'])
        return
      end
      %get list of plots to check
      file_list=dir(dircheck);
      %loop over all plots
      for i=1:numel(file_list)
        if file_list(i).isdir
          continue
        end
        [~,~,e]=fileparts(file_list(i).name);
        if ~strcmp(e,ext)
          if gswarmTest.debug
            str.say(['NOTICE: ignoring file ',file_list(i).name,' inside ',dircheck,...
            ' because it is not a file with extension ''',ext,'''.'])
          end
          continue
        end
        same_file=file.str_equal(...
          fullfile(dirtest, file_list(i).name),...
          fullfile(dircheck,file_list(i).name)...
        )==0;
        %show the plots if they don't match
        if ~same_file
          file.system(['diff ',...
            fullfile(dirtest, file_list(i).name),' ',...
            fullfile(dircheck,file_list(i).name)...
          ]);
        end
        %save test results
        testCase.verifyTrue(same_file,file_list(i).name)
      end
    end
  end
  methods(Test)
    function deepoceanmaskTest(testCase)
      gswarm.plotdeepoceanmask(4,fullfile(file.orbdir('plot'),'deepoceanmask','deepoceanmask.png'));
      close(gcf)
      gswarmTest.compare_plots(testCase,'deepoceanmask');
    end
    function c20modelTest(testCase)
      for i=1:numel(gswarmTest.project_list)
        %set expected values in PROJECT
        gswarmTest.update_project_name(gswarmTest.project_list{i});
        %delete it
        gswarm.c20model('clear',file.orbdir('plot'));
        %plot it
        gswarm.c20model('plot',file.orbdir('plot'));
        close(gcf);
        testCase=gswarmTest.compare_plots(testCase,'c20model');
        %create latex tables
        gswarm.c20model('latex',file.orbdir('plot'));
        testCase=gswarmTest.compare_ascii(testCase,'c20model','tex');
      end
    end
    function grace_modelTest(testCase)
      global PROJECT
      %NOTICE: the start/stop times are not taken from PROJECT.start/stop_date
      gswarm.grace_model(...
        'max_degree',gswarmTest.max_degree,...
        'force',true,...
        'start',datetime(PROJECT.start_date),...
        'stop', datetime(PROJECT.stop_date)...
      );
      close all
      gswarmTest.compare_plots(testCase,'grace_model');
    end
    %This needs to come after grace_modelTest
    function precombvalTest(testCase)
      %set expected values in PROJECT
      global PROJECT
      gswarmTest.update_project_name(gswarmTest.project_list{1});
      %call main routine, forcing re-computing everything
      [d,p]=gswarm.precombval(...
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
        testCase=gswarmTest.compare_plots(testCase,p{h}.name);
      end
    end
  end
end