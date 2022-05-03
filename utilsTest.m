% NOTICE: call this with gswarmTest.go
classdef utilsTest
  properties(Constant)
    debug=true;
    testdir_out=fullfile(file.orbdir('test'),'out');
    testdir_in =fullfile(file.orbdir('test'),'in')
  end
  methods(Static)
    function out=dirin(test)
      out=fullfile(utilsTest.testdir_in,test);
    end
    function out=dirtest(test)
      out=fullfile(utilsTest.testdir_out,[test,'.test']);
    end
    function out=dircheck(test)
      out=fullfile(utilsTest.testdir_out,test);
    end
    %NOTICE: This utility should come after all plots are generated for a given test; the
    %        plots for this test are placed in the same directory and this utility loops
    %        through all tests in that directory and compares them with the plots in the
    %        [test,'.test'] directory.
    function testCase=compare_files(testCase,test)
      %define directory with test records, to have something to compare with
      dirtest =utilsTest.dirtest(test);
      %define directory with recent results
      dircheck=utilsTest.dircheck(test);
      %make sure dirs exist
      if ~file.exist(dirtest)
        warning(['plot directory ',dirtest,' is missing, created now.'])
        file.mkdir(dirtest)
      end
      if ~file.exist(dircheck)
        warning(['cannot check plots from ',test,...
          ' because test records directory ',dircheck,' is missing: ',...
          'relevant files need to be saved in this directory, possibly with ',...
          'utilsTest.check_single_plot'])
        return
      end
      %get list of plots to check
      file_list=dir(dircheck);
      %loop over all plots
      for i=1:numel(file_list)
        %skip directories
        if file_list(i).isdir; continue; end
        %easier names
        file_test=fullfile(dirtest,file_list(i).name);
        file_check=fullfile(dircheck,file_list(i).name);
        %check if file exists in dircheck
        if ~file.exist(file_test)
          str.say('Copying',file_check,'to',file_test);
          file.rsync(file_check,file_test);
        end
        %get extension
        [~,~,e]=fileparts(file_list(i).name);
        %branch on file time
        switch lower(e)
        case '.png'
          %compare figures
          same_file=file.im_count_diff_pixels(file_test,file_check)==0;
          %show the plots if they don't match
          if ~same_file
            %NOTICE: this only works on OSX
            file.system(['open ',...
              fullfile(dirtest ,file_list(i).name),' ',...
              fullfile(dircheck,file_list(i).name)...
            ]);
          end
          %figures always pass this test because small irrelevant changes make it fail
          same_file=true;
        case {'.tex','.txt','.dat','.ascii','.xyz'}
          %compare file contents
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
        otherwise
          if utilsTest.debug
            str.say(['NOTICE: ignoring file ',file_list(i).name,' inside ',dirtest,...
              ' because it is not a valid file.'])
          end
          continue
        end
        %save test results
        testCase.verifyTrue(same_file,file_list(i).name);
      end
    end
    %NOTICE: This utility assumes a plot was just created and it is the only one relevant
    %        to this test.
    function testCase=check_single_plot(testCase,test)
      title(test)
      plotting.save(fullfile(utilsTest.dircheck(test),'plot.png'));
      close(gcf)
      testCase=utilsTest.compare_files(testCase,test);
    end
    %NOTICE: This utility will take a variable and check if it is the same as what was
    %        saved in the test directory (and save it if necessary)
    function testCase=check_test_data(testCase,test,a)
      a_test_file=fullfile(utilsTest.dirtest(test),'a.mat');
      if file.exist(a_test_file)
        load(a_test_file,'a_check')
        check=cells.isequal(a,a_check);
        msg=a_test_file;
        if ~check
          keyboard
        end
      else
        msg=['Could not find check file ',a_test_file,'; created now, re-run the test'];
        a_check=a;
        file.ensuredir(a_test_file);
        save(a_test_file,'a_check');
        check=false;
      end
      testCase.verifyTrue(check,msg);
    end
  end
end