classdef gravityTest < matlab.unittest.TestCase
  properties(Constant)
    l=4;
    t=[datetime('now'),datetime('now')+seconds(1)];
    rep_list={'y','mat','cs','tri','mod'};
  end
  methods(Static)
    function out=nr_coeff(l)
      out=(l+1).^2;
    end
    function out=y(l)
      s_prev=rng(0);
      out=randn(1,gravityTest.nr_coeff(l));
      rng(s_prev);
    end
    function go(testMethod)
      t=gravityTest;
      if ~exist('testMethod','var')
        disp(table(t.run))
      else
        disp(table(t.run(testMethod)))
      end
    end
  end
  methods(Test)
    function testCase=reps(testCase)
      y=gravityTest.y(gravityTest.l);
      rep_val={...
        y,...
        gravity.y2mat(y),...
        gravity.mat2cs(gravity.y2mat(y)),...
        gravity.cs2tri(gravity.mat2cs(gravity.y2mat(y))),...
        gravity.cs2mod(gravity.mat2cs(gravity.y2mat(y)))...
      };
      reps=gravityTest.rep_list;
      for i=1:numel(reps)
        for j=1:numel(reps)
          out=gravity.dtc(reps{i},reps{j},rep_val{i});
          switch reps{j}
          case 'cs'
            c=any(any([out.C,out.S] ~= [rep_val{j}.C,rep_val{j}.S]));
          otherwise
            c=any(any(out~=rep_val{j}));
          end
          testCase.verifyFalse(c,['failed data type conversion between ''',reps{i},''' and ''',reps{j},'''.'])
        end
      end
    end
    function testCase=unit(testCase)
      test='gravity.unit';
      ti=gravityTest.t;
      a.unit=gravity.unit_amplitude(gravityTest.l,'t',ti);
      disp('- C')
      a.cs_C=a.unit.cs(numel(ti)).C;
      disp(a.cs_C)
      disp('- S')
      a.cs_S=a.unit.cs(numel(ti)).S;
      disp(a.cs_S)
      disp('- tri')
      a.tri=a.unit.tri{numel(ti)};
      disp(a.tri)
      disp('- mod')
      a.mod=a.unit.mod{numel(ti)};
      disp(a.mod)
      disp('- das')
      a.das_dat=a.unit.das;
      disp(a.das_dat)
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=unit_rms(testCase)
      test='gravity.unit_rms';
      ti=gravityTest.t;
      a.unit=gravity.unit_amplitude(gravityTest.l,'t',ti);
      disp('- tri')
      a.tri=a.unit.tri{numel(ti)};
      disp(a.tri)
      disp('- drms')
      a.drms=a.unit.at(ti(numel(ti))).drms;
      disp(a.drms)
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=r(testCase)
      test='gravity.r';
      ti=gravityTest.t;
      a.unit_amplitude=gravity.unit_amplitude(gravityTest.l,'t',ti);
      disp('- tri: start')
      a.start=a.unit_amplitude.tri{numel(ti)};
      disp(a.start)
      disp('- tri: 2*R')
      a.end=a.unit_amplitude.scale(a.unit_amplitude.R*2,'R').tri{numel(ti)};
      disp(a.end)
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=gm(testCase)
      test='gravity.gm';
      ti=gravityTest.t;
      a.unit_amplitude=gravity.unit_amplitude(gravityTest.l,'t',ti);
      disp('- tri: start')
      a.start=a.unit_amplitude.tri{numel(ti)};
      disp(a.start)
      disp('- tri: 2*GM')
      a.end=a.unit_amplitude.scale(a.unit_amplitude.GM*2,'GM').tri{numel(ti)};
      disp(a.end)
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=minus(testCase)
      test='gravity.minus';
      ti=gravityTest.t;
      a.unit_amplitude=gravity.unit_amplitude(gravityTest.l,'t',ti);
      disp('- tri: a=')
      a.start=a.unit_amplitude.tri{numel(ti)};
      disp(a.start)
      disp('- tri: a-a.scale(2)=')
      a.end=a.unit_amplitude-a.unit_amplitude.scale(2);
      a.end=a.end.tri{numel(ti)};
      disp(a.end)
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=grid(testCase)
      test='gravity.grid';
      plotting.figure('plot_visible','off');
      s_prev=rng(0);
      gravity.unit_randn(120,'t',gravityTest.t).grid.imagesc;
      rng(s_prev);
      testCase=utilsTest.check_single_plot(testCase,test);
    end
    function testCase=mascons(testCase)
      a.grav=gravity.CSR_Mascons;
      disp('- print gravity')
      a.grav.print
      disp('- print grid')
      a.grid=a.grav.grid;
      a.grid.print
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,'gravity.mascons',a);
    end
    function testCase=stats(testCase)
      test='gravity.stats';
      a.start=gravity.static('ggm05g').setC(0,0,0).setC(2,0,0);
      stats_list={'dmean','cumdmean','drms','cumdrms','dstd','cumdstd','das','cumdas'};
      for i=1:numel(stats_list)
        a.(stats_list{i})=a.start.(stats_list{i});
      end
      testCase=utilsTest.check_test_data(testCase,test,a);
    end
    function testCase=c(testCase)
      ti=gravityTest.t;
      li=gravityTest.l;
      if numel(ti)<3
        now=juliandate(datetime('2000-01-01'),'modifiedjuliandate');
        ti=datetime(now,        'convertfrom','modifiedjuliandate'):...
           datetime(now+li*10-1,'convertfrom','modifiedjuliandate');
      end
      s_prev=rng(0);
      a.start=gravity.unit_randn(li,'t',ti);
      d=round(rand*li);
      o=round(rand*2*d)-d;
      rng(s_prev);
      disp('- tri: a=')
      disp(a.start.tri{numel(ti)})
      disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(ti(numel(ti))),')=',num2str(a.start.C(d,o,ti(numel(ti))))])
      disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(ti(1)        ),')=',num2str(a.start.C(d,o,ti(1))        )])
      v=9.9999;
      a.end=a.start.setC(d,o,v,ti(numel(ti)));
      disp(['- tri: a.setC(',num2str(d),',',num2str(o),')=',num2str(v)])
      disp(a.end.tri{numel(ti)})
      disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(ti(numel(ti))),')=',num2str(a.end.C(d,o,ti(numel(ti))))])
      disp(['- a.C(',num2str(d),',',num2str(o),',',datestr(ti(1)        ),')=',num2str(a.end.C(d,o,ti(1))        )])
      %check if test passes
      testCase=utilsTest.check_test_data(testCase,'gravity.c',a);
    end
    function testCase=smoothing(testCase)
      test='gravity.smoothing';
      ti=gravityTest.t;
      li=gravityTest.l;
      a=gravity.unit(li);
      if isdatetime(ti)
        ti=round(li/2);
      end
      methods={'gauss','spline','trunc'};
      plotting.figure('plot_visible','off');
      for i=1:numel(methods)
        a.scale_plot(ti,methods{i});
      end
      legend(methods);
      testCase=utilsTest.check_single_plot(testCase,test);
    end
    function testCase=deepoceanmaskplot(testCase)
      test='gravity.deepoceanmaskplot';
      li=gravityTest.l;
      a=simplegrid.unit(li*90,li*45);
      a=a.spatial_mask('deep ocean');
      c=gray(li*16);
      c=flipud(c(li*4:li*12,:));
      plotting.figure('plot_visible','off');
      a.imagesc;
      colormap(c)
      colorbar off
      plotting.enforce('plot_legend_location','none');
      testCase=utilsTest.check_single_plot(testCase,test);
    end
    function testCase=functionals(testCase)
      test='gravity.functionals';
      funcs=gravity.functionals;
      a=gravity.static('ggm05g');
      plotting.figure('plot_visible','off');
      for i=1:numel(funcs)
        a.plot('method','drms','functional',funcs{i});
      end
      plotting.enforce(...
        'plot_legend_location','eastoutside',...
        'plot_legend',cellfun(@(i) gravity.functional_names(i),funcs,'UniformOutput',false),...
        'plot_ylabel','see legend'...
      );
      testCase=utilsTest.check_single_plot(testCase,test);
    end
    function testCase=funct_circ(testCase)
      funcs=gravity.functionals;
      a=gravity.unit_randn(gravityTest.l);
      for i=1:numel(funcs)
        m0=a.scale(funcs{i},'functional');
        for j=i+1:numel(funcs)
          m1=m0.scale(funcs{j},'functional').scale(funcs{i},'functional');
          dm=m1-m0;
          dm_das=dm.das/m0.das;
          disp(str.tablify(14,'Checking ',funcs{i},' -> ',funcs{j},' -> ',funcs{i},'rel. error = ',num2str(max(abs(dm_das)))))
          tol=1e-15;
          testCase.verifyTrue(all(abs(dm_das)<tol),['Circular check failed for ',...
            funcs{i},' -> ',funcs{j},' -> ',funcs{i},...
            ': relative error (',num2str(max(abs(dm_das))),...
            ' larger than tolerance (',num2str(tol),')']);
        end
      end
    end
    function testCase=sh2grid(testCase)
      %define the test names
      s=struct(...
        'test',{'gravity.sh2grid_multiple','gravity.sh2grid_single'},...
        'file',{'GSWARM_GF_SABC_COMBINED_2020-1*_09','GSWARM_GF_SABC_COMBINED_2020-01_09'}...
      );
      functional='eqwh';
      for i=1:numel(s)
        test=s(i).test;
        modelfile=fullfile(utilsTest.dirin(test),[s(i).file,'.gfc']);
        %call the routine to be tested
        m=gravity.sh2grid(...
          modelfile,...
          'static',            'ggm05c',...
          'smoothing_radius',     750e3,...
          'functional',      functional,...
          'spatial_step',             1,...
          'lmax',                    20,...
          'C20_replacement',  'GSFC5x5' ...
        );
        %need to move the resulting files to the test directory because they are created in 
        %in the same dir as the input files
        file.ensuredir(utilsTest.dircheck(test),false);
        file.system(['mv -fv ',...
          fullfile(utilsTest.dirin(test),[s(i).file,'.*.xyz']),' ',...
          utilsTest.dircheck(test)],...
          'stop_if_error',true,...
          'disp',true...
        );
        for j=1:m.length
          %plot it
          plotname=fullfile(utilsTest.dircheck(test),[strrep(s(i).file,'2020-1*',datestr(m.t(j),'yyyy-mm-dd')),'.png']);
          plotting.figure('plot_visible','off');
          m.at(m.t(j)).imagesc('cb_title',gravity.functional_label(functional));
          plotting.enforce(...
            'plot_legend_location','none',...
            'plot_ylabel','none','plot_xlabel','none',...
            'plot_colormap','jetzero',...
            'plot_caxis',[-0.5,0.5],...
            'plot_title',['G-Swarm RL01 ',datestr(m.t(j),'yyyy-mm-dd')]...
          );
          plotting.no_ticks
          plotting.save(plotname)
        end
        close all
        %test it
        testCase=utilsTest.compare_files(testCase,test);
      end
    end
  end
end