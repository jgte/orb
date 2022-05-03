classdef simpledataTest < matlab.unittest.TestCase
  properties(Constant)
    l=100;
    w=4;
    debug=false;
    y_randn_scale=0.2;
    y_sin_scale=[0.8 0.5];
    y_cos_scale=[0.2 1.2];
end
  methods(Static)
    function out=args(w)
      out={...
        'units',     strcat(cellstr(repmat('y_unit-',w,1)),cellstr(num2str((1:w)'))),...
        'labels',    strcat(cellstr(repmat('label-', w,1)),cellstr(num2str((1:w)'))),...
        'x_units',    'x_unit',...
        'descriptor','original'...
      };
    end
    function out=x(l)
      out=transpose(0:l-1);
    end
    function out=T(l)
      out=l./[3 5];
    end
    function out=y_poly_scale(l)
      out=[1 1 1]./[1 l l^2];
    end
    function out=y_randn(l,w)
      out=simpledataTest.y_randn_scale*randn(l,w);
    end
    function out=col_scale(w)
      out=linspace(1,0,w+1);
      out=out(1:end-1);
    end
    function out=y_poly(l,w)
      p1=simpledataTest.y_poly_scale(l);
      out=pardecomp(...
        simpledataTest.x(l),...
        [],...
        'np',numel(p1),...
        'T',[],...
        'x',p1...
      ).y_sum*simpledataTest.col_scale(w);
        % % The code below is the same as the code above
        % c=simpledataTest.y_poly_scale;
        % x1=simpledataTest.x;
        % out=zeros(l,w);
        % for i=1:numel(c)
        %   out=out+c(i)*(x1).^(i-1)*ones(1,w);
        % end
    end
    function out=y_sin_T(l,w)
      out=pardecomp(...
        simpledataTest.x(l),...
        [],...
        'np',0,...
        'T',simpledataTest.T(l),...
        'x',[simpledataTest.y_sin_scale(l,w),zeros(size(simpledataTest.y_cos_scale(l,w)))]...
      ).y_sum*simpledataTest.col_scale(w);
      % % The code below is the same as the code above
      % c=simpledataTest.y_sin_scale;
      % T=simpledataTest.T;
      % x1=simpledataTest.x;
      % out=zeros(l,w);
      % for i=1:numel(T)
      %   out=out+c(i)*sin( ...
      %     2*pi/T(i)*x1*ones(1,w)...
      %   );
      % end
    end
    function out=y_cos_T(l,w)
      out=pardecomp(...
        simpledataTest.x(l),...
        [],...
        'np',0,...
        'T',simpledataTest.T(l),...
        'x',[zeros(size(simpledataTest.y_sin_scale(l,w))),simpledataTest.y_cos_scale(l,w)]...
      ).y_sum*simpledataTest.col_scale(w);
      % % The code below is the same as the code above
      % c=simpledataTest.y_cos_scale;
      % T=simpledataTest.T;
      % x1=simpledataTest.x',l);
      % out=zeros(l,w);
      % for i=1:numel(T)
      %   out=out+c(i)*cos( ...
      %     2*pi/T(i)*x1*ones(1,w)...
      %   );
      % end
    end
    function out=y_all(l,w)
      out=pardecomp(...
        simpledataTest.x(l),...
        [],...
        'np',numel(simpledataTest.y_poly_scale(l)),...
        'T',simpledataTest.T(l),...
        's',simpledataTest.y_sin_scale,...
        'c',simpledataTest.y_cos_scale...
      ).y_sum*simpledataTest.col_scale(w);
    end
    function out=mask(l)
      s_prev=rng(0);
      out=rand(l,1)<0.9;
      rng(s_prev);
    end
    function out=no_mask(l)
      out=true(l,1);
    end
    function out=columns(w)
%       %NOTICE: this return from 1 to w columns
%       out=[1,find(randn(1,w-1)>0)+1];
      out=1:2:w;
    end
    function out=obj(l,w)
      args=simpledataTest.args(w);
      out=simpledata(...
        simpledataTest.x(l),...
        simpledataTest.y_all(l,w),...
        'mask',simpledataTest.mask(l),...
        args{:}...
      );
    end
    function go(testMethod)
      t=simpledataTest;
      if ~exist('testMethod','var')
        disp(table(t.run))
      else
        disp(table(t.run(testMethod)))
      end
    end
  end
  methods(Test)
    function testCase=plot(testCase)
      test='simpledata.plot';
      a=simpledataTest.obj(simpledataTest.l,simpledataTest.w);
      plotting.figure('plot_visible','off');
      a.plot('columns',simpledataTest.columns(simpledataTest.w));
      testCase=utilsTest.check_single_plot(testCase,test);
    end
    function testCase=append(testCase)
      test='simpledata.append';
      l0=simpledataTest.l;
      w0=simpledataTest.w;
      a=simpledataTest.obj(l0,w0);
      args=simpledataTest.args(w0);
      columns=simpledataTest.columns(w0);
      b=a.append(...
        simpledata(-simpledataTest.l-1:-1,simpledataTest.y_all(l0+1,w0),...
          'mask',simpledataTest.mask(l0+1),...
          args{:}...
        )...
      );
      plotting.figure('plot_visible','off');
      subplot(2,1,1)
      a.plot('columns',columns)
      title('original')
      subplot(2,1,2)
      b.plot('columns',columns)
      testCase=utilsTest.check_single_plot(testCase,test);
    end
  end
end