classdef kp
  properties(Constant)
    data_dir=file.resolve_home('~/data/kp')
  end
  methods(Static)
    %% https://github.com/mattkjames7/kpindex
    % these methods make use of the python package kpindex
    % the data is taken from ftp://ftp.gfz-potsdam.de/pub/home/obs/kp-ap/tab/
    % the data is saved in '~/data/kp/kpindex'
    function out=kpindex_install
      out=file.system('pip3 install DateTimeTools kpindex --user','stop_if_error',true,'disp',true);
    end
    function out=kpindex_update
      com={
        'import kpindex';...
        'kpindex.UpdateLocalData()';...
      };
      out=file.python3(strjoin(com,';'),'cd',fullfile(fp.data_dir,'kpindex'),'stop_if_error',true,'disp',true);
    end
    function out=kpindex_list
      com={
        'import kpindex';...
        'import pandas as pd';...
        'pd.set_option(''display.max_rows'', None)';...
        'pd.set_option(''display.max_columns'', None)';...
        'pd.set_option(''display.width'', None)';...
        'pd.set_option(''display.max_colwidth'',None)';...
        'print(pd.DataFrame(kpindex.GetKp(None)))';...
      };
      [~,out]=file.python3(strjoin(com,';'),'cd',fullfile(fp.data_dir,'kpindex'),'stop_if_error',true,'disp',false);
    end
    function out=kpindex_parse
      % Please set KPDATA_PATH environment variable
      % Reading Data 100.00%
      %           Date  Index   ut0   ut1        Kp        Sum  Ap   Cp Activity
      % 0      19940101      0   0.0   3.0  4.000000  31.000000  26  1.2      D5
      % 1      19940101      0   3.0   6.0  4.000000  31.000000  26  1.2      D5
      data=textscan(kp.kpindex_list,...
        '%*d %{yyyyMMdd}D %*d %f %f %f %f %f %f %s',...
        'HeaderLines',3,...
        'MultipleDelimsAsOne',true,...
        'ReturnOnError',false...
      );
      out.ut0=data{1}+hours(data{2});
      out.ut1=data{1}+hours(data{3});
      out.Kp=data{4};
      out.Sum=data{5};
      out.Ap=data{6};
      out.Cp=data{7};
      out.Activity=data{8};
    end
    function kpindex_plot
      out=kp.kpindex_parse;
      plot(out.ut0,smooth(out.Sum,81*8))
    end
    %% ftp.gfz-potsdam.de:~/pub/home/obs/Kp_ap_Ap_SN_F107
    function download
      ftpobj = ftp('ftp.gfz-potsdam.de');
      cd(ftpobj,'~/pub/home/obs/Kp_ap_Ap_SN_F107');
      mget(ftpobj,'Kp_ap_Ap_SN_F107_[0-9]*.txt',fullfile(kp.data_dir,'Kp_ap_Ap_SN_F107'));
      close(ftpobj)
    end
    function out=raw_filelist(varargin)
      out=file.unwrap(fullfile(kp.data_dir,'Kp_ap_Ap_SN_F107','Kp_ap_Ap_SN_F107_[0-9]*.txt'),'disp',false);
    end
    function out=load(varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addParameter( 'disp' , false,                                 @(i) isscalar(i) && islogical(i));
      p.addParameter( 'start', datetime(2000,1,1),                    @(i) isscalar(i) && isdatetime(i));
      p.addParameter( 'stop' , datetime(inf,'convertfrom','datenum'), @(i) isscalar(i) && isdatetime(i));
      p.parse(varargin{:})
      filelist=kp.raw_filelist(varargin{:});
      out=struct(...
        't',      [],...
        'Ap',     [],...
        'SN',     [],...
        'F107obs',[],...
        'F107adj',[],...
        'D',      []...
      );
      for i=1:numel(filelist)
        txt=filelist{i};
        mat=str.rep(txt,'.txt','.mat');
        if file.exist(mat) && file.datenum(mat)>file.datenum(txt)
          if p.Results.disp; disp(['Loading data from ',mat]); end
          load(mat,'S')
        else
          if p.Results.disp; disp(['Loading data from ',txt]); end
          fid=fopen(txt);
          S=textscan(fid,...
            '%d %d %d %*d %*f %*d %*d %*f %*f %*f %*f %*f %*f %*f %*f %*d %*d %*d %*d %*d %*d %*d %*d %f %f %f %f %f',...
            'HeaderLines',40,...
            'MultipleDelimsAsOne',true,...
            'ReturnOnError',false...
          );
          fclose(fid);
          save(mat,'S');
        end
        out.t      =[out.t;datetime(S{1},S{2},S{3})];
        out.Ap     =[out.Ap     ;S{4}];
        out.SN     =[out.SN     ;S{5}];
        out.F107obs=[out.F107obs;S{6}];
        out.F107adj=[out.F107adj;S{7}];
        out.D      =[out.D      ;S{8}];
      end
      idx=out.t>p.Results.start & out.t<p.Results.stop;
      for i={'t','Ap','SN','F107obs','F107adj','D'}
        out.(i{1})=out.(i{1})(idx);
      end
    end
    function plot(varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addParameter( 'smooth'  , days(81),   @(i) isscalar(i) && isduration(i));
      p.addParameter( 'variable', 'F107adj', @char);
      p.parse(varargin{:})
      %get the data
      d=kp.load(varargin{:});
      %determine smoothing window
      n=ceil(p.Results.smooth/mean(diff(d.t)));
      %open figure
      plotting.figure;
      %plot it  
      plot(d.t,smooth(d.(p.Results.variable),n))
      ylabel(p.Results.variable)
      plotting.enforce(...
        'plot_legend_location','none',...
        'plot_line_width',4 ...
      );
    end
  end
end