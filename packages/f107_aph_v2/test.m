function out=test(t)
  if ~exist('t','var') || isempty(t)
    t=datetime(2003,4,1,0,0,0):days(7):datetime(2017,6,1,0,0,0);
  end

  [out.F107A, out.F107, out.APH]=f107_aph(t);
  plotting.figure;
  plot(t,out.F107A)
%   legend('F10.7A')
  ylabel('F10.7 [sfu]')
  plotting.enforce(...
    'plot_legend_location','none',...
    'plot_line_width',4 ...
  );
  
  %TODO: implement solarflux class from ftp://ftp.gfz-potsdam.de/pub/home/obs/Kp_ap_Ap_SN_F107
  
end