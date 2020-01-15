%% header
disp([datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' - startup.m started'])

%% parameters

local_here=fileparts(mfilename('fullpath'));

%% collect dir list

local_dir_list={local_here,fullfile(local_here,'packages','yamlmatlab')};

if datetime(version('-date'))<=datetime('2016-02-11')
  local_dir_list{end+1}=fullfile(local_here,'version_patching');
end

%% add to path

for i=local_dir_list
  if exist(i{1},'dir') && isempty(strfind(path,i{1}))
    addpath(i{1})
  end
end

%% cleanup

clear i local_dir_list local_here local_machines

%% options

dbstop if error
% dbstop if warning
recycle('on')

%% footer

disp([datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' - startup.m ended'])
