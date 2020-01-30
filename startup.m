%% header
disp([datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' - startup.m started'])

%% setup path

%get dir of this file
local_here=fileparts(mfilename('fullpath'));
disp(['NOTICE: startup from : ',local_here])
% builddir list
local_dir_list={...
  local_here;...
  fullfile(local_here,'packages','yamlmatlab');...
};
% add legacy support
if datetime(version('-date'))<=datetime('2016-02-11')
  local_dir_list{end+1}=fullfile(local_here,'version_patching');
end
%add to path
for i=1:numel(local_dir_list)
  if exist(local_dir_list{i},'dir')
    addpath(local_dir_list{i})
  end
end
%inform user
pathhead=strsplit(path,':');
disp(['NOTICE: top 3 entries in path are:',newline,strjoin(pathhead(1:3),newline)])

%% determine project

%define global project dir
global PROJECT;
%define project yaml filename
project_filename='project.yaml';
%check if any project has been defined
if ~exist(fullfile(local_here,'project.yaml'),'file')
  disp([...
    'NOTICE: could not find any ''',project_filename,''' file in the current directory (',local_here,'). ',newline,...
    'This is most likely because this instance has just been checkout from git.'
  ])
  project_filename='default.yaml';
  disp([...
    'Loading defaults project options from ''',project_filename,'''.'
  ])
end
%loading project metadata
PROJECT=yaml.ReadYaml(project_filename);

%% options

%get metadata entries
project_entries=fieldnames(PROJECT);
%make room for user feebdack
project_msg=cell(0);c=0;
%loop over all project entries
for i=1:size(project_entries,1)
  %propagate 
  switch project_entries{i}
    case 'name'
      %do nothing
    case 'stopiferror'
      if PROJECT.stopiferror; dbstop if error; end
    case 'stopifwarning'
      if PROJECT.stopifwarning; dbstop if warning; end
    case 'recycle'
      recycle(str.logical(PROJECT.recycle,'onoff'))
    otherwise
      c=c+1;
      project_msg(c,:)={project_entries{i},':',str.show(PROJECT.(project_entries{i}))};
  end
end
%user feedback, if any
if c>0
  disp('NOTICE: the following PROJECT fields are available:')
  if size(project_msg,1)==1
    disp(strjoin(project_msg,' '))
  else
    disp(strjoin(str.columns(project_msg,'left'),newline))
  end
end

%% project dirs 

%define relevant directories
dir_list={'data','metadata'};
%loop over them
for i=1:numel(dir_list)
  %make sure this directory is there
  if ~file.exist(fullfile(local_here,dir_list{i}))
    %build possible dir location name
    linked_dirs={};
    %check for PROJECT-defined dirs
    if isfield(PROJECT,[dir_list{i},'_dir'])
      linked_dirs{end+1}=PROJECT.([dir_list{i},'_dir']); %#ok<SAGROW>
    end
    %append frequent locations
    linked_dirs{end+1}={fullfile(getenv('HOME'),dir_list{i})}; %#ok<SAGROW>
    %init success flag
    linked_flag=false;
    %loop over possible locations
    for j=1:numel(linked_dirs)
      %link obvious locations if they exist
      linked_flag=file.ln(linked_dirs{j},local_here,true);
      if linked_flag; break; end
    end
    %create dir if linking didn't work
    if ~linked_flag
      file.mkdir(fullfile(local_here,dir_list{i}));
    end
  end
  %make sure there's a project dir in the metadata
  if strcmp(dir_list{i},'metadata')
    file.mkdir(fullfile(local_here,dir_list{i},PROJECT.name),true);
  end
  % user feedback
  disp(['NOTICE: ',str.just(dir_list{i},max(cellfun(@length,dir_list))),...
    ' dir is : ',file.ls(fullfile('.',dir_list{i}),'-ldFh')])
end
%relevant user feedback
disp(['NOTICE: available metadata is : ',newline,file.ls([fullfile('.','metadata',PROJECT.name)],'-lFh')])
disp(['NOTICE: current project is : ',PROJECT.name])

%% cleanup

clear ans project_entries i c dir_list local_dir_list local_here pathhead project_filename project_msg

%% footer

disp([datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' - startup.m ended'])