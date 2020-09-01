function startup

%% header
disp([datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' - startup.m started'])

%% setup path

%get dir of this file
dirnow=fileparts(mfilename('fullpath'));
disp(['NOTICE: startup from : ',dirnow])
addpath(dirnow)
% add legacy support
if datetime(version('-date'))<=datetime('2016-02-11')
  addpath(fullfile(dirnow,'version_patching'));
end
%inform user
pathhead=strsplit(path,':');
disp(['NOTICE: top 3 entries in path are:',newline,strjoin(pathhead(1:3),newline)])

%% add packages

%make room for path dirs to be added
path_dir_list={};
%define packages dir
package_dir=fullfile(dirnow,'packages');
%get list of packages
package_list=file.find(package_dir,'-mindepth 1 -maxdepth 1 -type d');
%loop over them
for i=1:numel(package_list)
  %get package names
  [~,package_name]=fileparts(package_list{i});
  %check if this is a matlab package (starts with '+'), see:
  % https://www.mathworks.com/help/matlab/matlab_oop/scoping-classes-with-packages.html
  if package_name(1)=='+'
    %matlab will automaticall add the +<package name> dir to the path
    path_dir_list=[path_dir_list(:);{fullfile(dirnow,'packages')}];
  else
    %look for all subdirs with .m files
    mfile_list=file.find(package_list{i},'-type f -name *.m');
    mfile_dir_list=unique(cellfun(@fileparts,mfile_list,'UniformOutput',false));
    %add them to the list of dirs
    path_dir_list=[path_dir_list(:);mfile_dir_list(:)];
  end
end
%loop over all collected dirs and add them to the path
cellfun(@(i)addpath(i,'-end'),path_dir_list);
cellfun(@(i)disp(['NOTICE: added to path dir: ',i]),path_dir_list);

%% determine project

%define global project dir
global PROJECT;
%define project yaml filename
project_filename='project.yaml';
%check if any project has been defined
if ~exist(fullfile(dirnow,'project.yaml'),'file')
  disp([...
    'NOTICE: could not find any ''',project_filename,''' file in the current directory (',dirnow,'). ',newline,...
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

%define relevant directories
dir_list={'data','plot','metadata'};

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
      %do nothing (avoids the message on the project fields
    otherwise
      if ~ismember(project_entries{i},cellfun(@(i) [i,'_dir'],dir_list,'UniformOutput',false))
        c=c+1;
        project_msg(c,:)={project_entries{i},':',str.show(PROJECT.(project_entries{i}))};
      end
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

%loop over them
for i=1:numel(dir_list)
  %define target
  target_dir=fullfile(dirnow,dir_list{i});
  %make sure this directory is there
  if ~file.exist(target_dir)
    %init loop vars
    linked_flag=false;
    source_dir={};
    %check for PROJECT-defined dirs
    if isfield(PROJECT,[dir_list{i},'_dir'])
      source_dir{end+1}=PROJECT.([dir_list{i},'_dir']); %#ok<AGROW>
    end
    %try linking frequent locations
    source_dir{end+1}=fullfile(getenv('HOME'),dir_list{i}); %#ok<AGROW>
    %loop over source locations
    for j=1:numel(source_dir)
      %link them, unless already done
      if ~linked_flag && file.exist(source_dir{j})
        linked_flag=file.ln(source_dir{j},dirnow,true);
      end
    end
    %create dir if linking didn't work
    if ~linked_flag
      file.mkdir(target_dir);
    end
  end
  %make sure there's a project dir in the metadata
  if strcmp(dir_list{i},'metadata')
    file.mkdir(fullfile(target_dir,PROJECT.name),true);
  end
  % user feedback
  disp(['NOTICE: ',str.just(dir_list{i},max(cellfun(@length,dir_list))),...
    ' dir is : ',file.ls(fullfile('.',dir_list{i}),'-ldFh')])
end
%relevant user feedback
disp(['NOTICE: available metadata is : ',newline,file.ls(fullfile('.','metadata',PROJECT.name),'-lFh')])
disp(['NOTICE: current project is : ',PROJECT.name])

%% footer

disp([datestr(datetime('now'),'yyyy-mm-dd HH:MM:SS'),' - startup.m ended'])