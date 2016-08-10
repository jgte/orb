clear
close all

%gather all object files
scripts=dir('*.m');

for i=1:numel(scripts)
    o=strrep(scripts(i).name,'.m','');
    %no test for this script
    if strcmp(o,mfilename)
      continue
    end
    disp(['--- testing object ',o,' ---'])
    eval([o,'.test']);
end