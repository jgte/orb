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
  %TODO: fix this, not sure why it triggers the following error:
  %   The class file has no Constant property or Static method named 'orbidr'.
  if any(cells.isstrequal({...
      'orbit','nrtdm','nrtdm_metadata','nrtdm_product',...
    },o))
    continue
  end
  %only test if it makes sense
  if cells.respondto({o},'test')
    disp(['--- testing object ',o,' ---'])
    eval([o,'.test']);
  end
end