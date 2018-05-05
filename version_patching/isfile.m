function out=isfile(listFiles)
% ISFILE(LISTFILES) tests if files exist returns a list of boleans for each file.

if ~iscell(listFiles)
   listFiles={listFiles};
end
out=false(size(listFiles));
for i=1:numel(listFiles)
    out(i)= exist(listFiles{i},'file') && ~exist(listFiles{i}, 'dir');
end
