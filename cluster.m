classdef cluster
  properties(Constant)
    scratch=getenv('SCRATCH');
  end
  methods(Static)
    function out=istextonly
      out=usejava('jvm') && ~feature('ShowFigureWindows');
    end
    function out=stage(filelst,rootdir,sinkdir,more_args)
      if ischar(filelst)
        filelst={filelst};
      end
      if ~exist('rootdir','var') || isempty(rootdir)
        rootdir=file.fullpath('.');
      end
      if ~exist('sinkdir','var') || isempty(sinkdir)
        sinkdir=cluter.scratch;
      end
      if ~exist('more_args','var')
        more_args='';
      end
      %assert assumed types
      assert(iscellstr(filelst),['input ''filelst'' must be of class ''cellstr'', not of class ''',class(filelst),'''.'])
      assert(   ischar(rootdir),['input ''rootdir'' must be of class ''char'', not of class ''',   class(rootdir),'''.'])
      assert(   ischar(sinkdir),['input ''sinkdir'' must be of class ''char'', not of class ''',   class(sinkdir),'''.'])
      %make sure trailing slashes are there
      filelst=file.trailing_filesep(file.fullpath(filelst));
      rootdir=file.trailing_filesep(file.fullpath(rootdir ));
      sinkdir=file.trailing_filesep(file.fullpath(sinkdir ));
      %sanity
      assert(all(cellfun(@(i) ~strcmp(strrep(i,rootdir,''),i),filelst)),'Some file(s)/dir(s) are not in specified root dir.')
      %outputs
      out=cell(numel(filelst));
      %loop over all specified files
      for i=1:numel(filelst)
        if exist(filelst{i},'file')
          source=filelst{i};
          sink=strrep(filelst{i},rootdir,sinkdir);
          [out{i},s]=file.rsync(source,sink,more_args);
          assert(s==0,['failed to rsync file/dir ',source,' to ',sink])
        else
          warning(['Could not stage file/dir ',filelst{i},' because it does not exist.'])
        end
      end
      %flatten output
      out=cells.flatten(out);
    end
    function out=unstage(filelist,rootdir,sinkdir,more_args)
      if ~exist('sinkdir','var') || isempty(sinkdir)
        sinkdir=cluter.scratch;
      end
      if ~exist('more_args','var')
        more_args='';
      end
      
    end
  end
end

