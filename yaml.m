classdef yaml
  properties(Constant)
    package='yamlmatlab';
  end
  methods(Static)
    function out=scriptdir
      persistent scriptdir_value
      if isempty(scriptdir_value)
        scriptdir_value=fileparts(which(mfilename));
        if isempty(scriptdir_value)
          scriptdir_value='.';
        end
      end
      out=scriptdir_value;
    end
    function out=packagedir
      persistent packagedir_value
      if isempty(packagedir_value)
        packagedir_value=fullfile(yaml.scriptdir,yaml.package);
      end
      out=packagedir_value;
    end
    function out=addpath
      persistent isaddpath
      if isempty(isaddpath) || ~isaddpath
        if ~cells.isincluded(strsplit(path,':'),yaml.package)
          out=addpath(genpath(yaml.packagedir));
        end
        isaddpath=true;
      end
    end
    function out=read(filename)
      yaml.addpath;
      if exist([filename,'.mat'],'file') ~= 0
        load([filename,'.mat'],'S');
        out=S; %#ok<NODEF>
      elseif exist(filename,'file') ~= 0
        out=ReadYaml(filename);
        S=out; %#ok<NASGU>
        save([filename,'.mat'],'S')
      else
        out=struct([]);
      end
    end
    function write(filename,in)
      yaml.addpath;
      WriteYaml(filename,in);
      S=in; %#ok<NASGU>
      save([filename,'.mat'],'S')
    end
    function data=update(filename,data,mode)
      if ~exist('mode','var') || isempty(mode)
        mode='both';
      end
      %load file
      file=yaml.read(filename);
      switch lower(mode)
      case 'file'
        %copy data to file
        file=structs.copy(data,file);
        %write to file
        write_to_file=true;
      case 'data'
        %copy file to data
        data=structs.copy(file,data);
        %no need to write to file
        write_to_file=false;
      case 'both'
        %copy file to data to file
        file=structs.copy(structs.copy(file,data),file);
        %write to file
        write_to_file=true;
      otherwise
        error(['Unknown mode ''',mode,'''.'])
      end
      if write_to_file
        yaml.write(filename,file);
      end
    end
  end
end