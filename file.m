classdef file
  properties(Constant)
    %NOTICE: this used to be called 'DATE_PLACE_HOLDER'
    dateplaceholder='DATE_PLACEHOLDER';
    archivedfilesext={'.gz','.gzip','.z','.zip','.tgz','.tar.gz','.tar'};
    %NOTICE: these need to be hard-coded because the point is to translate paths with any
    %        of these sub-strings into $HOME. Append as needed.
    %TODO: move this to project.yaml, if it exists
    homes={...
      '~';...
      '/home1/00767/byaa676';...
      '/Users/teixeira';...
    };
    build_particle_char='_';
    build_element_char='.';
  end
  methods(Static)
    %% dir management
    %returns the directory this script, or 'scriptname', sits in
    %uses the which command so it must be in matlab's path to work as expected
    function out=scriptdir(scriptname)
      if ~exist('scriptname','var') || isempty(scriptname)
        scriptname = mfilename;
      end
      out=fileparts(which(scriptname));
      if isempty(out)
        out='.';
      end
    end
    %appends file.scriptdir (with no arguments, i.e. the dir of file.m) to scriptname
    function out=orbscriptpath(scriptname)
      out=fullfile(file.scriptdir,scriptname);
    end
    %returns the path of orb directories, which are by default sitting on the same dir
    %as this script; typical dirs are: auxiliary, packages, metadata, data, plot (the latter 2 are not in git)
    function out=orbdir(type,ensure_is_in_project)
      if ~exist('ensure_is_in_project','var') || isempty(ensure_is_in_project)
        ensure_is_in_project=false;
      end
      global PROJECT
      if isfield(PROJECT,[type,'_dir'])
        out=PROJECT.([type,'_dir']);
        if isfield(PROJECT,'dir')
          %NOTICE: file.absolutepath handles the case out is already absolute
          out=file.absolutepath(out,PROJECT.dir);
        else
          %NOTICE: file.absolutepath defaults to file.scriptdir if no 2nd arg given
          out=file.absolutepath(out);
        end
      else
        assert(~ensure_is_in_project,['Need the directory ''',type,'_dir',...
          ''' to be defined in the project configuration file ''',PROJECT.source,'''.'])
        out=file.orbscriptpath(type);
      end
      %special cases
      switch type
      case 'metadata'
        if isfield(PROJECT,'name')
          out=fullfile(out,PROJECT.name);
        end
      end
    end
    %% file IO utils
    function valid=isfid(fid)
      %https://www.mathworks.com/matlabcentral/newsreader/view_thread/85047
      valid=fid>2;
      if valid
        try
          valid=isfinite(ftell(fid));
        catch
          valid=false;
        end
      end
    end
    function empty=isempty(filename)
      %open the file if needed
      [fid,~,close_file]=file.open(filename);
      %try to read one byte
      empty=(fseek(fid,1,'bof') ~= 0);
      %close the file (if filename give, otherwise rewind)
      if close_file
        fclose(fid);
      else
        frewind(fid);
      end
    end
    function [out,faction]=im_count_diff_pixels(img1,img2)
      A1=imread(img1);
      A2=imread(img2);
      if numel(A1)==numel(A2)
        out=sum(A1(:)~=A2(:));
        faction=out/numel(A1);
      else
        out=max([numel(A1),numel(A2)]);
        faction=1;
      end
    end
    function [fid,filename,close_file,msg]=open(filename,perm,varargin)
      if ~exist('perm','var') || isempty(perm)
        perm = 'r';
      end
      %open the file (if a filename is given)
      if file.isfid(filename)
        %rename input
        fid=filename;
        %get name of the already open file
        filename=fopen(fid);
        %don't close this file afterwards
        close_file=false;
        %patch outputs
        msg=['fid ',num2str(fid),' attributed to file ''',filename,'''.'];
      else
        %resolve mat/compressed file names
        filename=file.unwrap(filename,...
          'op_sequence', {...
            @file.resolve_wildcards;...
            @file.resolve_compressed;...
            @file.resolve_ext;...
          },...
          varargin{:},...
          'scalar_as_strings',true...
        );
        %make sure dir exists
        d=fileparts(filename);
        if ~contains(perm,'r')
          file.mkdir(d);
        end
        %open the file
        [fid,msg]=fopen(filename,perm);
        assert(fid>0,['error opening ',filename,': ',msg])
        %close this file afterwards
        close_file=true;
      end
    end
    function [header,nlines]=header(filename,max_header_len)
      %maximum number of lines to scan for header
      if ~exist('max_header_len','var') || isempty(max_header_len)
        max_header_len = 50;
      end
      %open the file if needed
      [fid,filename,close_file]=file.open(filename);
      %trivial call
      if file.isempty(fid)
        nlines = 0;
        header='';
        return
      end
      %init header
      numeric_flag = false(max_header_len,1);
      nr_columns=zeros(size(numeric_flag));
      header=cell(size(numeric_flag));
      %looping over the first lines of this file, to gather info
      for i=1:max_header_len
        header{i} = fgetl(fid);
        %error checking
        if isempty(header{i}) && ~isempty(ferror(fid))
          error(['error reading header of file ',filename,'. Error message:',newline,ferror(fid)])
        end
        %end-of-file checking
        if header{i}==-1
          header{i}='';
          numeric_flag(i) = true;
          nr_columns(i)=0;
          break
        else
          numeric_flag(i) = isnumstr(header{i});
          nr_columns(i)=numel(strsplit(header{i}));
        end
      end
      %checking we'e reached the part with numeric data
      if ~isnumstr(fgetl(fid))
        nlines=-1;
        %maybe this file has some columns with non-numeric data, so checking for changes in the number of columns
        for i=max_header_len-1:-1:1
          if nr_columns(i)~=nr_columns(end)
            nlines = i;
            break
          end
        end
        %sanity
        assert(nlines>0,['Cannot find any numeric data in ',filename])
        %if that only happens on the last line, then there's nothing I can do
        if nlines==(max_header_len-1)
          error([...
            ': file ',filename,' has more header lines than max search value (',num2str(max_header_len),').'...
          ]);
        end
      else
        %getting number of header lines (numeric data interleaved with header lines is treated as part of the header)
        nlines = find(numeric_flag,1,'first')-1;
      end
      %check if there's a header
      if isempty(nlines)
        nlines = 0;
      end
      %collapse the string cell array with the header to a plain string
      header=strjoin(header(1:nlines),'\n');
      %close the file (if fid not given)
      if close_file, fclose(fid); end
    end
    function [data,header] = textscan(filename,format,nlines,CommentStyle)
      % [DATA,HEADER] = TEXTSCAN(FILENAME,NLINES) reads numerical data from an ASCII text
      % file with name FILENAME, ignoring header lines, which are considered so
      % if holding any form of non-numerical data.
      %
      %   Input FILENAME can be a valid file ID (as returned from fopen).
      %
      %   The optional input <nlines> determines how many lines are read. If ommited,
      %   empty or -1, the whole file is returned.
      %
      %   This function makes use of <textscan>, which determines automatically the
      %   number of columns and lines (if empty <nlines>) of the output <data>.
      %
      %   Example:
      %       a=textscanh(<some file>)

      % Created by J.Encarnacao <J.G.deTeixeiradaEncarnacao@tudelft.nl>

      %input sanity
      if ~exist('CommentStyle','var') || isempty(CommentStyle)
        CommentStyle = '';
      end
      if ~exist('nlines','var') || isempty(nlines)
        nlines = -1;
      end
      if ~exist('format','var')
        format = '';
      end
      if ~exist('filename','var') || isempty(filename)
        error('need one input file.')
      end
      %open the file if needed
      [fid,~,close_file]=file.open(filename);
      %get header
      [header,nr_header_lines]=file.header(fid);
      %trivial call
      if file.isempty(fid)
        data = [];
        return
      end
      %read the data
      data  = textscan(fid,format,nlines,...
        'headerlines',nr_header_lines,...
        'returnonerror',0,...
        'emptyvalue',0,...
        'CommentStyle',CommentStyle...
      );
      %need to be sure that all columns have equal length
      len = cellfun(@(i) size(i,1),data);
      %cropping so that is not the case
      if any(len~=len(1))
        for i=1:numel(data)
          str.say('say_stack_delta',1,['NOTICE: removed ',num2str(size(data{i},1)-min(len)),' entries from the ',...
            num2str(i),'-th column of the data in file ''',filename,'''.'])
          data{i} = data{i}(1:min(len),:);
        end
      end
      %numerical output, so transforming into numerical array
      data=[data{cellfun(@(i) isnumeric(i),data)}];
      %close the file (if fid not given)
      if close_file, fclose(fid); end
    end
    %quick-saving str into filename
    function out=strsave(filename,str)
      %input sanity
      if ~exist('filename','var') || isempty(filename)
        error('need filename.')
      end
      if ~exist('str','var') %str can be empty
        error('need str.')
      end
      %open the file if needed
      [fid,~,close_file]=file.open(filename,'w');
      %save the string
      out=fprintf( fid, '%s', str);
      %close the file (if fid not given)
      if close_file, fclose(fid); end
    end
    %quick-loading str from filename
    function str=strload(filename)
      %input sanity
      if ~exist('filename','var') || isempty(filename)
        error('need filename.')
      end
      %open the file if needed
      [fid,~,close_file]=file.open(filename,'r');
      %load the string
      str=fscanf( fid, '%s');
      %close the file (if fid not given)
      if close_file, fclose(fid); end
    end
    %check if text files are the same
    function out=str_equal(f1,f2)
      s1=file.strload(f1);
      s2=file.strload(f2);
      out=numel(s1)==numel(s2) && ~any((s1.^2-s2.^2)~=0);
    end
    %% resolves filenames that exist in multiple machines, each one with a home directory listed in file.homes
    function io=resolve_home(io)
      for i=1:numel(file.homes)
        io=str.rep(io,file.homes{i},getenv('HOME'));
      end
    end
    function io=unresolve_home(io)
      io=str.rep(io,'~',getenv('HOME'));
    end
    %% mat file saving/loading
    function [out,loaded_flag]=load_mat(filename,varargin)
      p=machinery.inputParser;
      p.addRequired( 'filename',          @(i) ischar(i) || iscellstr(i));
      p.addParameter('data_var'   ,'out', @ischar);
      p.parse(filename,varargin{:})
      %check if mat file is available
      datafile=file.ext(filename,'.mat','set');
      if file.exist(datafile)
        load(datafile,p.Results.data_var)
        %sanity on the loaded data
        if ~exist(p.Results.data_var,'var')
          error(['expecting to load variable ''',p.Results.data_var,''' from file ',datafile,'.'])
        end
        if ~strcmp(p.Results.data_var,'out')
          %propagate p.Results.data_var to 'out'
          out=eval(p.Results.data_var);
        end
        loaded_flag=true;
      else
        out=[];
        loaded_flag=false;
      end
    end
    % wrapper for saving mat file (NOTICE: the non-mat file is given in 'filename')
    function save_mat(out,filename,varargin)
      p=machinery.inputParser;
      p.addRequired( 'filename',       @(i) ischar(i));
      p.addParameter('save_mat', true, @(i) isscalar(i) && islogical(i))
      p.addParameter('data_var','out', @ischar);
      p.parse(filename,varargin{:})
      %save mat file if requested
      if p.Results.save_mat && ~isempty(out)
        if ~strcmp(p.Results.data_var,'out')
          %propagate 'out' to p.Results.data_var
          eval([p.Results.data_var,'=out']);
        end
        %save the data with the requested variable name (properly propagated above)
        try
          save(file.ext(filename,'.mat','set'),p.Results.data_var)
        catch ME
          disp(['Could not save ''',mat_filename,'''; error: ',ME.message])
        end
      end
    end
    %% mat-file resolving (also works with any other extension)
    %checks if is mat file
    function out=isext(in,ext) %used to be called ismat
      if iscellstr(in)
        out=cellfun(@(i) file.isext(i,ext),in);
      else
        [~,~,e]=fileparts(in);
        out=strcmp(e,ext);
      end
    end
    %add (remove) '.mat' to (from) a file name unless it already has it (not)
    %NOTICE: always returns a cellstr
    function io=ext(io,ext,direction) %used to be declared as io=mat(io,direction,ext)
      if ~exist('direction','var') || isempty(direction)
        direction='get';
      end
      if iscellstr(io)
        io=cellfun(@(i) file.ext(i,ext,direction),io,'UniformOutput', false);
      else
        if ext(1)~='.'
          warning(['Expecting input ''ext'' to start with a period, not ''',ext,'''.'])
        end
        switch lower(direction)
        case 'set'; if ~file.isext(io,ext); io=[io,ext]; end
        case 'get'; if  file.isext(io,ext); io=io(1:end-length(ext));          end
        otherwise; error(['Cannot handle input ''direction'' with value ''',direction,'''.'])
        end
      end
    end
    function io=replace_ext(io,new_ext,varargin)
      p=machinery.inputParser;
      p.addParameter('default_dir',  '.', @ischar);
      p.parse(varargin{:})
      if iscellstr(io)
        io=cellfun(@(i) file.replace_ext(i,new_ext,p.Results.default_dir),io,'UniformOutput', false);
      else
        [d,f]=fileparts(io);
        %plug default dir, if no dir is given
        if isempty(d) || strcmp(d,'.')
          d=p.Results.default_dir;
        end
        %patch leading dot
        if new_ext(1)~='.', new_ext=['.',new_ext]; end
        %build filename with extension replaced
        io=fullfile(d,[f,new_ext]);
      end
    end
    % prefers mat or non-mat file in a file list (that may mix mat with non-mat files)
    % a mat file has extension .mat, a non-mat file has a different extension (strangely enough)
    % this means that if in={'file1.dat','file1.dat.gz'} and file1.dat.mat exists, this function will
    % return {'file1.dat.mat','file1.dat.gz'} because it looks for file1.dat.gz.mat as the corresponding mat file.
    function [out,no_change_flag]=resolve_ext(in,ext,varargin)
      p=machinery.inputParser;
      p.addRequired( 'in',                            @(i) ischar(i)   || iscellstr(i));
      p.addRequired( 'ext',                           @ischar);
      p.addParameter('prefer_non_ext_files',   false, @(i) islogical(i) && isscalar(i));
      p.addParameter('prefer_ext_files',       false, @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',      false, @(i) islogical(i) && isscalar(i));
      p.addParameter('debug',                  false, @islogical);
      p.addParameter('warn',                    true, @islogical);
      p.parse(in,ext,varargin{:})
      %sanity
      assert(~(p.Results.prefer_non_ext_files && p.Results.prefer_ext_files),'Cannot both prefer non-ext and ext files')
      if iscellstr(in)
        out=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_ext(i,ext,varargin{:}),in,'UniformOutput', false)...
        )));
        %convert to string if requested
        out=file.ensure_scalar(out,p.Results.scalar_as_strings);
      else
        %NOTICE: it is assumed that in/out are scalar below
        %first get non-ext file
        nonext=file.ext(in,ext,'get'); nonextexists=file.exist(nonext);
        if p.Results.debug; str.say('nonext:',nonext); end
        %get ext file name
        yesext=file.ext(nonext,ext,'set'); yesextexists=file.exist(yesext);
        if p.Results.debug; str.say('yesext:',yesext); end
        %if there's only the ext file, then preference is irrelevant
        if ~nonextexists && yesextexists
          out=yesext;
          if p.Results.debug; str.say('picked yesext:',yesext,'(nonext file non-existing:',nonext,')'); end
        %if there's only the non-ext file, then preference is irrelevant
        elseif nonextexists && ~yesextexists
          out=nonext;
          if p.Results.debug; str.say('picked nonext:',nonext,'(yesext file non-existing:',yesext,')'); end
        %if there's a ext file and that is prefered, take it
        elseif yesextexists && p.Results.prefer_ext_files
          out=yesext;
          if p.Results.debug; str.say('picked yesext:',yesext,'(preferred)'); end
        %if there's a non-ext file and that is prefered, take it
        elseif nonextexists && p.Results.prefer_non_ext_files
          out=nonext;
          if p.Results.debug; str.say('picked nonext:',nonext,'(preferred)'); end
        %if there's no preferences and both files exists, return the newest
        elseif yesextexists && nonextexists
          if p.Results.debug
            str.say('mat:',datetime(file.datenum(yesext),'convertfrom','datenum'),yesext   )
            str.say('out:',datetime(file.datenum(nonext),'convertfrom','datenum'),nonext)
          end
          if file.datenum(yesext)>file.datenum(nonext)
            out=yesext;
            if p.Results.debug; str.say('picked yesext:',yesext,'(newer)'); end
          else
            out=nonext;
            if p.Results.debug; str.say('picked nonext:',nonext,'(newer)'); end
          end
        %if neither file exists, make some noise
        else
          out=in;
          if p.Results.warn; str.say(['WARNING: cannot find neither ''',yesext,''' nor ''',nonext,'''.']); end
        end
      end
      %handle additional outputs
      if nargout>1
        no_change_flag=str.iseq(out,in);
      end
    end
    %% resolves wildcarded filenames
    function out=iswildcarded(in)
      out=contains(in,'*');
    end
    function [io,wildcarded_flag] = translate_wildcard(io)
      wildcarded_flag=file.iswildcarded(io);
      if wildcarded_flag
        io=str.rep(io,'.*','<dotasterisc>','*','.*','<dotasterisc>','.*');
      end
      %handle special characters, reduce them to literals
      special_chars={'+'};
      for i=1:numel(special_chars)
        if contains(io,special_chars{i})
          io=str.rep(io,special_chars{i},['\',special_chars{i}]);
        end
      end
    end
    function [out,wildcarded_flag]=resolve_wildcards(in,varargin)
      p=machinery.inputParser;
      p.addRequired( 'in',                         @(i) ischar(i)   || iscellstr(i));
      p.addParameter('disp',                true,  @(i) islogical(i) && isscalar(i));
      p.addParameter('directories_only',    false, @(i) islogical(i) && isscalar(i));
      p.addParameter('files_only',          false, @(i) islogical(i) && isscalar(i));
      p.addParameter('stop_if_empty',       false, @(i) islogical(i) && isscalar(i));
      p.addParameter('passtrough_if_empty', false, @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',   false, @(i) islogical(i) && isscalar(i));
      p.addParameter('say_stack_delta',find(arrayfun(@(i) strcmp(i.file,'file.m'),dbstack),1,'last'), @num.isscalar);
      % parse it
      p.parse(in,varargin{:});
      %reduce scalar cell to char (avoids infinite loops with vector mode)
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        out=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_wildcards(i,varargin{:}),in,'UniformOutput',false)...
        )));
      else
        %NOTICE: it is NOT assumed that out are scalar below
        %trivial call
        if isempty(in)
          if p.Results.scalar_as_strings
            out='';
          else
            out={''};
          end
          return
        end
        %wildcard character '*' needs to be translated to '.*' (if not already)
        [in,wildcarded_flag]=file.translate_wildcard(file.resolve_home(in));
        %finding relevant directory
        [a,f,e]=fileparts(in);
        if isempty(a)
            a='.';
        end
        %rebuilding complete filename
        f_in=[f,e];
        %check if a simple directory has been passed
        if isempty(f_in) && ~wildcarded_flag
          if p.Results.scalar_as_strings
            out=a;
          else
            out={a};
          end
          return
        end
        %fetching directory listing
        f=dir(a);
        f=struct2cell(f);
        %branch on options
        if p.Results.directories_only
          idx=[f{end-1,:}];
        elseif p.Results.files_only
          idx=~[f{end-1,:}];
        else
          idx=true(1,size(f,2));
        end
        f=f(1,idx);
        %greping
        f1=reduce_cell_array(regexp(f(:),['^',f_in,'$'],'match'));
        %inform user if no files found
        if isempty(f1)
          msg=['found no files named ''',in,'''.'];
          if p.Results.stop_if_empty
            error(msg)
          end
          if p.Results.passtrough_if_empty
            out=fullfile(a,f_in);
            return
          end
          if p.Results.disp
            str.say('say_stack_delta',p.Results.say_stack_delta,['WARNING: ',msg])
          end
          out='';
          return
        end
        %add path
        root=fileparts(in);
        out=cellfun(@(i) fullfile(root,i),f1,'UniformOutput',false);
        %inform user
        if p.Results.disp && wildcarded_flag
          str.say('say_stack_delta',p.Results.say_stack_delta,['found ',num2str(numel(out)),' filenames in wildcarded string ''',in,'''.'])
        end
      end
      %convert to string if requested
      out=file.ensure_scalar(out,p.Results.scalar_as_strings);
    end
    %% resolves compressed files (expanding them as needed)
    function out=iscompressed(in)
      if iscellstr(in)
        out=cellfun(@(i) file.iscompressed(i),in,'UniformOutput', false);
      else
        out=false;
        for i=file.archivedfilesext
          out=file.isext(in,i{1});
          if out; return; end
        end
      end
    end
    function out=resolve_compressed(in,varargin)
      p=machinery.inputParser;
      p.addRequired( 'in',                            @(i) ischar(i)   || iscellstr(i));
      p.addParameter('prefer_compressed_files',false, @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',      false, @(i) islogical(i) && isscalar(i));
      p.addParameter('say_stack_delta',find(arrayfun(@(i) strcmp(i.file,'file.m'),dbstack),1,'last'), @num.isscalar);
      p.addParameter('debug'            , false, @islogical);
      p.parse(in,varargin{:})
      %reduce scalar
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        out=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_compressed(i,varargin{:},'scalar_as_strings',true),in,'UniformOutput',false)...
        )));
        %convert to string if requested
        out=file.ensure_scalar(out,p.Results.scalar_as_strings);
      else
        %NOTICE: it is assumed that in/out are scalar below
        %resolve homes
        in=file.resolve_home(in);
        %branch on type of file
        if file.iscompressed(in)
          %need fileparts
          [~,f,e]=fileparts(in);
          %now handling compressed files: prepend tar to extension if there
          if strcmp(f(max([1,end-3]):end),'.tar')
            e=['.tar',e];
          end
          %check if uncompressed is already available
          in=file.resolve_ext(in,e,'prefer_ext_files',p.Results.prefer_compressed_files,varargin{:});
          %need fileparts again
          [d,~,e]=fileparts(in);
          %try to uncompress archives
          try
            switch lower(e)
            case {'.z','.zip'}
              arch=true;
              out=unzip(in,d);
            case {'.tgz','.tar.gz','.tar'}
              arch=true;
              out=untar(in,d);
              %get rid of PaxHeaders
              out(cells.cellstrin(out,'PaxHeaders'))=[];
            case {'.gz','.gzip'}
              arch=true;
              out=gunzip(in,d);
            otherwise
              arch=false;
              out=in;
            end
            if arch && p.Results.debug
              str.say('say_stack_delta',p.Results.say_stack_delta,['From archive ''',in,''' extracted the following files:'],newline,...
                strjoin(out,newline))
            end
          catch
            %if the zip file is corrupted, assume data file is missing
            str.say('say_stack_delta',p.Results.say_stack_delta,['WARNING: error extracting archive ''',in,'''.'])
            out=in;
            return
          end
          %handle zipped files
          if arch
            %some sanity
            if ~iscellstr(out)
              error(['expecting variable ''unzipped_filename'' to be a cellstr, not a ''',...
                class(out),'''.'])
            end
            if numel(out)~=1
              error(['expecting zip archive ''',out,''' to contain one file only, not ',...
                num2str(numel(out)),':',10,strjoin(out,'\n')])
            end
          end
        else
          %look for the compressed file (resolve_mat resolved according to file creation date)
          for i=file.archivedfilesext
            [out,no_change_flag]=file.resolve_ext(in,i{1},'prefer_ext_files',p.Results.prefer_compressed_files,varargin{:});
            %check if things changes
            if ~no_change_flag
              %run again this routine, so that compressed files can be expanded
              out=file.resolve_compressed(out,varargin{:});
              %we're done, no need to resolve other archive types (assumes there is only one type per expanded file, which is reasonable)
              break
            end
          end
        end
      end
    end
    function delete_uncompressed(filename,varargin)
      p=machinery.inputParser;
      p.addRequired( 'filename',       @(i) ischar(i));
      p.addParameter('del_arch', true, @(i) isscalar(i) && islogical(i))
      p.parse(filename,varargin{:})
      %delete uncompressed file if compressed file is there
      if p.Results.del_arch
        %loop over all supported compressed formats
        for i=1:numel(file.archivedfilesext)
          compressed_file=file.replace_ext(filename,file.archivedfilesext{i},varargin{:});
          if file.exist(compressed_file)
            delete(filename)
            disp(['Deleted uncompressed file ''',filename,''' because compressed file exists: ''',compressed_file,'''.'])
          end
        end
      end
    end
    %% build a list of filenames from 'start' to 'stop' with 'period' step, if file.dateplaceholder is present in 'in'
    function out=resolve_timestamps(in,varargin)
      v=varargs.wrap('sources',{....
        {...
          'start',       time.zero_date, @(i) isscalar(i) && isdatetime(i);...
          'stop',        time.inf_date,  @(i) isscalar(i) && isdatetime(i);...
          'period',      days(1),        @(i) isscalar(i) && isduration(i);...
          'date_fmt',    'yyyy-mm-dd',   @ischar;...
          'scalar_as_strings',false,     @(i) islogical(i) && isscalar(i);...
        },...
      },varargin{:});
      %reduce scalar
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        out=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_timestamps(i,varargin{:}),in,'UniformOutput',false)...
        )));
      else
        %trivial calll
        if ~contains(in,file.dateplaceholder)
          out=in;
          return
        end
        %sanity
        assert(v.start<v.stop,[...
          'input ''start'' (',datestr(v.start),') is later than ',...
          'input ''stop'' (', datestr(v.stop),').'...
        ])
        %NOTICE: it is NOT assumed that out are scalar below
        % build a cell string with the dates
        date_list=time.list(v.start,v.stop,v.period);
        out=cell(size(date_list));
        for i=1:numel(date_list)
          out{i}=strrep(in,file.dateplaceholder,datestr(date_list(i),v.date_fmt));
        end
      end
      %convert to string if requested
      out=file.ensure_scalar(out,v.scalar_as_strings);
    end
    %% retrive remote file if needed
    %NOTICE: this function will always return scalar strings
    function out=resolve_remote(in,varargin)
      p=machinery.inputParser;
      p.addRequired( 'in',                          @(i) ischar(i)   || iscellstr(i));
      p.addParameter('remote_url'       , ''      , @(i) url.is(i) || str.none(i));
      p.addParameter('expiration_period', days(30), @(i) isscalar(i) && isduration(i));
      p.addParameter('data_age'         , days(0) , @(i) isscalar(i) && isduration(i));
      p.addParameter('debug'            , false   , @islogical);
      p.parse(in,varargin{:})
      %trivial call (most of the times there's not need to retrieve remote data)
      if isempty(p.Results.remote_url)
        out=in;
        return
      end
      %reduce scalar cell to char (avoids infinite loops with vector mode)
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        out=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_remote(i,varargin{:}),in,'UniformOutput',false)...
        )));
        %convert to string if requested
        out=file.ensure_scalar(out,p.Results.scalar_as_strings);
        return
      end
      %NOTICE: it is assumed that in/out are scalar below
      %find reasons to download the data
      download_flag=false;
      %check if data age is higher than expiration period
      if p.Results.data_age > p.Results.expiration_period
        if p.Results.debug, str.say('Downloading because data_age > expiration_period'); end
        download_flag=true;
      end
      %check if all files exist
      if ~all(file.exist(in))
        if p.Results.debug, str.say('Downloading because local files are missing'); end
        download_flag=true;
      end
      %check if all files are not too old
      if ~all(file.age(in)<p.Results.expiration_period)
        if p.Results.debug, str.say('Downloading because local files are too old'); end
        download_flag=true;
      end
      %check if not downloading is explicily asked
      if str.none(p.Results.remote_url)
        download_flag=false;
      end
      %download if needed
      if download_flag
        [d,f,e]=fileparts(in);
        if isempty(d); d='./'; end
        f=[f,e];
        file.ensuredir(in);
        if file.iswildcarded(in)
          if p.Results.debug, str.say('Using wget to download wildcarded file:',in); end
          file.system(...
            ['wget ',p.Results.remote_url,'/',f],...
            'cd',d,...
            'disp',true,...
            'stop_if_error',true...
          );
          %patch output (wildcards must be resolved externally, possibly with file.resolve_wildcards)
          out=in;
        else
          if url.is_web(p.Results.remote_url)
            if p.Results.debug, str.say('Using websave to download non-wildcarded file:',in); end
            out=websave(in,[p.Results.remote_url,'/',f]);
          elseif url.is_ftp(p.Results.remote_url)
            if p.Results.debug, str.say('Using wget to download non-wildcarded file:',in); end
            file.system(['wget -O ',f,' ',p.Results.remote_url,'/',f],'disp',true,'cd',d);
            out=in;
          else
            error(['Cannot support URL ',p.Results.remote_url])
          end
        end
        assert(file.exist(out),['BUG TRAP: cannot find downloaded file: ',out])
        if p.Results.debug, str.say('Finished downloading file:',out); end
      else
        %patch output
        out=in;
        if p.Results.debug, str.say('No need to re-download file:',out); end
      end
    end
    %% resolves anything
    %NOTICE: this function only deals with resolving filenames:
    % - replaces file.dateplaceholder with list of dates (from 'start','stop','period')
    % - checks if local files are older than 'outdated_period' and downloads them from 'remote_url' (assumed to be a directory)
    % - retrieve directory contents according to wildcarded filenames ('passtrough_if_empty' controls what happens if not file found)
    % - looks for compressed files and extracts them if they are newer than the uncompressed (i.e., the input)
    % - looks for .mat files and returns them if they exist
    %NOTICE: this is done elsewhere:
    % - loading/saving in mat format (see file.load_mat and file.save_mat
    % - loading/saving in ascii/other format (each class has its own)
    function io=unwrap(io,varargin) %used to be simpletimeseries.unwrap_datafiles
      p=machinery.inputParser;
      p.addRequired( 'io',          @(i) ischar(i) || iscellstr(i));
      p.addParameter('op_sequence', {...
        @file.resolve_timestamps;...
        @file.resolve_remote;...
        @file.resolve_wildcards;...
        @file.resolve_compressed;...
        @file.resolve_ext;...
      }, @(i) isa(i,'function_handle') || (iscell(i) && all(cellfun(@(j) isa(j,'function_handle'), i))));
      p.addParameter('debug'            , false, @islogical);
      p.addParameter('scalar_as_strings', false, @(i) islogical(i) && isscalar(i));
      p.parse(io,varargin{:})
      %loop over all operations
      for i=1:numel(p.Results.op_sequence)
        f=p.Results.op_sequence{i};
        if p.Results.debug
          if ischar(io)
            str.say('Starting ',f,':',newline,io);
          elseif numel(io)<10
            str.say('Starting ',f,':',newline,strjoin(io,newline));
          else
            str.say('Starting ',f,':',newline,strjoin(io(1:5),newline),newline,strjoin(io(end-4:end),newline));
          end
        end
        switch func2str(f)
          case 'file.resolve_ext';       io=f(io,'.mat',varargin{:});
          %NOTICE: passtrough_if_empty deals with cases that the non-mat is absent but the mat file is available
          case 'file.resolve_wildcards'; io=f(io,'passtrough_if_empty',true,varargin{:});
          otherwise;                     io=f(io,varargin{:});
        end
        if p.Results.debug
          if ischar(io)
            str.say('Finished ',f,':',newline,io);
          elseif numel(io)<10
            str.say('Finished ',f,':',newline,strjoin(io,newline));
          else
            str.say('Finished ',f,':',newline,strjoin(io(1:5),newline),newline,strjoin(io(end-4:end),newline));
          end
        end
      end
      %always return cell arrays
      io=cells.scalar(io,'set');
      %... unless explicitly requested to get a scalar
      io=file.ensure_scalar(io,p.Results.scalar_as_strings);
    end
    %% CLI interfaces
    function [out,result]=system(com,varargin)
      p=machinery.inputParser;
      p.addRequired('com',                   @ischar);
      p.addParameter('cd',            '',    @ischar);
      p.addParameter('disp',          false, @islogical);
      p.addParameter('stop_if_error', false, @islogical);
      p.parse(com,varargin{:})
      if ~isempty(p.Results.cd)
        cdnow=p.Results.cd;
      else
        cdnow=pwd;
      end
      [status,result]=system(['cd ',strrep(cdnow,' ','\ '),'; ',com]);
      result=str.chomp(result);
      out=(status==0);
      if ~out
        if p.Results.stop_if_error
          error([str.dbstack,result])
        else
          disp(['WARNING: problem issuing command ''',com,''' from directory ',cdnow,newline,'error message:',newline,result])
        end
      else
        if p.Results.disp
          disp(result)
        end
      end
    end
    function [out,result]=python3(com,varargin)
      [out,result]=file.system(['echo "',com,'" | python3'],varargin{:});
    end
    function [out,s]=find(dirs,varargin)
      if iscell(dirs)
        assert(all(cellfun(@ischar,dirs)),'Need all elements of input ''dirs'' to be char.')
      elseif ischar(dirs)
        dirs={dirs};
      else
        error(['Need input ''dirs'' to be cellstr or char, not ',class(dirs)])
      end
      com=['find "',strjoin(str.clean(dirs,'regex'),'" "'),'"'];
      if numel(varargin)>0
        com=[com,' ',strjoin(str.clean(varargin,'regex'),' ')];
      end
      [s,r]=file.system(com);
      if s
        out=cells.rm_empty(strsplit(r,newline));
      else
        out={};
      end
    end
    function [out,s]=rsync(from,to,more_args,args)
      if ~exist('more_args','var')
        more_args='';
      end
      if ~exist('args','var') || isempty(args)
        args='--archive --hard-links --copy-unsafe-links --recursive --itemize-changes --exclude=._* --exclude=.*';
      end
      assert(~isempty(from) && ~isempty(to),'Inputs ''from'' and ''to'' cannot be empty.')
      [s,r]=file.system(['rsync ',args,' ',more_args,' ',from,' ',to]);
      if s
        out=cells.rm_empty(strsplit(r,newline));
        out=cellfun(@(i) strsplit(i),out,'UniformOutput',false);
        out=cellfun(@(i) i{2:end},   out,'UniformOutput',false);
      else
        out={};
      end
    end
    function out=mkdir(dirname,disp_flag)
      if ~exist('disp_flag','var') || isempty(disp_flag)
        disp_flag=false;
      end
      dirname=file.resolve_home(dirname);
      if ~exist(dirname,'dir')
        [out,msg]=mkdir(dirname);
        if disp_flag
          if out
            disp(['Created directory: ',dirname])
          else
            disp(['WARNING: could not create directory: ',dirname,':',newline,msg])
          end
        end
      else
        out=true;
        if disp_flag
          disp(['Directory already available: ',dirname])
        end
      end
    end
    function out=ln(source,target,disp_flag)
      if ~exist('disp_flag','var') || isempty(disp_flag)
        disp_flag=false;
      end
      source=file.resolve_home(source);
      assert(file.exist(source),['Cannot link to ',source,' because it does not exist'])
      [~,sn]=fileparts(source);
      target=file.resolve_home(fullfile(target,sn));
      if ~exist(target,'file')
        %avoid circular links
        if strcmp(file.fullpath(source),file.fullpath(target))
          disp(['WARNING: circular linking ',source,' to ',target,' is ilegal'])
          out=false;
        else
          out=file.system(['ln -sv ',source,' ',target],'disp',disp_flag);
        end
      else
        out=false;
        if disp_flag
          disp(['WARNING: cannot link ',source,' to ',target,'; the latter exists'])
        end
      end
    end
    function result=ls(in,ls_flags)
      if ~exist('in','var') || isempty(in)
        in=pwd;
      end
      if ~exist('ls_flags','var') || isempty(ls_flags)
        ls_flags='';
      end
      in=file.resolve_home(in);
      [~,result]=file.system(['ls ',ls_flags,' ',in,],'disp',false);
    end
    function result=md5(in,md5_flags)
      if ~exist('ls_flags','var') || isempty(md5_flags)
        md5_flags='-q';
      end
      in=file.resolve_home(in);
      [status,result]=system(['md5 ',md5_flags,' ',in,]);
      out=(status==0);
      if ~out
        disp(['WARNING: problem issuing command ''md5 ',md5_flags,' ',in,''':',newline,result])
      else
        result=result(1:end-1);
      end
    end
    function out=git(mode,filename,commit_msg)
      if ~exist('filename','var') || isempty(filename)
        filename=[mfilename('fullpath'),'.m'];
      end
      if ~exist('commit_msg','var') || isempty(commit_msg)
        commit_msg=['updated ',filename];
      end
      %check if this is a directory
      if ~exist(filename,'dir')~=0
        [d,f,e]=fileparts(filename);
        if isempty(d); d='./'; end
        f=[f,e];
      else
        %this is actually a directory
        d=filename;
        f='';
      end
      git_com=['git -C ',d];
      switch mode
        case 'date'
          git_com=[git_com,' log -1 --format=%cd --date=iso-local ',f];
        case {'status','st','add'}
          git_com=[git_com,' ',mode,' ',fullfile(d,f)];
        case {'commit','ci'}
          git_com=[git_com,' ',mode,' -m "',commit_msg,'"'];
        case 'push'
          git_com=[git_com,' ',mode];
        case 'isuptodate'
          out=str.contains(...
            file.git('status',filename),'nothing to commit, working tree clean'...
          );
          return
        otherwise
          error(['Cannot handle git command ''',mode,'''.'])
      end
      %git it
      [~,out]=file.system(git_com,'disp',false,'stop_if_error',false);
      %pre-pend directory report
      out=['Result of "',git_com,'" is:',newline,out];
    end
    %% general utils
    function io=ensure_scalar(io,scalar_as_strings)
      %convert to string if requested
      if iscellstr(io)
        if scalar_as_strings
          assert(numel(io)==1,['If ''scalar_as_strings'' is true, then results must be scalar, not with length ',...
            num2str(numel(io)),'.'])
          io=cells.scalar(io,'get');
        end
      elseif ischar(io)
        %do nothing
      else
        error(['Cannot handle inputs of class ''',class(io),'''.'])
      end
    end
    function out=ensuredir(filename,file_flag)
      if ~exist('file_flag','var') || isempty(file_flag)
        file_flag=true;
      end
      if file_flag
        d=fileparts(filename);
      else
        d=filename;
      end
      if ~isempty(d)
        out=file.mkdir(d);
      else
        out=true;
      end
    end
    function count=length(filename)
      [fid,~,close_file]=file.open(filename,'r');
      count = 0;
      while true
        if ~ischar( fgetl(fid) ); break; end
        count = count + 1;
      end
      if close_file, fclose(fid); end
    end
    function out=fullpath(filename)
      if iscellstr(filename)
        out=cellfun(@(i) file.fullpath(i),filename,'UniformOutput',false);
      	return
      end
      if contains(filename,['~',filesep])
        filename=strrep(filename,['~',filesep],[getenv('HOME'),filesep]);
      elseif numel(filename)==1 && strcmp(filename,'~')
        filename=strrep(filename,'~',[getenv('HOME'),filesep]);
      elseif contains(filename,'~')
        filename=strrep(filename,'~',[getenv('HOME'),filesep,'..',filesep]);
      end
      out = char(java.io.File(filename).getCanonicalPath());
    end
    function out=absolutepath(relativepath,root)
      if ~exist('root','var') || isempty(root)
        root=file.scriptdir;
      end
      %check if this is really a relative path
      if file.isabsolute(relativepath)
        out=file.fullpath(relativepath);
      else
        out=fullfile(root,relativepath);
      end
    end
    function out=isabsolute(filename)
      out=filename(1)==filesep || filename(1)=='~';
    end
    function io=trailing_filesep(io)
      if iscellstr(io)
        isdir=cellfun(@(i) exist(i,'dir')~=0,io);
        io(isdir)=cellfun(@(i) file.trailing_filesep(i),io(isdir),'UniformOutput',false);
      	return
      end
      if io(end)~=filesep
        io=[io,filesep];
      end
    end
    %NOTICE: if 'in' is a dangling link, out is false
    function out=exist(in,force,disp_flag)
      if ~exist('force','var') || isempty(force)
        force=false;
      end
      if ~exist('disp_flag','var') || isempty(disp_flag)
        disp_flag=false;
      end
      %vector mode
      if iscellstr(in)
        out=cellfun(@(i)file.exist(i,force,disp_flag),in);
      elseif ischar(in)
        %NOTICE: exist returns 7 if in 1st input is a dir and 2nd input is 'file' or 'dir'
        %        but returns 2 if 1st input is a dir only if 2nd input is 'file' (returns 0 if it is 'dir')
        %        so using 'file' captures both cases
        out= ~str.logical(force) && exist(file.resolve_home(in),'file')~=0;
        if disp_flag
          if exist(file.resolve_home(in),'file')~=0
            if str.logical(force)
              str.say('say_stack_delta',1,'file exists but force is true',in)
            else
              str.say('say_stack_delta',1,'file exists',in)
            end
          else
            str.say('say_stack_delta',1,'file does not exist',in)
          end
        end
      else
        error(['Cannot handle inputs of class ',class(in),'.'])
      end
    end
    function out=datenum(in)
      %reduce scalar cell to char (avoids infinite loops with vector mode)
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        out=arrayfun(@(i) file.datenum(i),in);
        return
      end
      fileinfo=dir(file.resolve_home(in));
      if exist(in,'dir')
        assert(fileinfo(1).name=='.',['Expecting the first entry of ''fileinfo'' to be relative to ''.'', not to ''',fileinfo(1).name,'''.'])
        out=fileinfo(1).datenum;
      elseif isempty(fileinfo)
        out=nan;
      else
        assert(isscalar(fileinfo),['Expecting ''fileinfo'' to be scalar, since ''',in,''' doesn not appears to be a directory.'])
        out=fileinfo.datenum;
      end
    end
    function out=datetime(in)
      out=datetime(file.datenum(in),'ConvertFrom','datenum');
    end
    function out=datestr(in)
      out=datestr(file.datetime(in));
    end
    function [newest_file,newest_date]=newest(in,varargin)
      newest_date=0;
      if ~iscellstr(in)
        in=file.unwrap(in,varargin{:});
      end
      newest_file='';
      for i=1:numel(in)
        if file.datenum(in{i})>newest_date
          newest_file=in{i};
          newest_date=file.datenum(newest_file);
        end
      end
    end
    function out=age(in)
      out=datetime('now')-file.datetime(in);
    end
    function out=up(path,n)
      if ~exist('n','var') || isempty(n)
        n=1;
      end
      path=file.resolve_wildcards(path,'disp',false,'scalar_as_strings',true);
      switch exist(path,'dir')
      case 2
        path=fileparts(path);
      case 7
        %do nothing
      otherwise
        error(['Cannot find directory ''',path,'''.'])
      end
      cdnow=pwd;
      cd(path)
      for i=1:n
        cd ..
      end
      out=pwd;
      cd(cdnow);
    end
    function out=build(varargin)
      if nargin>1
        ext=cells.scalar(varargin{end},'get');
        %if the last argument is char, then assume it is an extension
        %(force using '.' as separator)
        if ischar(ext)
          out=strjoin(file.build_particle(varargin{1:end-1}),file.build_element_char);
          out=[out,'.',cells.scalar(varargin{end},'get')];
        else
          out=strjoin(file.build_particle(varargin{:}),file.build_element_char);
        end
      else
        out=cells.scalar(varargin{1},'get');
      end
      %make sure filename is not too long
      if length(out)>512
        %get the elements
        elements=file.build_particle(varargin{:});
        %get length of elements
        element_len=cellfun(@(i) numel(i),elements);
        %find the element that is too long
        idx=find(element_len==max(element_len));
        %get 10% of the longest element length
        trim_len=floor(element_len(idx)*0.1);
        %trim down longest element
        elements{idx}=[...
          elements{idx}(1:trim_len),...
          '...',...
          elements{idx}(end-trim_len+1:end)...
        ];
        %try again
        out=file.build(elements{:});
      end
      %resolve multi-machine paths
      out=file.resolve_home(out);
    end
    function out=build_particle(varargin)
      out=cellfun(...
        @(i) str.show(...
          cells.rm_empty(cells.flatten(cells.scalar(i,'set'))),'join_char',file.build_particle_char...
        ),...
        cells.rm_empty(varargin),'UniformOutput',false...
      );
    end
    function io=remove_ext(io)
      [io,changed]=cells.scalar(io,'set');
      for i=1:numel(io)
        [p,f]=fileparts(io{i});
        io{i}=fullfile(p,f);
      end
      if changed
        io=cells.scalar(io,'get');
      end
    end
    function out=basename(in)
      [~,f,e]=fileparts(in);
      out=[f,e];
    end
    function out=hostname
      out=file.system('hostname');
    end
  end
end

function out = isnumstr(in)
  %Determines if the input string holds only numerical data. The test that is
  %made assumes that numerical data includes only the characters defined in
  %<num_chars>.
  if isnumeric(in)
      out=true;
      return
  end
  %characters that are allowed in numerical strings
  num_chars = ['1234567890.-+EeDd:',9,10,13,32];
  %deleting the numerical chars from <in>
  for i=1:length(num_chars)
      num_idx = strfind(in,num_chars(i));
      if ~isempty(num_idx)
          %flaging this character
          in(num_idx)='';
      end
  end
  %now <in> must be empty (if only numerical data)
  out = isempty(in);
  %warning
  %if ~out
  %    disp(['debug: non-numerical chars: ',in])
  %end
end

function out = reduce_cell_array(in)
  out=cell(size(in));
  for i=1:numel(in)
      if ~isempty(in{i})
          out{i}=in{i}{1};
      else
          out{i}='';
      end
  end
  % clean up empty entries
  out(cellfun('isempty',out))=[];
end