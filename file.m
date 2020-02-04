classdef file
  properties(Constant)
    %NOTICE: this used to be called 'DATE_PLACE_HOLDER'
    dateplaceholder='DATE_PLACEHOLDER';
    archivedfilesext={'.gz','.gzip','.z','.zip','.tgz','.tar.gz','.tar'};
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
    %as this script; typical dirs are: aux, packages, metadata, data, plot (the latter 2 are not in git)
    function out=orbdir(type)
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
    function [fid,filename,close_file]=open(filename,perm)
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
      else
        %make sure dir exists
        d=fileparts(filename);
        if isempty(strfind(perm,'r'))
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
        if isempty(header{i})
          error([mfilename,': error reading header of file ',filename,'. Error message:',10,ferror(fid)])
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
          error([mfilename,...
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
        error([mfilename,': need one input file.'])
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
          str.say('stack_delta',1,['NOTICE: removed ',num2str(size(data{i},1)-min(len)),' entries from the ',...
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
        error([mfilename,': need filename.'])
      end
      if ~exist('str','var') || isempty(str)
        error([mfilename,': need str.'])
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
        error([mfilename,': need filename.'])
      end
      %open the file if needed
      [fid,~,close_file]=file.open(filename,'r');
      %save the string
      str=fscanf( fid, '%s');
      %close the file (if fid not given)
      if close_file, fclose(fid); end
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
        switch lower(direction)
        case 'set'; if ~file.isext(io,ext); io=[io,ext]; end
        case 'get'; if  file.isext(io,ext); io=io(1:end-length(ext));          end
        otherwise; error(['Cannot handle input ''direction'' with value ''',direction,'''.'])
        end
      end
    end
    % prefers mat or non-mat file in a file list (that may mix mat with non-mat files)
    % a mat file has extension .mat, a non-mat file has a different extension (strangely enough)
    % this means that if in={'file1.dat','file1.dat.gz'} and file1.dat.mat exists, this function will
    % return {'file1.dat.mat','file1.dat.gz'} because it looks for file1.dat.gz.mat as the corresponding mat file.
    function [filenames,no_change_flag]=resolve_ext(in,ext,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'in',                            @(i) ischar(i)   || iscellstr(i)); 
      p.addRequired( 'ext',                           @(i) ischar(i)                  ); 
      p.addParameter('prefer_non_ext_files',   false, @(i) islogical(i) && isscalar(i));
      p.addParameter('prefer_ext_files',       false, @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',      false, @(i) islogical(i) && isscalar(i));
      p.addParameter('debug',                  false, @islogical);
      p.parse(in,ext,varargin{:})
      %sanity
      assert(~(p.Results.prefer_non_ext_files && p.Results.prefer_ext_files),'Cannot both prefer non-ext and ext files')
      if iscellstr(in)
        filenames=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_ext(i,ext,varargin{:}),in,'UniformOutput', false)...
        )));
        %convert to string if requested
        filenames=file.ensure_scalar(filenames,p.Results.scalar_as_strings);
      else
        %first get non-ext file
        nonext=file.ext(in,ext,'get'); nonextexists=file.exist(nonext);
        %get ext file name
        yesext=file.ext(nonext,ext,'set'); yesextexists=file.exist(yesext);
        %if there's only the ext file, then preference is irrelevant
        if ~nonextexists && yesextexists
          filenames=yesext;
          if p.Results.debug; str.say('picked yesext:',yesext,'(nonext file non-existing:',nonext,')'); end
        %if there's only the non-ext file, then preference is irrelevant
        elseif nonextexists && ~yesextexists
          filenames=nonext;
          if p.Results.debug; str.say('picked nonext:',nonext,'(yesext file non-existing:',yesext,')'); end
        %if there's a ext file and that is prefered, take it
        elseif yesextexists&& p.Results.prefer_ext_files
          filenames=yesext;
          if p.Results.debug; str.say('picked yesext:',yesext,'(preferred)'); end
        %if there's a non-ext file and that is prefered, take it
        elseif nonextexists&& p.Results.prefer_non_ext_files
          filenames=nonext;
          if p.Results.debug; str.say('picked nonext:',nonext,'(preferred)'); end
        %if there's no preferences and both files exists, return the newest
        else
          if p.Results.debug
            str.say('mat:',datetime(file.datenum(yesext),'convertfrom','datenum'),yesext   )
            str.say('out:',datetime(file.datenum(nonext),'convertfrom','datenum'),nonext)
          end
          if file.datenum(yesext)>file.datenum(nonext)
            filenames=yesext;
            if p.Results.debug; str.say('picked yesext:',yesext,'(newer)'); end
          else
            filenames=nonext;
            if p.Results.debug; str.say('picked nonext:',nonext,'(newer)'); end
          end
        end
      end
      %handle additional outputs
      if nargout>1
        no_change_flag=str.iseq(filenames,in);
      end
    end
    %% resolves wildcarded filenames
    function [filenames,wildcarded_flag]=resolve_wildcards(in,varargin)
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'in',                         @(i) ischar(i)   || iscellstr(i)); 
      p.addParameter('disp',                true,  @(i) islogical(i) && isscalar(i));
      p.addParameter('directories_only',    false, @(i) islogical(i) && isscalar(i));
      p.addParameter('files_only',          false, @(i) islogical(i) && isscalar(i));
      p.addParameter('stop_if_empty',       false, @(i) islogical(i) && isscalar(i));
      p.addParameter('passtrough_if_empty', false, @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',   false, @(i) islogical(i) && isscalar(i));
      p.addParameter('stack_delta',find(arrayfun(@(i) strcmp(i.file,'file.m'),dbstack),1,'last'), @(i) isnumeric(i) && isscalar(i));
      % parse it
      p.parse(in,varargin{:});
      %reduce scalar cell to char (avoids infinite loops with vector mode)
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        filenames=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_wildcards(i,varargin{:}),in,'UniformOutput',false)...
        )));
        %convert to string if requested
        filenames=file.ensure_scalar(filenames,p.Results.scalar_as_strings);
        return
      end
      %trivial call
      if isempty(in)
        if p.Results.scalar_as_strings
          filenames='';
        else
          filenames={''};
        end
        return
      end
      %wildcard character '*' needs to be translated to '.*' (if not already)
      [in,wildcarded_flag]=translate_wildcard(file.resolve_home(in));
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
          filenames=a;
        else
          filenames={a};
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
          filenames=fullfile(a,f_in);
          return
        end
        if p.Results.disp
          str.say('stack_delta',p.Results.stack_delta,['WARNING: ',msg])
        end
        filenames='';
        return
      end
      %add path
      root=fileparts(in);
      filenames=cellfun(@(i) fullfile(root,i),f1,'UniformOutput',false);
      %inform user
      if p.Results.disp && wildcarded_flag
        str.say('stack_delta',p.Results.stack_delta,['found ',num2str(numel(filenames)),' filenames in wildcarded string ''',in,'''.'])
      end
      %convert to string if requested
      if p.Results.scalar_as_strings
        assert(numel(filenames)==1,['If ''scalar_as_strings'' is true, then results must be scalar, not with length ',...
          num2str(numel(filenames)),'.'])
        filenames=filenames{1};
      end
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
    function filenames=resolve_compressed(in,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'in',                            @(i) ischar(i)   || iscellstr(i)); 
      p.addParameter('prefer_compressed_files',false, @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',      false, @(i) islogical(i) && isscalar(i));
      p.addParameter('stack_delta',find(arrayfun(@(i) strcmp(i.file,'file.m'),dbstack),1,'last'), @(i) isnumeric(i) && isscalar(i));
      p.parse(in,varargin{:})
      %reduce scalar
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        filenames=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_compressed(i,varargin{:},'scalar_as_strings',true),in,'UniformOutput',false)...
        )));
        %convert to string if requested
        filenames=file.ensure_scalar(filenames,p.Results.scalar_as_strings);
        return
      end
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
            filenames=unzip(in,d);
          case {'.tgz','.tar.gz','.tar'}
            arch=true;
            filenames=untar(in,d);
            %get rid of PaxHeaders
            filenames(cells.cellstrin(filenames,'PaxHeaders'))=[];
          case {'.gz','.gzip'}
            arch=true;
            filenames=gunzip(in,d);
          otherwise
            arch=false;
            filenames=in;
          end
          if arch
            str.say('stack_delta',p.Results.stack_delta,['From archive ''',in,''' extracted the following files:'],newline,...
              strjoin(filenames,newline))
          end
        catch
          %if the zip file is corrupted, assume data file is missing
          str.say('stack_delta',p.Results.stack_delta,['WARNING: error extracting archive ''',in,'''.'])
          return
        end
        %handle zipped files
        if arch
          %some sanity
          if ~iscellstr(filenames)
            error([mfilename,': expecting variable ''unzipped_filename'' to be a cellstr, not a ''',...
              class(filenames),'''.'])
          end
          if numel(filenames)~=1
            error([mfilename,': expecting zip archive ''',filenames,''' to contain one file only, not ',...
              num2str(numel(filenames)),':',10,strjoin(filenames,'\n')])
          end
        end
      else
        %look for the compressed file (resolve_mat resolved according to file creation date)
        for i=file.archivedfilesext
          [filenames,no_change_flag]=file.resolve_ext(in,i{1},'prefer_ext_files',p.Results.prefer_compressed_files,varargin{:});
          %check if things changes
          if ~no_change_flag
            %run again this routine, so that compressed files can be expanded
            filenames=file.resolve_compressed(filenames,varargin{:});
            %we're done, no need to resolve other archive types (assumes there is only one type per expanded file, which is reasonable)
            break
          end
        end
      end
    end
    %% resolves timestamps
    function filenames=resolve_timestamps(in,varargin)
      p=inputParser; p.KeepUnmatched=true;
      p.addRequired( 'in',                                         @(i) ischar(i)   || iscellstr(i)); 
      p.addParameter( 'start',       time.ToDateTime(0,'datenum'), @(i) isscalar(i) && isdatetime(i));
      p.addParameter( 'stop',        time.ToDateTime(0,'datenum'), @(i) isscalar(i) && isdatetime(i));
      p.addParameter( 'period',      days(1),                      @(i) isscalar(i) && isduration(i));
      p.addParameter( 'date_fmt',    'yyyy-mm-dd',                 @(i) ischar(i));
      p.addParameter( 'only_existing',true,                        @(i) islogical(i));
      p.addParameter( 'scalar_as_strings',false,                   @(i) islogical(i) && isscalar(i));
      p.parse(in,varargin{:})
      %reduce scalar
      in=cells.scalar(in);
      %vector mode
      if iscellstr(in)
        filenames=unique(cells.rm_empty(cells.flatten(...
          cellfun(@(i) file.resolve_timestamps(i,varargin{:}),in,'UniformOutput',false)...
        )));
        %convert to string if requested
        filenames=file.ensure_scalar(filenames,p.Results.scalar_as_strings);
        return
      end
      %trivial calll
      if isempty(strfind(in,file.dateplaceholder))
        filenames=in;
        return
      end
      %sanity
      assert(p.Results.start<p.Results.stop,[...
        'input ''start'' (',datestr(p.Results.start),') is later than ',...
        'input ''stop'' (',datestr(p.Results.stop),').'...
      ])
      % build a cell string with the dates
      date_list=time.list(p.Results.start,p.Results.stop,p.Results.period);
      filenames=cell(size(date_list));
      for i=1:numel(date_list)
        filenames{i}=strrep(in,file.dateplaceholder,datestr(date_list(i),p.Results.date_fmt));
      end
    end
    %% resolves anything
    function io=unwrap(io,varargin) %used to be simpletimeseries.unwrap_datafiles
      p=inputParser;
      p.KeepUnmatched=true;
      p.addRequired( 'io',          @(i) ischar(i) || iscellstr(i));
      p.addParameter('op_sequence', {...
        @file.resolve_timestamps;...
        @file.resolve_wildcards;...
        @file.resolve_compressed;...
        @file.resolve_ext;...
      },                            @(i) isfun(i) || iscell(i) && all(cellfun(@(j) isa(j,'function_handle'), i))); 
      p.addParameter('debug', false,  @islogical);
      p.parse(io,varargin{:})
      %loop over all operations
      for i=1:numel(p.Results.op_sequence)
        f=p.Results.op_sequence{i};
        if p.Results.debug;
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
        if p.Results.debug;
          if ischar(io)
            str.say('Finished ',f,':',newline,io);
          elseif numel(io)<10
            str.say('Finished ',f,':',newline,strjoin(io,newline));
          else
            str.say('Finished ',f,':',newline,strjoin(io(1:5),newline),newline,strjoin(io(end-4:end),newline));
          end
        end
      end
    end
    %% CLI interfaces
    function [out,result]=system(com,disp_flag,stop_if_error)
      if ~exist('disp_flag','var') || isempty(disp_flag)
        disp_flag=false;
      end
      if ~exist('stop_if_error','var') || isempty(stop_if_error)
        stop_if_error=false;
      end
      [status,result]=system(com);
      result=str.chomp(result);
      out=(status==0);
      if ~out 
        if stop_if_error
          error([str.dbstack,result])
        else
          disp(['WARNING: problem issuing command ''',com,''':',newline,result])
        end
      end
      if disp_flag
        disp(result)
      end
    end
    function [out,s]=find(varargin)
      com=['find ',strjoin(str.clean(varargin,'regex'),' ')];
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
        out=mkdir(dirname,disp_flag,true);
      else
        out=true;
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
          out=file.system(['ln -sv ',source,' ',target],disp_flag);
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
      [~,result]=file.system(['ls ',ls_flags,' ',in,],false);
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
    %% general utils
    function io=ensure_scalar(io,force)
      %convert to string if requested
      if iscellstr(io)
        if force
          assert(numel(io)==1,['If ''scalar_as_strings'' is true, then results must be scalar, not with length ',...
            num2str(numel(io)),'.'])
          io=io{1};
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
        out=relativepath;
      else
        out=fullfile(root,relativepath);
      end
    end
    function out=isabsolute(filename)
      out=filename(1)==filesep;
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
    function out=exist(in)
      %vector mode
      if iscellstr(in)
        out=cellfun(@file.exist,in);
      elseif ischar(in)
        %NOTICE: exist returns 7 if in 1st input is a dir and 2nd input is 'file' or 'dir'
        %        but returns 2 if 1st input is a dir only if 2nd input is 'file' (returns 0 if it is 'dir')
        %        so using 'file' captures both cases
        out=exist(file.resolve_home(in),'file')~=0;
      else
        error(['Cannot handle inputs of class ',class(in),'.'])
      end
    end
    function out=datenum(in)
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
        in=cells.scalar(file.unwrap(in,varargin{:}),'set');
      end
      newest_file='';
      for i=1:numel(in)
        if file.datenum(in{i})>newest_date
          newest_file=in{i};
          newest_date=file.datenum(newest_file);
        end
      end
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
          cells.rm_empty(cells.flatten(cells.scalar(i,'set'))),'',file.build_particle_char...
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
  %    disp([mfilename,':debug: non-numerical chars: ',in])
  %end
end
function [io,wildcarded_flag] = translate_wildcard(io)
  wildcarded_flag=contains(io,'*');
  if wildcarded_flag
    io=str.rep(io,'.*','<dotasterisc>','*','.*','<dotasterisc>','.*');
%     idx=strfind(io,'*');
%     for i=1:numel(idx)
%       if idx(i) > 1 && ~strcmp(io(idx(i)-1),'.')
%         io=translate_wildcard([io(1:idx(i)-1),'.',io(idx(i):end)]);
%         return
%       end
%       if idx(i) == 1
%         io=translate_wildcard(['.',io(idx(i):end)]);
%         return
%       end
%     end
  end
  %handle special characters, reduce them to literals
  special_chars={'+'};
  for i=1:numel(special_chars)
    if ~isempty(strfind(io,special_chars{i}))
      io=str.rep(io,special_chars{i},['\',special_chars{i}]);
    end
  end
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

%% obsolete functions

function filenames=unwrap_old(in,varargin) %used to be simpletimeseries.unwrap_datafiles
  p=inputParser;
  p.KeepUnmatched=true;
  p.addRequired(  'in',                                        @(i) ischar(i) || iscellstr(i));
  p.addParameter( 'start',       time.ToDateTime(0,'datenum'), @(i) isscalar(i) && isdatetime(i));
  p.addParameter( 'stop',        time.ToDateTime(0,'datenum'), @(i) isscalar(i) && isdatetime(i));
  p.addParameter( 'period',      days(1),                      @(i) isscalar(i) && isduration(i));
  p.addParameter( 'date_fmt',    'yyyy-mm-dd',                 @(i) ischar(i));
  p.addParameter( 'only_existing',true,                        @(i) islogical(i));
  p.addParameter( 'stack_delta',find(arrayfun(@(i) strcmp(i.file,'file.m'),dbstack),1,'last'), @(i) isnumeric(i) && isscalar(i));
  p.parse(in,varargin{:})
  %parse wildcards and handle .mat files
  in=file.wildcard(in,varargin{:});
  %reduce scalar
  in=cells.scalar(in);
  %vector mode
  if iscellstr(in)
    filenames=cells.rm_empty(cells.flatten(...
      cellfun(@(i) file.unwrap(i,varargin{:},'scalar_as_strings',true),in,'UniformOutput',false)...
    ));
    return
  end
  %need fileparts
  [d,f,e]=fileparts(in);
  %if argument filename has a date place holder, build a cell string with those dates
  if ~isempty(strfind(in,file.dateplaceholder))
    %sanity
    if p.Results.start>=p.Results.stop
      error([mfilename,': input ''start'' (',datestr(p.Results.start),...
        ') is not after input ''stop'' (',datestr(p.Results.stop),').'])
    end
    date_list=time.list(p.Results.start,p.Results.stop,p.Results.period);
    filenames=cell(size(date_list));
    for i=1:numel(date_list)
      filenames{i}=strrep(in,file.dateplaceholder,datestr(date_list(i),p.Results.date_fmt));
    end
    %unwrap these files again
    filenames=file.unwrap(filenames,varargin{:});
    %done
    return
  end
  %now handling compressed files: prepend tar to extension if there
  if strcmp(f(max([1,end-3]):end),'.tar')
    e=['.tar',e];
  end
  %try to uncompress archives
  try
    switch lower(e)
    case {'.z','.zip'}
      arch=true;
      filenames=file.resolve_home(unzip(in,file.unresolve_home(d)));
    case {'.tgz','.tar.gz','.tar'}
      arch=true;
      filenames=file.resolve_home(untar(in,file.unresolve_home(d)));
      %get rid of PaxHeaders
      filenames(cells.cellstrin(filenames,'PaxHeaders'))=[];
    case {'.gz','.gzip'}
      arch=true;
      filenames=file.resolve_home(gunzip(in,file.unresolve_home(d)));
    otherwise
      arch=false;  
    end
    if arch
      str.say('stack_delta',p.Results.stack_delta,['Extracted archive ''',in,'''.'])
    end
  catch
    %if the zip file is corrupted, assume data file is missing
    str.say('stack_delta',p.Results.stack_delta,['WARNING: error extracting archive ''',in,'''.'])
    return
  end
  %handle zipped files
  if arch
    %some sanity
    if ~iscellstr(filenames)
      error([mfilename,': expecting variable ''unzipped_filename'' to be a cellstr, not a ''',...
        class(filenames),'''.'])
    end
    if numel(filenames)~=1
      error([mfilename,': expecting zip archive ''',filenames,''' to contain one file only, not ',...
        num2str(numel(filenames)),':',10,strjoin(filenames,'\n')])
    end
    %and we're done (no recursive unwrapping!)
    return
  end
  %if none of the conditions above were met, this is the name of a single file (return char!)
  if exist(in,'file')==0 && p.Results.only_existing
    filenames='';
  else
    filenames=in;
  end
end

%handling wildcards and .mat files
function [filenames,wildcarded_flag]=wildcard_old(in,varargin)
  % FILENAMES = FILENAME_WILDCARD(STR) looks for files which fit the wilcard
  % input string STR and returns a cell array FILENAMES with the
  % corresponding filenames.
  %
  % Optional 'strict_mat_search' as true (default is false) makes script look for mat files named:
  %
  % ['^',f_in,'\.mat$']
  %
  % Instead of the standard way which is:
  %
  % ['^',f_in,'.*\.mat$']
  %
  % Optional 'prefer_non_mat_files' as true will skip looking for corresponding .mat
  % files (default is false).
  %
  % Optional 'scalar_as_strings' returns strings (instead of cell arrays) when the
  % result is scalar (default is false). In case the result is non-scalar and this
  % parameter is true, an error is triggered.

  % Created by J.Encarnacao <J.deTeixeiradaEncarnacao@tudelft.nl>
  p=inputParser;
  p.KeepUnmatched=true;
  p.addRequired( 'in',                         @(i) ischar(i)   || iscellstr(i)); 
  p.addParameter('strict_mat_search'   ,false, @(i) islogical(i) && isscalar(i));
  p.addParameter('prefer_non_mat_files',false, @(i) islogical(i) && isscalar(i));
  p.addParameter('disp',                true,  @(i) islogical(i) && isscalar(i));
  p.addParameter('scalar_as_strings',   false, @(i) islogical(i) && isscalar(i));
  p.addParameter('directories_only',    false, @(i) islogical(i) && isscalar(i));
  p.addParameter('files_only',          false, @(i) islogical(i) && isscalar(i));
  p.addParameter('stop_if_empty',       false, @(i) islogical(i) && isscalar(i));
  p.addParameter('stack_delta',find(arrayfun(@(i) strcmp(i.file,'file.m'),dbstack),1,'last'), @(i) isnumeric(i) && isscalar(i));
  % parse it
  p.parse(in,varargin{:});
  %reduce scalar cell to char (avoids infinite loops with vector mode)
  in=cells.scalar(in);
  %vector mode
  if iscellstr(in)
    filenames=cells.rm_empty(cells.flatten(...
      cellfun(@(i) file.wildcard(i,varargin{:}),in,'UniformOutput',false)...
    ));
    return
  end
  %trivial call
  if isempty(in)
    if p.Results.scalar_as_strings
      filenames='';
    else
      filenames={''};
    end
    return
  end
  %resolve home dirs
  in=file.resolve_home(in);
  %turn off advanced features: trailing wildcard messes up big time the .mat and zipped files handling
  simple_file_search=(in(end)=='*');
  %wildcard character '*' needs to be translated to '.*' (if not already)
  [in,wildcarded_flag]=translate_wildcard(in);
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
      filenames=a;
    else
      filenames={a};
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
  %check if advanced features are available
  if simple_file_search
    f1c=f1;
  else
    %taking into account also mat files, find any relevant .mat files in current dir
    if p.Results.strict_mat_search
      f2=reduce_cell_array(regexp(f(:),['^',f_in,'\.mat$'],'match'));
    else
      f2=reduce_cell_array(regexp(f(:),['^',f_in,'.*\.mat$'],'match'));
    end
    %strip mat extension from f2 list, these are the relevant .mat files that may be in f1
    f2c_mat=cellfun(@(i) file.ext(i,'.mat','get'),f2,'UniformOutput',false);
    %delete mat files that are in also in f1
    f1c=setdiff(f1,f2c_mat);
    %add these mat files to the final file list
    f1c=unique([f1c(:);f2(:)]);
    %handle zipped files (pretty much as is done for mat files, but prefer what is in f1c)
    for i=1:numel(file.archivedfilesext)
      %find any relevant zipped files in current dir
      f2=reduce_cell_array(regexp(f(:),['^',f_in,strrep(file.archivedfilesext{i},'.','\.'),'$'],'match'));
      %avoiding the rest of the loop if there are not archived files with this ext
      if isempty(f2); continue; end
      %strip zip extension from f2 list, these are the relevant zipped files that may be in f1
      f2c=cellfun(@(j) strrep(j,file.archivedfilesext{i},''),f2,'UniformOutput',false);
      %clean this list of files of those that are in f1
      tmp=setdiff(f2c,f1);
      %restore the zip extension
      tmp=cellfun(@(j) [j,file.archivedfilesext{i}],tmp,'UniformOutput',false);
      %add these zipped files to the final file list
      f1c=unique([tmp(:);f1c(:)]);
    end
  end
  %one filename is frequently empty
  f1c=file.resolve_ext(cells.rm_empty(f1c),'.mat',p.Results.prefer_non_mat_files);
  %inform user if no files found
  if isempty(f1c)
    msg=['found no files named ''',in,'''.'];
    if p.Results.stop_if_empty
      error(msg)
    end
    if p.Results.disp
      str.say('stack_delta',p.Results.stack_delta,['WARNING: ',msg])
    end
    filenames='';
    return
  end
  %add path
  root=fileparts(in);
  filenames=cell(1,length(f1c));
  for i=1:length(f1c)
    filenames{i}=fullfile(root,f1c{i});
  end
  %inform user
  if p.Results.disp && wildcarded_flag
    str.say('stack_delta',p.Results.stack_delta,['found ',num2str(numel(filenames)),' filenames in wildcarded string ''',in,'''.'])
  end
  %convert to string if requested
  if p.Results.scalar_as_strings
    assert(numel(filenames)==1,['If ''scalar_as_strings'' is true, then results must be scalar, not with length ',...
      num2str(numel(filenames)),'.'])
    filenames=filenames{1};
  end
end