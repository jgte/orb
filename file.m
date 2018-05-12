classdef file
  properties(Constant)
    %NOTICE: this used to be called 'DATE_PLACE_HOLDER'
    dateplaceholder='DATE_PLACEHOLDER';
    archivedfilesext={'.gz','.gzip','.z','.zip','.tgz','.tar.gz','.tar'};
  end
  methods(Static)
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
        if ~exist(d,'dir') && isempty(strfind(perm,'r'))
          [st,msg]=mkdir(d);
          assert(st,['error creating directory ''',d,''': ',msg])
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
        %maybe this file has some columns with non-numeric data, so checking for changes in the number of columns
        for i=max_header_len-1:-1:1
          if nr_columns(i)~=nr_columns(end)
            nlines = i;
            break
          end
        end
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
    %checks if is mat file
    function out=ismat(in)
      if iscellstr(in)
        out=cellfun(@file.ismat,in);
      else
        [~,~,e]=fileparts(in);
        out=strcmp(e,'.mat');
      end
    end
    %add (remove) '.mat' to (from) a file name unless it already has it (not)
    function io=mat(io,direction)
      if ~exist('direction','var') || isempty(direction)
        direction='get';
      end
      if iscellstr(io)
        io=cellfun(@(i) file.mat(i,direction),io,'UniformOutput', false);
      else
        switch lower(direction)
        case 'set'; if ~file.ismat(io); io=[io,'.mat']; end
        case 'get'; if  file.ismat(io); io=io(1:end-4);          end
        otherwise; error(['Cannot handle input ''direction'' with value ''',direction,'''.'])
        end
      end
    end
    % prefers mat or non-mat file in a file list (that may mix mat with non-mat files)
    function out=mat_preference(in,prefer_non_mat_files)
      if ~exist('prefer_non_mat_files','var') || isempty(prefer_non_mat_files)
        prefer_non_mat_files=false;
      end
      %first get list of non-mat files
      out=unique(file.mat(in,'get'));
      %loop over non-mat files
      for i=1:numel(out)
        %get indexes of original files with and without mat extension
        idx=[cells.strcmp(in,out{i}),cells.strcmp(in,file.mat(out{i},'set'))];
        %if nothing was found, then there are only mat files, so preference is irrelevant
        if isempty(idx);
          % add the mat extension to this entry
          out{i}=file.mat(out{i},'set');
        else
          %check if there is any .mat file in this subset
          if any(file.ismat(in(idx))) && ~prefer_non_mat_files
            %if so, then can prefer the mat file, so add the mat extension to this entry
            out{i}=file.mat(out{i},'set');
          end
        end
      end
    end
    %handling wildcards and .mat files
    function [filenames,wildcarded_flag]=wildcard(in,varargin)
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
        f2c_mat=cellfun(@(i) file.mat(i,'get'),f2,'UniformOutput',false);
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
      f1c=file.mat_preference(cells.rm_empty(f1c),p.Results.prefer_non_mat_files);
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
    %handling file.dateplaceholder and compressed files
    %NOTICE: compressed files can only be properly handled if a trailing * is part of 'in'
    function filenames=unwrap(in,varargin) %used to be simpletimeseries.unwrap_datafiles
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
          cellfun(@(i) file.unwrap(i,varargin{:}),in,'UniformOutput',false)...
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
    function out=ensuredir(filename,file_flag)
      if file_flag
        d=fileparts(filename);
      else
        d=filename;
      end
      if ~exist(d,'dir')
        out=mkdir(d);
      else
        out=true;
      end
    end
    function [out,s]=find(varargin)
      com=['find ',strjoin(str.clean(varargin,'regex'),' ')];
      str.say('stack_delta',1,com)
      [s,r]=system(com);
      if s~=0
        str.say('stack_delta',1,r)
        out={};
      else
        out=cells.rm_empty(strsplit(r,char(10)));
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
      [s,r]=system(['rsync ',args,' ',more_args,' ',from,' ',to]);
      if s~=0
        out={};
      else
        out=cells.rm_empty(strsplit(r,char(10)));
        out=cellfun(@(i) strsplit(i),out,'UniformOutput',false);
        out=cellfun(@(i) i{2:end},   out,'UniformOutput',false);
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
      if ~isempty(strfind(filename,['~',filesep]))
        filename=strrep(filename,['~',filesep],[getenv('HOME'),filesep]);
      elseif numel(filename)==1 && strcmp(filename,'~')
        filename=strrep(filename,'~',[getenv('HOME'),filesep]);
      elseif ~isempty(strfind(filename,'~'))
        filename=strrep(filename,'~',[getenv('HOME'),filesep,'..',filesep]);
      end
      out = char(java.io.File(filename).getCanonicalPath());
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
    function io=exist(io)
      io=str.rep(io,...
        '/Users/teixeira/cloud/Work/projects','~/projects',...
        '/Users/teixeira','~',...
        '/home1/00767/byaa676','~');
      io=~isempty(file.unwrap(io,...
        'disp',             false,...
        'directories_only', false,...
        'files_only',       false...
      ));
    end
    function out=datenum(in)
      fileinfo=dir(in);
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
      path=file.wildcard(path,'disp',false,'scalar_as_strings',true);
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
  wildcarded_flag=~isempty(strfind(io,'*'));
  if wildcarded_flag
    idx=strfind(io,'*');
    for i=1:numel(idx)
      if idx(i) > 1 && ~strcmp(io(idx(i)-1),'.')
        io=translate_wildcard([io(1:idx(i)-1),'.',io(idx(i):end)]);
        return
      end
      if idx(i) == 1
        io=translate_wildcard(['.',io(idx(i):end)]);
        return
      end
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
