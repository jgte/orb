classdef file
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
          disp(['NOTICE: removed ',num2str(size(data{i},1)-min(len)),' entries from the ',...
            num2str(i),'-th column of the data in file ',filename,'.'])
          data{i} = data{i}(1:min(len),:);
        end
      end
      %numerical output, so transforming into numerical array
      data=[data{cellfun(@(i) isnumeric(i),data)}];
      %close the file (if fid not given)
      if close_file, fclose(fid); end
    end
    function out=wildcard(in,varargin)
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
      p.addRequired( 'in',@(i) ischar(i));
      p.addParameter('strict_mat_search'   ,false, @(i) islogical(i) && isscalar(i));
      p.addParameter('prefer_non_mat_files',false, @(i) islogical(i) && isscalar(i));
      p.addParameter('disp',                true,  @(i) islogical(i) && isscalar(i));
      p.addParameter('scalar_as_strings',   false, @(i) islogical(i) && isscalar(i));
      p.addParameter('directories_only',    false, @(i) islogical(i) && isscalar(i));
      p.addParameter('files_only',          false, @(i) islogical(i) && isscalar(i));
      p.addParameter('stop_if_empty',       false, @(i) islogical(i) && isscalar(i));
      % parse it
      p.parse(in,varargin{:});
      %wildcard character '*' needs to be translated to '.*' (if not already)
      in=translate_wildcard(in);
      %finding relevant directory
      [a,f,e]=fileparts(in);
      if isempty(a)
          a='.';
      end
      %rebuilding complete filename
      f_in=[f,e];
      %fetching directory listing
      f=dir(a);
      f=struct2cell(f);
      %branch on options
      if p.Results.directories_only
        idx=[f{4,:}];
      elseif p.Results.files_only
        idx=~[f{4,:}];
      else
        idx=true(1,size(f,2));
      end
      f=f(1,idx);
      %greping
      f1=reduce_cell_array(regexp([f(:)],['^',f_in,'$'],'match'));
      %branch on options
      if p.Results.prefer_non_mat_files
        filenames=unique(f1);
      else
        %taking into account also mat files
        if p.Results.strict_mat_search
          f2=reduce_cell_array(regexp([f(:)],['^',f_in,'\.mat$'],'match'));
        else
          f2=reduce_cell_array(regexp([f(:)],['^',f_in,'.*\.mat$'],'match'));
        end
        %preferentially load mat files and not original files
        for i=1:length(f2)
          f2_now=strrep(f2{i},'.mat','');
          for j=1:length(f1)
            if strcmp(f2_now,f1{j})
              %delete non-mat file
              f1{j}='';
            end
          end
        end
        filenames=unique([f1(:);f2(:)]);
      end
      %one filename is frequently empty
      filenames(cellfun('isempty',filenames))=[];
      %inform user if no files found
      if isempty(filenames)
        msg=['found no files named ',in,'.'];
        if p.Results.stop_if_empty
          error(msg)
        end
        if p.Results.disp
          disp([mfilename,': warning: ',msg])
        end
        out='';
        return
      end
      %add path
      root=fileparts(in);
      out=cell(1,length(filenames));
      for i=1:length(filenames)
        out{i}=fullfile(root,filenames{i});
      end
      %inform user
      if p.Results.disp
        disp([mfilename,': found ',num2str(numel(out)),' filenames in wildcarded string ''',in,'''.'])
      end
      %convert to string if requested
      if p.Results.scalar_as_strings
        assert(numel(out)==1,['If ''scalar_as_strings'' is true, then results must be scalar, not with length ',numel(out),'.'])
        out=out{1};
      end
    end
    function out=ensuredir(filename)
      d=fileparts(filename);
      if ~exist(d,'dir')
        out=mkdir(d);
      else
        out=true;
      end
    end
    function [out,s]=find(varargin)
      com=['find ',strjoin(str.clean(varargin,'regex'),' ')];
      disp(com)
      [s,r]=system(com);
      if s~=0
        disp(r)
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
    function f=trailing_filesep(f)
      if iscellstr(f)
        isdir=cellfun(@(i) exist(i,'dir')~=0,f);
        f(isdir)=cellfun(@(i) file.trailing_filesep(i),f(isdir),'UniformOutput',false);
      	return
      end
      if f(end)~=filesep
        f=[f,filesep];
      end
    end
    function out=exist(in)
      out=~isempty(file.wildcard(in,...
        'disp',             false,...
        'directories_only', false,...
        'files_only',       false...
      ));
    end
    function out=datenum(in)
      fileinfo=dir(in);
      if exist(in,'dir')
        assert(fileinfo(1).name=='.',['Excepting the first entry of ''fileinfo'' to be relative to ''.'', not to ''',fileinfo(1).name,'''.'])
        out=fileinfo(1).datenum;
      else
        assert(isscalar(fileinfo),['Expecting ''fileinfo'' to be scalar, since ''',in,''' doesn not appears to be a directory.'])
        out=fileinfo.datenum;
      end
    end
    function out=newest(in,varargin)
      newest_date=0;
      file_list=file.wildcard(in,varargin{:});
      for i=1:numel(file_list)
        if file.datenum(file_list{i})>newest_date
          out=file_list{i};
        end
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
function out = translate_wildcard(in)
  if ~isempty(strfind(in,'*'))
    idx=strfind(in,'*');
    for i=1:numel(idx)
      if idx(i) > 1 && ~strcmp(in(idx(i)-1),'.')
        out=translate_wildcard([in(1:idx(i)-1),'.',in(idx(i):end)]);
        return
      end
      if idx(i) == 1
        out=translate_wildcard(['.',in(idx(i):end)]);
        return
      end
    end
  end
  out=in;
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
