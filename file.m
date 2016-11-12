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
    function [fid,filename,close_file]=open(filename)
      %open the file (if a filename is given)
      if file.isfid(filename)
        %rename input
        fid=filename;
        %get name of the already open file
        filename=fopen(fid);
        %don't close this file afterwards
        close_file=false;
      else
        %'filename' is an actual filename, so open it
        fid = fopen(filename,'r');
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
        numeric_flag(i) = isnumstr(header{i});
        nr_columns(i)=numel(strsplit(header{i}));
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
        nlines = find(numeric_flag,1,'last');
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
    function [data,header] = textscan(filename,format,nlines)
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
        'emptyvalue',0 ...
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