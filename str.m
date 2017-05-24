classdef str
  methods(Static)
    function test
      disp(' - str.rand -')
      for i={'u','l','n'}
        disp(str.rand(10,1,i{1}))
      end
      disp(' - str.show -')
      for i={1,'1',true,datetime('now'),seconds(1),{2,'2',false,datetime('now'),seconds(2)}}
        disp(str.show(i{1}))
      end
      disp(' - str.tabbed -')
      for j={true,false}
        disp(['right_justified:',str.show(j{1})])
        for i=1:10
          disp(['<',str.tabbed('str',i,j{1}),'>'])
        end
      end
      for i=4:4:16
        disp([' - str.tablify: tab=',num2str(i),' -'])
        disp(str.tablify(i,'1','2','3','4'))
        disp(str.tablify(i,[1,2,3,randn(1)]))
        disp(str.tablify(i,randn(1),'2',{'3',4}))
      end
    end
    function out=rand(n,l,mode)
      if ~exist('l','var') || isempty(l)
        l=1;
      end
      if ~exist('mode','var') || isempty(mode)
        mode='l';
      end
      switch lower(mode)
        case {'u','upper','caps'}
          ascii_start=65;
          ascii_stop=90;
        case {'l','lower','noncaps'}
          ascii_start=97;
          ascii_stop=122;
        case {'n','numeric'}
          ascii_start=48;
          ascii_stop=57;
        otherwise
          error([mfilename,': unknown mode ''',mode,'''.'])
      end
      out=char(floor((ascii_stop-ascii_start)*rand(l,n)) + ascii_start);
    end
%     function out=show(in)
%       if isnumeric(in)
%         out=num2str(in);
%       elseif ischar(in)
%         out=in;
%       elseif islogical(in)
%         if in
%           out='T';
%         else
%           out='F';
%         end
%       elseif isdatetime(in)
%         out=datestr(in,'yyyy-mm-dd HH:MM:SS.FFF');
%       elseif isduration(in)
%         out=char(in);
%       elseif iscell(in)
%         out=strjoin(cellfun(@(i)([str.show(i),'; ']),in,'UniformOutput',false));
%       else
%         error([mfilename,': cannot handle variables of class ',class(in),'.'])
%       end
%     end
    function out=show(in,fmt)
      %trivial call
      if ischar(in)
        out=in;
        return
      end
      %handle non-scalar quantities
      if ~isscalar(in)
        out=cell(size(in));
        for i=1:numel(in)
          out{i}=str.show(in(i));
        end
        out=strjoin(out,' ');
        return
      end
      %branch on s
      %branch on scalar type
      switch class(in)
      case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64','single','double'}
        if exist('fmt','var')
          out=num2str(in,fmt);
        else
          out=num2str(in);
        end
      case 'logical'
        if in, out='T';
        else   out='F';
        end
      case 'cell'
        out=str.show(in{1}); %non-scalar cells already handled above
      case 'datetime'
        if isnat(in)
          out='NaT';
        elseif ~isfinite(in)
          out='Inf';
        else
          out=datestr(in);
        end
      case 'duration'
        out=time.str(seconds(in));
      case 'struct'
        out=[10,'  ',strjoin(strsplit(structs.str(in),char(10)),[10,'  ']),10];
      otherwise
        try
          out=in.str;
        catch
          error([mfilename,': class ''',class(in),''' is not supported.'])
        end
      end
    end
    function out=tabbed(in,tab,right_justified)
      if ~exist('right_justified','var') || isempty(right_justified)
        right_justified=false;
      end
      if right_justified
        out=str.just(in,tab,'just','right');
      else
        out=str.just(in,tab,'just','left');
      end
    end
    function s=clean(s,mode)
      if iscellstr(mode)
        for i=1:numel(mode)
        	s=str.clean(s,mode{i});
        end
        return
      end
      %trivial call
      if isempty(mode)
        return
      end
      %branch on mode
      switch lower(mode)
      case 'basename'
        [~,s,e]=fileparts(s);
        s=[s,e];
      case 'file'
        [~,s]=fileparts(s);
      case 'succ_blanks'
        while ~isempty(strfind(s,'  '))
          s=strrep(s,'  ',' ');
        end
        s=strtrim(s);
      case '_'
        s=str.clean(strrep(s,'_',' '),'succ_blanks');
      case '.'
        s=str.clean(strrep(s,'.',' '),'succ_blanks');
      case 'title'
        s=strrep(s,'_','\_');
      case 'grace'
        sats={'A','B'};
        names={'gr<SAT>.','gr<SAT>.','G<SAT>_','G<SAT>_'};
        for i=1:numel(sats)
          for j=1:numel(names)
            part=strrep(names{j},'<SAT>',sats{i});
            s=strrep(s,part,['GRACE-',sats{i},' ']);
          end
        end
      otherwise
        s=strrep(s,mode,'');
      end
    end
    function out=ispresent(parser,field)
      %sanity
      if iscell(field)
        out=cellfun(@str.ispresent,parser,field);
        return
      end
      % look for existence
      if any(strcmp(parser.Parameters,field))
        out=~any(strcmp(parser.UsingDefaults,field));
      else
        out=isfield(parser.Unmatched,field);
      end
    end
    function out=just(in,len,varargin)
      % parse mandatory arguments
      p=inputParser;
      p.addRequired( 'in',           @(i) ischar(i));
      p.addRequired( 'len',          @(i) isfinite(i));
      p.addParameter('just','center',@(i) any(strcmp(i,{'left','center','right'})));
      p.parse(in,len,varargin{:});
      %truncate 'in' if needed
      if numel(in)>len; in=in(1:len); end
      %call mother routine
      out=strjust([repmat(' ',1,len-numel(in)),in],p.Results.just);
    end
    function out=th(i)
      switch i
      case 1
        out='1st';
      case 2
        out='2nd';
      case 3
        out='3rd';
      otherwise
        out=[num2str(i),'th'];
      end
    end
    function sizetrap(var1,var2)
      if numel(var1)~=numel(var2)
        throwAsCaller(MException([mfilename,':SizeTrap'],[...
          'size of variable ''',inputname(1),''' (',num2str(numel(var1)),') different than ',...
          'size of variable ''',inputname(2),''' (',num2str(numel(var2)),'). This is ilegal.'...
        ]))
      end
    end
    %first argument is field width, all remaining inputs are values to print.
    function out=tablify(w,varargin)
      %justification scheme
      mode='center';
      %expand scalar widths
      if isscalar(w)
        w=ones(size(varargin))*w;
      else
        assert(numel(varargin)==numel(w),[mfilename,': ',...
          'number of elements of vector input ''w'' (',numel(w),') must be the same as the ',...
          'number of additional input arguments (',numel(varargin),').'])
      end
      out = cell(size(varargin));
      c=0;
      for i=1:numel(varargin)
        for j=1:numel(varargin{i})
          c=c+1;
          out{c} = str.just(str.show(varargin{i}(j),['%',num2str(w(i)),'g']),w(i),'just',mode);
        end
%         switch class(varargin{i})
%         case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64','single','double','logical'}
%           if isscalar(varargin{i})
%             out{i} = str.just(num2str(varargin{i},['%',num2str(w(i)),'g']),w(i),'just',mode);
%           else
%             out{i} = str.tablify(w(i),num2cell(varargin{i}));
%           end
%         case 'char'
%           out{i} = str.just(varargin{i},w(i),'just',mode);
%         case 'cell'
%           out_now=cell(size(varargin{i}));
%           for j=1:length(varargin{i})
%             out_now{j} = str.tablify(w(i),varargin{i}{j});
%           end
%           out{i}=out_now;
%         case 'datetime'
%           out{i}=str.just(datestr(varargin{i}),w(i),'just',mode);
%         case 'duration'
%           out{i}=str.just(time.str(seconds(varargin{i})),w(i),'just',mode);
%         otherwise
%           error([mfilename,': class ''',class(varargin{i}),''' of the ',str.th(i),' input argument is not supported.'])
%         end
      end
      out=strjoin(out,' ');
    end
    function out=num(in)
      out1=regexprep(in,'(\d)[Dd]([-+\d])','$1e$2');  %replace D and d exponents in scientific notation with e
      out2=strsplit(out1);                            %split string along blank characters
      out3=out2(cellfun(@(i) ~isempty(i),out2));      %remove empty cells
      out=cellfun(@(i) sscanf(i,'%f'),out3);          %convert to numeric cells
      %sanity
      if ~isnumeric(out)
        error([mfilename,': convertion from string to vector failed.'])
      end
    end
    function out=titlecase(in,separator)
      if ~exist('separator','var') || isempty(separator)
        separator=' ';
      end
      out=strsplit(in,separator);
      for i=1:numel(out)
        if ~isempty(out{i})
          out{i}=[upper(out{i}(1)),lower(out{i}(2:end))];
        end
      end
      out=strjoin(out,separator);
    end
  end
end