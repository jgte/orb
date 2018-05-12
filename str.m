classdef str
  properties(Constant)
  latex_keywords={...
'\alpha','\upsilon','\sim','\angle','\phi','\leq','\ast','\chi','\infty','\beta','\psi','\clubsuit','\gamma',...
'\omega','\diamondsuit','\delta','\Gamma','\heartsuit','\epsilon','\Delta','\spadesuit','\zeta','\Theta','\leftrightarrow',...
'\eta','\Lambda','\leftarrow','\theta','\Xi','\Leftarrow','\vartheta','\Pi','\uparrow','\iota','\Sigma','\rightarrow',...
'\kappa','\Upsilon','\Rightarrow','\lambda','\Phi','\downarrow','\mu','\Psi','\circ','\nu','\Omega','\pm','\xi','\forall',...
'\geq','\pi','\exists','\propto','\rho','\ni','\partial','\sigma','\cong','\bullet','\varsigma','\approx','\div','\tau',...
'\Re','\neq','\equiv','\oplus','\aleph','\Im','\cup','\wp','\otimes','\subseteq','\oslash','\cap','\in','\supseteq','\supset',...
'\lceil','\subset','\int','\cdot','\o','\rfloor','\neg','\nabla','\lfloor','\times','\ldots','\perp','\surd','\prime','\wedge',...
'\varpi','\0','\rceil','\rangle','\mid','\vee','\langle','\copyright'....
  };
  latex_chars={'{','}'};
  end
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
    function out=show(in,fmt,join_char)
      %trivial call
      if ischar(in)
        out=in;
        return
      end
      %handle non-scalar quantities
      if ~isscalar(in)
        out=cell(size(in));
        for i=1:numel(in)
          if exist('fmt','var') && ~isempty(fmt)
            out{i}=str.show(in(i),fmt);
          else
            out{i}=str.show(in(i));
          end
        end
        if ~exist('join_char','var')
          join_char=' ';
        end
        out=strjoin(out,join_char);
        return
      end
      %branch on s
      %branch on scalar type
      switch class(in)
      case {'int8','uint8','int16','uint16','int32','uint32','int64','uint64','single','double'}
        if exist('fmt','var') && ~isempty(fmt)
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
        if ~exist('join_char','var')
          out=['  ',strjoin(strsplit(structs.str(in,'','',false),char(10)),[10,'  '])];
        else
          if strcmp(join_char,'_')
            out=structs.str(in,'','_',false);
          else
            out=strjoin(strsplit(structs.str(in,'','',false),char(10)),join_char);
          end
        end
      case 'function_handle'
        out=func2str(in);
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
    function out=columns(in,just,splitchar)
      if ~exist('just',     'var') || isempty(just     );      just='right'; end
      if ~exist('splitchar','var') || isempty(splitchar); splitchar=' ';     end
      assert(iscellstr(in),     ['Input ''in'' must be a cell string, not a ',class(in),'.'])
      assert(numel(size(in))==2,['Input ''in'' must be a 1D or 2D cell string array, not ',numel(size(in)),'D.'])
      if any(size(in)==1)
        %entries in s will become rows
        rows=numel(in);
        %find out the number of columns
        cols=cellfun(@(i) numel(strsplit(i,splitchar,'CollapseDelimiters',true)),in);
        %make room for outputs
        out=cell(rows,max(cols));out(:)={''};
        %loop over all rows
        for i=1:rows
          %get columns for this row
          row_now=strsplit(strtrim(in{i}),splitchar,'CollapseDelimiters',true);
          %distribute data
          switch lower(just)
          case {'left','center'}
            out(i,end-numel(row_now)+1:end)=row_now;
          case {'right'}
            out(i,1:numel(row_now))=row_now;
          otherwise
            error(['Cannot handle ''just'' with value ''',just,'''.'])
          end
        end
        %propagate
        in=out;
      end
      %loop over columns
      for i=1:size(in,2)
        %get maximum width of this column
        tab=max(cellfun(@(i) length(strtrim(i)),in(:,i)));
        %justify
        in(:,i)=cellfun(@(i) str.just(i,tab,'just',just),in(:,i),'UniformOutput',false);
      end
      %make room for outputs
      out=cell(rows,1);
      %glue columns together, loop over lines
      for i=1:size(in,1)
        out{i}=strjoin(in(i,1:cols(i)),splitchar);
      end
      %re-align everything, get maximum width of trimmed column, account for latex string
      tab_latex=max(cellfun(@(i) str.latex_length(strtrim(i)),out));
      tab      =max(cellfun(@(i)           length(strtrim(i)),out));
      %justify, only if not latex; if latex justify with care
      for i=1:numel(out)
        if str.islatex(out{i})
          out{i}=str.just(strtrim(out{i}),tab,'just',just);
        else
          out{i}=str.just(str.latex_trim(out{i}),tab_latex,'just',just);
        end
      end
    end
    %% latex stuff
    function out=islatex(in)
      out=~isempty(strfind(in,'\'));
    end
    %account for latex keywords and chars when calculating the length of a string
    function out=latex_length(in)
      %usual way to compute the length
      out=length(in);
      %reduce the search, by checking if there's any '\'
      if ~str.islatex(in); return; end
      %loop over all latex keyword
      for i=1:numel(str.latex_keywords)
        if ~isempty(strfind(in,str.latex_keywords{i}))
          out=out-length(str.latex_keywords{i})+1;
        end
      end
      %handle special chars too
      for i=1:numel(str.latex_chars)
        out=out-numel(strfind(in,str.latex_chars{i}));
      end
    end
    %only trim if it's not a latex string
    function out=latex_trim(in)
      if str.islatex(in)
        out=strtrim(in);
      else
        out=in;
      end
    end
    %transforms a matrix into the data part of a latex table
    function out=latex_table(in,fmt)
      if ~exist('fmt','var'); fmt='%g'; end
      assert(isnumeric(in) && numel(size(in)<=2),'Need numeric matrix/vector/scales input')
      out=cell(size(in));
      for i=1:size(in,1)
        for j=1:size(in,2)
          if j==size(in,2)
            out{i,j}=[num2str(in(i,j),fmt),' \\',newline];
          else
            out{i,j}=[num2str(in(i,j),fmt),' & '];
          end
        end
      end
      out=strjoin(transpose(out));
    end
    %% string manipulation
    function s=clean(s,mode,alt_char)
      if iscellstr(mode)
        for i=1:numel(mode)
        	s=str.clean(s,mode{i});
        end
        return
      end
      if iscell(s)
        if isempty(cells.rm_empty(s))
          s='';
          return
        end
        assert(iscellstr(s),'Can only handle cells of strings.')
      end
      %trivial calls
      if isempty(mode);return;end
      if isempty(s   );return;end
      
      %branch on mode
      switch lower(mode)
      case 'basename'
        [~,s,e]=fileparts(s);
        s=[s,e];
      case 'file'
        [~,s]=fileparts(s);
      case 'succ_blanks'
        s=str.clean(s,'succ',' ');
      case 'succ'
        assert(exist('alt_char','var')~=0,['If mode is ''',mode,''', then need ''alt_char'' input arg.'])
        while ~isempty(strfind(s,[alt_char,alt_char]))
          s=strrep(s,[alt_char,alt_char],alt_char);
        end
        if alt_char==' '
          s=strtrim(s);
        end
      case 'title'
        s=strrep(s,'_',' ');
      case 'fieldname'
        if ~exist('alt_char','var') || isempty(alt_char)
          alt_char='';
        end
        s=str.rep(s,...
          '/',alt_char,...
          '-',alt_char,...
          ' ',alt_char,...
          '.',alt_char,...
          '+',alt_char);
      case {'regex','regexp'}
        s=str.rep(s,...
          '[','\[',...
          ']','\]',...
          '+','\+',...
          '?','\?',...
          '*','\*');
      case {'noregex','noregexp'}
        s=str.clean(s,{'[',']','?','+','*'});
      case '<grace>'
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
    function s=rep(s,varargin)
      assert(mod(numel(varargin),2)==0,'Need pairs of input arguments following ''s''.')
      for i=1:numel(varargin)/2
        s_out=strrep(s,varargin{2*i-1},varargin{2*i});
        s=s_out;
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
      p.addRequired( 'in',           @(i) ischar(i) || is);
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
    %pads all entries of in with blanks until they all have the same length
    function out=cellpad(in)
      assert(iscellstr(in),['Need cell string, not ',class(in),'.'])
      n=max(cellfun(@numel,in));
      out=cellfun(@(i) [i,repmat(' ',1,n-numel(i))],in,'UniformOutput',false);
    end
    function out=common(in,splitchar)
      assert(iscellstr(in),['Input ''in'' has to be of class cellstr, not ',class(in),'.'])
      in=cellfun(@(i) strsplit(i,splitchar),in,'UniformOutput',false);
      out=in{1};
      for i=2:numel(in)
        out=intersect(out,in{i},'stable');
      end
      out=strjoin(out,splitchar);
    end
    function out=unique(in,splitchar)
      if ~exist('splitchar','var')
        out=unique(in);
        return
      end
      common=strsplit(str.common(in,splitchar),splitchar);
      in=cellfun(@(i) strsplit(i,splitchar),in,'UniformOutput',false);
      out=cell(size(in));
      for i=1:numel(in)
        o=setdiff(in{i},common,'stable');
        if isempty(o)
          out{i}='';
        elseif iscell(o)
          out{i}=strjoin(o,splitchar);
        else
          out{i}=o;
        end
      end
    end
    function out=contains(str,pattern)
      out=~isempty(regexp(str,['.*',pattern,'.*'],'once'));
    end
    function io=trunc(io,n,default)
      if isempty(io)
        io=default;
      end
      if length(io)>n
        io=io(1:min([n,length(io)]));
      end
    end
    %% conversion
    %first argument is field width, all remaining inputs are values to print.
    function out=tablify(w,varargin)
      %justification scheme
      mode='center';
      %propagate all arguments (may need to edit them)
      in=varargin;
      %expand scalar widths
      if isscalar(w)
        w=ones(size(in))*w;
      else
        assert(numel(in)==numel(w),[mfilename,': ',...
          'number of elements of vector input ''w'' (',num2str(numel(w)),') must be the same as the ',...
          'number of additional input arguments (',num2str(numel(in)),').'])
      end
      %make room for outputs (only an estimate, maybe there are non-scalar arguments)
      out = cell(size(in)); c=0;
      %loop over every argument
      for i=1:numel(in)
        %need to encapsulate strings in a cell array to avoid looping over all individual characters
        if ischar(in{i})
          in{i}=in(i);
        end
        %loop over every element of the current argument (even if scalar)
        for j=1:numel(in{i})
          c=c+1; out{c} = str.just(str.show(in{i}(j),['%',num2str(w(i)),'g']),w(i),'just',mode);
        end
      end
      out=strjoin(out,' ');
    end
    function out=num(in)
      %trivial call
      if isnumeric(in); out=in; return; end
      %handle keywords
      switch lower(in)
      case  'inf'; out= inf;return
      case '-inf'; out=-inf;return
      case  'nan'; out= NaN;return
      end
      out1=regexprep(in,'(\d)[Dd]([-+\d])','$1e$2');  %replace D and d exponents in scientific notation with e
      out2=strsplit(out1);                            %split string along blank characters
      out3=out2(cellfun(@(i) ~isempty(i),out2));      %remove empty cells
      out=cellfun(@(i) sscanf(i,'%f'),out3);          %convert to numeric cells
      %sanity
      assert(isnumeric(out),'Convertion from string to vector failed.')
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
    function out=logical(in,mode)
      if ~exist('mode','var') || isempty(mode)
        mode='truefalse';
      end
      %first make sure it's logical
      switch class(in)
      case 'logical'
        %do nothing
      case 'cell'
        out=cellfun(@str.logical,in);
      case 'datetime'
        out=(in~=datetime(0,0,0));
      case 'duration'
        out=(in~=seconds(0));
      otherwise
        try
          out=(in~=0);
        catch
          error([mfilename,': class ''',class(in),''' is not supported.'])
        end
      end
      %then convert to requested mode
      switch lower(mode)
      case 'truefalse'; if in; out='true'; else out='false';end
      case 'tf';        if in; out='T';    else out='F';    end
      case 'onoff';     if in; out='on';   else out='off';  end
      case 'yesno';     if in; out='yes';  else out='no';   end
      otherwise
        idx=strfind(mode,'-');
        assert(~isempty(idx),['Unknown mode ''',mode,'''.'])
        mode=strsplit(mode,'-');
        if in; out=mode{1}; else out=mode{2}; end
      end
    end
    %% user feedback
    function say(varargin)
      %default value for internal parameters
      stack_delta=0;
      %loop control
      start_idx=1;
      start_idx_old=0;
      %pass parameters to this method as the first two arguments
      while start_idx~=start_idx_old
        %update loop controls
        start_idx_old=start_idx;
        %check if this is an internal parameter
        switch varargin{start_idx}
        case 'stack_delta'
          stack_delta=varargin{start_idx+1};
          start_idx=start_idx+2;
        end
      end
      s=dbstack(1+stack_delta);
      if isempty(s)
        disp(str.show(varargin(start_idx:end)))
      else
        disp([s(1).name,':',num2str(s(1).line),': ',str.show(varargin(start_idx:end))])
      end
    end
    function log(filename,msg,varargin)
      % parse mandatory arguments
      p=inputParser;
      p.addRequired( 'filename',    @(i) ischar(i));
      p.addRequired( 'msg',         @(i) iscellstr(i) || ischar(i));
      p.addParameter('clear',false, @(i) islogical(i));
      p.parse(filename,msg,varargin{:});
      if p.Results.clear && exist(filename,'file')
        delete(filename)
      end
      if isempty(msg)
        return
      end
      if ischar(msg)
        tmp=cell(size(msg,1),1);
        for i=1:size(msg,1)
          tmp{i}=msg(i,:);
        end
        msg=tmp;
      end
      fid = fopen(filename,'a');  
      fprintf(fid,[strjoin(msg,'\n'),'\n']);
      fclose(fid);
    end
    %% under development
%     function out=fmt_f2c(in)
%       %split input fortran format string into a cell
%       in=split(in,',');
%       out=cell(size(in));
%       %loop over all entries
%       for i=1:numel(in)
%         count=0;
%         %loop over all characters
%         for j=1:length(in{i})
%           switch in{i}(j)
%           %check if it's a digit
%           case arrayfun(@(i) {num2str(i)},1:9)
%             count=count*10^(j-1)+str2double(in{i}(j));
%           case 'A'; out{i}='s';
%           case ''; out{i}='s';
%           end
%         end
%       end
%     end
  end
end