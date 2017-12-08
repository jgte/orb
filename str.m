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
        s=strrep(s,'_','\_');
      case 'fieldname'
        if ~exist('alt_char','var') || isempty(alt_char)
          alt_char='';
        end
        s=str.rep(s,...
          '-',alt_char,...
          ' ',alt_char,...
          '.',alt_char,...
          '+',alt_char);
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
      s=dbstack(1);
      if isempty(s)
        disp(str.show(varargin(start_idx:end)))
      else
        si=min([numel(s),1+stack_delta]);
        disp([s(si).name,':',num2str(s(si).line),': ',str.show(varargin(start_idx:end))])
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
    %pads all entries of in with blanks until they all have the same length
    function out=cellpad(in)
      assert(iscellstr(in),['Need cell string, not ',class(in),'.'])
      n=max(cellfun(@numel,in));
      out=cellfun(@(i) [i,repmat(' ',1,n-numel(i))],in,'UniformOutput',false);
    end
  end
end