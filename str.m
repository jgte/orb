classdef str
  methods(Static)
    function test
      disp(' - str.rand -')
      for i={'u','l','n'}
        disp(str.rand(10,1,i{1}))
      end
      disp(' - str.show -')
      for i={1,'1',true,datetime('now'),seconds(1)}
        disp(str.show(i{1}))
      end
      disp(' - str.tabbed -')
      for j={true,false}
        disp(['right_justified:',str.show(j{1})])
        for i=1:10
          disp(['<',str.tabbed('str',i,j{1}),'>'])
        end
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
    function out=show(in)
      if isnumeric(in)
        out=num2str(in);
      elseif ischar(in)
        out=in;
      elseif islogical(in)
        if in
          out='true';
        else
          out='false';
        end
      elseif isdatetime(in)
        out=datestr(in,'yyyy-mm-dd HH:MM:SS.FFF');
      elseif isduration(in)
        out=char(in);
      elseif iscell(in)
        out=strjoin(cellfun(@(i)([str.show(i),'; ']),in,'UniformOutput',false));
      else
        error([mfilename,': cannot handle variables of class ',class(in),'.'])
      end
    end
    function out=tabbed(in,tab,right_justified)
      if ~exist('right_justified','var') || isempty(right_justified)
        right_justified=false;
      end
      if right_justified
        out=[repmat(' ',1,tab-numel(in)),in];
      else
        out=[in,repmat(' ',1,tab-numel(in))];
      end
    end
    function s=clean(s,mode)
      if iscellstr(mode)
        for i=1:numel(mode)
        	s=str.clean(s,mode{i});
        end
        return
      end
      switch lower(mode)
      case 'succ_blanks'
        while ~isempty(strfind(s,'  '))
          s=strrep(s,'  ',' ');
        end
        s=strtrim(s);
      case '_'
        s=str.clean(strrep(s,'_',' '),'succ_blanks');
      case 'title'
        s=strrep(s,'_','\_');
      otherwise
        error([mfilename,': unknown mode ''',mode,'''.'])
      end
    end
  end
end
    
