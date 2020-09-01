classdef url
  properties(Constant)
  type_list={...
'http','https','ftp','scp'
  };
  type_set='://';
  end
  methods(Static)
    function test
      u='http://domain.something.edu/dir1/page.html';
      disp(' - u -')
      disp(u)
      disp(' - url.type(u) -')
      disp(url.type(u))
      disp(' - url.isurl(u) -')
      disp(url.is(u))
      disp(' - url.address(u) -')
      disp(url.address(u))
    end
    function out=type(in)
      out=strsplit(in,url.type_set);
      out=out{1};
    end
    function out=is(in)
      out=any(contains(url.type_list,url.type(in)));
    end
    function out=address(in)
      out=strsplit(in,url.type_set);
      out=out{2};
    end
  end
end