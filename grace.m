classdef grace
  properties(Constant)
    l1bdir_options={...
      '~/data/grace';...
      './data/grace';...
    };
  end
  methods(Static)
    function out=dir(type)
      switch type
        case 'l1b'
          for i=1:numel(grace.l1bdir_options)
            if file.exist(grace.l1bdir_options{i})
              out=grace.l1bdir_options{i};
              return
            end
          end
        %add more directories here
      end
    end
  end
end