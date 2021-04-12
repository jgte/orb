classdef machinery
  methods(Static)
    function [out,success]=trycatch(success,errorID,func,args)
      if ~sucess
        try 
          out=func(args{:});
          success=true;
        catch ME
          switch ME.identifier
            case errorID
              %do nothing
            otherwise
              rethrow(ME)
          end
        end
      end
    end
  end
end