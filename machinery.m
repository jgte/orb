classdef machinery
  methods(Static)
    %% handy overloads
    function out=isa(var,types)
      if iscellstr(types)
        out=any(cellfun(@(i) isa(var,i),types));
      else
        out=isa(var,types);
      end
    end
    %% flow control
    function [out,success]=trycatch(success,errorID,func,args)
      if ~success
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
      else
        out=[];
      end
    end
  end
end