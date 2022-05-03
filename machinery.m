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
    %% parser with some defaults and quicker declaration
    function p=inputParser(varargin)
      %init the object      
      p=inputParser;
      %implement reasonable and safe input parser options
      p.KeepUnmatched=true;
      p.PartialMatching=false;
      %implement requested options (may overwrite the defaults forced above)
      for i=1:numel(varargin)/2
        p.(varargin{2*i-1})=varargin{2*i};
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