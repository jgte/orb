%This object defines details of the product
classdef nrtdm_product
  properties(GetAccess = 'private', SetAccess = 'private')
    sep_index
    debug
  end
  properties(Dependent)
    category
    name
    sat
    field
    file
  end
  properties(GetAccess = 'public', SetAccess = 'private')
    str
  end
  methods(Static)
    function test
      if isempty(dir(nrtdm.config_dir))
        disp([mfilename,':WARNING: cannot find NRTDM config dir: ',nrtdm.config_dir,'. Skipping test.'])
        return
      end
      %TODO: complete this test
    end
  end
  methods
    function obj=nrtdm_product(product_name,debug)  
      if ~exist('debug','var') || isempty(debug)
        debug=false;
      end
      if (debug); disp(['start:',current_method]); end
      %initialize
      obj.str=product_name;
      obj.debug=debug;
      obj.sep_index=strfind(product_name,'/');
      %sanity
      if (numel(obj.sep_index) < 1) || (numel(obj.sep_index) > 2)
        error([mfilename,'Can not handle product ''',obj.str,'''.'])
      end
      if numel(obj.sep_index)==1
        obj.sep_index=[obj.sep_index,length(product_name)+1];
      end
      if (obj.debug); disp(['end  :',current_method]); end
    end
    function out=get.category(obj)
      if (obj.debug); disp(['start:',current_method]); end
      out=obj.str(1:obj.sep_index(1)-1);
      if (obj.debug); disp(['end  :',current_method,':',out]); end
    end
    function out=get.name(obj)
      if (obj.debug); disp(['start:',current_method]); end
      out=obj.str(obj.sep_index(1)+1:obj.sep_index(2)-1);
      if (obj.debug); disp(['end  :',current_method,':',out]); end
    end
    function out=get.sat(obj)
      if (obj.debug); disp(['start:',current_method]); end
      out=obj.str(1:min(strfind(obj.str,'_')-1));
      if (obj.debug); disp(['end  :',current_method,':',out]); end
    end
    function out=get.field(obj)
      if (obj.debug); disp(['start:',current_method]); end
      if obj.sep_index(2)>length(obj.str)
        out=0;
      else
        out=obj.str(obj.sep_index(2)+1:end);
        field_index=strfind(out,'-');
        if isempty(field_index)
          out=str2double(out);
        else
          out=str2double(out(1:field_index-1)):str2double(out(field_index+1:end));
        end
      end
      if (obj.debug); disp(['end  :',current_method,':',num2str(out)]); end
    end
    function out=get.file(obj)
      if (obj.debug); disp(['start:',current_method]); end
      out=fullfile(nrtdm.config_dir,[obj.category,'.products.txt']);
      if isempty(dir(out))
        error([mfilename,': Cannot find metadata files ''',out,'''.'])
      end
      if (obj.debug); disp(['end  :',current_method,':',out]); end
    end
  end
end

function out=current_method
  s=dbstack;
  out=s(2).name;
end