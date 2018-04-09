%This object defines details of the product
classdef nrtdm_metadata
  properties(GetAccess = 'private', SetAccess = 'private')
    debug
  end
  properties(Dependent)
    units
    fields
    fields_no_units
    dimension
  end
  properties(GetAccess = 'public', SetAccess = 'private')
    product
    entries
  end
  methods(Static)
    function out=clean_entry_name(entry_name)
      out=strrep(entry_name,...
        '-','_'...
      );
    end
    function test
      if isempty(dir(nrtdm_product.config_dir))
        disp([mfilename,':WARNING: cannot find NRTDM config dir: ',nrtdm_product.config_dir,'. Skipping test.'])
        return
      end
      a=nrtdm_metadata('SC_Basic/Quaternion_Interpolated');
      %TODO: finish this test
    end
  end
  methods
    function obj=nrtdm_metadata(product_name,debug)  
      if ~exist('debug','var') || isempty(debug)
        debug=false;
      end
      %initialize
      obj.product=nrtdm_product(product_name);
      obj.debug=debug;
      %load all metadata
      text = fileread(obj.product.file);
      %discover newlines
      newlines=[0,strfind(text,newline)];
      %sanitize
      if isempty(newlines)
        error([mfilename,': could not discover any new lines in file ',obj.product.file,':',newline,text])
      end
      %search for this product
      line_nr=0;
      for i=1:numel(newlines)-1
        if str.contains(text(newlines(i)+1:newlines(i+1)),['Orbital: ',obj.product.name])
          line_nr=i;
          break
        end
      end
      %sanity
      if line_nr==0
        error([mfilename,': could not definition of product ',obj.product.name,' in file ',obj.product.file,'.'])
      end
      %init metadata structure
      obj.entries=struct('Orbital',obj.product.name);
      %parse the entries of this product
      for i=line_nr+1:numel(newlines)-1
        %get this line
        line=text(newlines(i)+1:newlines(i+1));
        %get separator index within this line
        sep_index=strfind(line,':');
        if ~isempty(sep_index)
          %extract entry name
          entry_name=strtrim(line(1:sep_index-1));
          %extract entry value
          entry_value=strtrim(line(sep_index+1:end-1));
          %add to structure 
          if isfield(obj.entries,nrtdm_metadata.clean_entry_name(entry_name))
            %if already there, append
            if iscell(obj.entries.(nrtdm_metadata.clean_entry_name(entry_name)))
              obj.entries.(nrtdm_metadata.clean_entry_name(entry_name)){end+1}=entry_value;
            else
              obj.entries.(nrtdm_metadata.clean_entry_name(entry_name))={obj.entries.(nrtdm_metadata.clean_entry_name(entry_name)),entry_value};
            end
          else
            %otherwise, create new entry
            obj.entries.(nrtdm_metadata.clean_entry_name(entry_name))=entry_value;
          end
        else
          %we're done
          break
        end
      end
    end
    function out=get(obj,entry_name)
      if ~isfield(obj.entries,nrtdm_metadata.clean_entry_name(entry_name))
        error([mfilename,': metadata of product ',obj.product.name,' does not include entry ',entry_name,'.'])
      else
        out=obj.entries.(nrtdm_metadata.clean_entry_name(entry_name));
      end
    end
    function out=full_entries(obj)
      c=cell(1,obj.dimension);
      for i=1:numel(c)
        c{i}=[obj.product.str,'/',num2str(i),' '];
      end
      out=[c{:}];
    end
    function out=get.dimension(obj)
      out=str2double(obj.entries.('dimension'));
    end
    function out=get.units(obj)
      out=cell(1,obj.dimension);
      for i=1:numel(out)
        if isfield(obj.entries,['field_',num2str(i,'%02i'),'_unit'])
          out{i}=obj.entries.(['field_',num2str(i,'%02i'),'_unit']);
        else
          desc=obj.entries.(['field_',num2str(i,'%02i')]);
          idx=[strfind(desc,'('),strfind(desc,')')];
          if numel(idx) ~= 2
            out{i}='?';
          else
            out{i}=desc(idx(1)+1:idx(2)-1);
          end
        end
      end
      for i=1:numel(out)
        %translate '/s/s' to /s^2'
        idx=strfind(out{i},'/s/s');
        if ~isempty(idx)
          out{i}=strrep(out{i},'/s/s','/s^2');
        end
      end
    end
    function out=get.fields_no_units(obj)
      out=cell(1,obj.dimension);
      for i=1:numel(out)
        if isfield(obj.entries,['field_',num2str(i,'%02i'),'_descr'])
          out{i}=obj.entries.(['field_',num2str(i,'%02i'),'_descr']);
        else
          desc=obj.entries.(['field_',num2str(i,'%02i')]);
          idx=strfind(desc,'(');
          switch numel(idx)
          case 0
            out{i}=desc;
          otherwise
            out{i}=desc(1:min(idx)-1);
          end
        end
      end
    end
    function out=data_filename(obj,t_now,extension,data_dir)
      if ~exist('extension','var') || isempty(extension)
        extension='orbit';
      end
      if ~exist('data_dir','var') || isempty(data_dir)
        data_dir=fullfile(nrtdm.data_dir,extension);
      end      
      parent_dir=fullfile(data_dir,obj.product.str);
      if isempty(dir(parent_dir))
        mkdir(parent_dir)
      end
      out=fullfile(parent_dir,[datestr(t_now,'yyyy'), '_',num2str(day(t_now,'dayofyear'),'%03i'),'.',extension]);
    end
  end
end











