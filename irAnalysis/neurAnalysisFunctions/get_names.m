function [all_names, num_names]=get_names(name_templates, workspace_variables)
%creates cell array of variables starting with name_templates.
%name_templates is a cell array of names e.g., {'sig', 's1','s2'}
%workspace variables is output of 'who' command in base workspace.

num_templates=length(name_templates);
num_possible=length(workspace_variables);
num_names=0;

for i= 1 : num_possible
    tmp_nm=workspace_variables{i};
    for j=1:num_templates
        tmp_template=name_templates{j};
        if strfind(tmp_nm, tmp_template)
            num_names=num_names+1;
            all_names{num_names}=workspace_variables{i};
            break; %breaks out of j loop into i loop once match is found
        end %if
    end %for j
end  %for i

%