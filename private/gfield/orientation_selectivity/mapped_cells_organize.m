function [dataruns,types_ids,types_inx,mapped_ids_by_types_interest]=mapped_cells_organize(dataruns,ndf_names,types_interest)
%see map_cells for info

%how many cells map across dataruns from using map analysis
mapped_ids=dataruns.(ndf_names{end}).cell_ids;
for i=1:length(ndf_names)-1
    mapped_ids=intersect(mapped_ids,dataruns.(ndf_names{i}).cell_ids); 
end

%record indices in each datarun by all mapped cells
for j=1:length(ndf_names)
    temp_datarun=dataruns.(ndf_names{j});
    dataruns.(ndf_names{j}).mapped_inx=get_cell_indices(temp_datarun,mapped_ids);
end

%how many of those are typed from the master datarun

%make sure type names are acceptable
for j=1:length(dataruns.(ndf_names{end}).cell_types)
    if ~isempty(dataruns.(ndf_names{end}).cell_types{j}.name)
        dataruns.(ndf_names{end}).cell_types{j}.name=genvarname(dataruns.(ndf_names{end}).cell_types{j}.name);
    end
end
for i=1:length(dataruns.(ndf_names{end}).cell_types)
    mapped_ids_by_type{i}=intersect(mapped_ids,dataruns.(ndf_names{end}).cell_types{i}.cell_ids);
end
%get names and indices of the types from datarun master
k=1;
for i=1:length(dataruns.(ndf_names{end}).cell_types)
    if ~isempty(dataruns.(ndf_names{end}).cell_types{i}.name)
        types_name{k}=dataruns.(ndf_names{end}).cell_types{i}.name;
        types_minx(k)=i;
        k=k+1;
    end
end
for i=1:length(types_interest)
    for j=1:length(types_name)
        if isequal(types_interest{i},types_name{j})
            types_interest_minx(i)=types_minx(j);
            continue
        end
    end
end
%ids that have mapped for each type
mapped_ids_by_types_interest=mapped_ids_by_type(types_interest_minx);
for i=1:length(mapped_ids_by_types_interest)
    temp(i)=length(mapped_ids_by_types_interest{i});
end
disp(strcat('total number of typed, mapped cells that I care about is ', num2str(sum(temp)))) 

%indices for those cells per datarun
for i=1:length(types_interest)
    types_ids(i)={strcat(types_interest{i},'_ids')}; %not really needed
    types_inx(i)={strcat(types_interest{i},'_inx')};
end

%record indices in each datarun by types_interest
for i=1:length(types_interest)
    for j=1:length(ndf_names)
        temp_datarun=dataruns.(ndf_names{j});
        dataruns.(ndf_names{j}).(types_inx{i})=get_cell_indices(temp_datarun,mapped_ids_by_types_interest{i});
    end
end



end
