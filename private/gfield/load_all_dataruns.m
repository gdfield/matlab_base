
function data=load_all_dataruns(base_path,data_path,ndf_names,wnruns,varargin)

% dataruns=load_all_dataruns(base_path,data_path,ndf_names,wnruns,marks_params);

%inputs: path, datarun names (ndfs), base path, data structure. data_path
%and names are cells. wnruns are indices within ndf_names that are white
%noise (and need stas loaded)
%outputs: dataruns within data structure

%get full data paths
full_paths=cell(1,length(data_path));
for i=1:length(data_path)
    full_paths(i)=strcat(base_path,data_path(i));    
end

k=1;
%load all info about RFs and build into a structure
for i=1:length(data_path)
    placeholder=load_data(full_paths{i});
    placeholder=load_neurons(placeholder);
    placeholder=load_params(placeholder);
    placeholder=load_ei(placeholder,'all'); 
    i
    if ~isempty(wnruns)
        if i==wnruns(k) 
            placeholder=load_sta(placeholder,'load_sta','all');
            if length(varargin) == 0
                placeholder=get_sta_summaries(placeholder,'all');
            elseif length(varargin{1}) == 1
                mp = varargin{1};
                placeholder=get_sta_summaries(placeholder,'all', 'marks_params', mp);
            end
            k=k+1;
        end
    end
    
    data.(ndf_names{i})=placeholder;    
end
clear placeholder

end

 