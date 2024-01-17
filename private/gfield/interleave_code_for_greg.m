% load stuff
base_path='/Volumes/lab/Experiments/Array/Analysis/2018-11-30-0';
types_interest={'off brisk sustained';'off brisk transient';'on brisk sustained';'on brisk transient'};
types_interest=genvarname(types_interest);
data_path={'/data001_KR-map/data001_KR-map';'/data002_KR-map/data002_KR-map';'/data003_KR/data003_KR'}; %master last
ndf_names={'NDF4';'NDF0';'NDF0m'};
wnruns=[3];%indices of wnruns
dataruns=load_all_dataruns(base_path,data_path,ndf_names,wnruns);
wn_runs_names={'NDF0m'};clc
[dataruns,types_ids,types_inx,mapped_ids_by_types_interest]=mapped_cells_organize(dataruns,ndf_names,types_interest);

for type = 1:length(types_interest)
    figure;plot_rf_summaries(dataruns.NDF0m, mapped_ids_by_types_interest{type})
end


% organize interleaves
%from stimulus file:
total_rep_time = 1125; %in seconds
total_nonrep_time = 4050;
%repeat lengths set to 125s
num_interleaves = total_rep_time/125; %should evenly divide
nonrep_length = total_nonrep_time/num_interleaves;
rep_length = 5; %in seconds, one repeat
num_reps_win = 125/rep_length; %evenly divide

%2430 total from nonrepeated (4050*60/100)
%675 total from repeated (1125*60/100)
%=3105 cool

for j = 1:2
    triggers(j,:) = dataruns.(ndf_names{j}).triggers;
    temp = triggers(j,:);
    nonreptrigs{j} = zeros(270,num_interleaves); %450*60/100 per interleave     1:(75+270):2070
    reptrigs{j} = zeros(75, num_interleaves); %125*60/100 per interleave
    for ni = 1:num_interleaves
        nonreptrigs{j}(:,ni) = temp( 345*(ni-1)+1 : 345*(ni-1)+1+269 );
        reptrigs{j}(:,ni) = temp( 345*(ni-1)+1+270 : 345*(ni-1)+1+270+74 );
    end     
end
%stim
seeds = [11111 22222 33333 44444 55555 66666 77777 88888 99999];
nm1s{1} = '/Volumes/lab/acquisition/movie-xml/BW-60-2-0.48-';
nm1s{2} = '/Volumes/lab/acquisition/movie-xml/BW-60-1-0.48-';
nm3s = '-13x10-60.35_xoffset.xml';
for j = 1:2
    tr = str2num(nm1s{j}(end-6));
    for ni = 1:num_interleaves
        nm2 = seeds(ni);
        [temp_mov,~,~,dur,refresh(j)] = get_movie([nm1s{j} num2str(nm2) nm3s],nonreptrigs{j}(:,ni),nonrep_length*60/tr); %num of unique frames
        mov{ni} = squeeze(temp_mov(:,:,1,:));
    end
    nx = size(temp_mov,2);
    ny = size(temp_mov,1);
    white_noise_movie = cat(3,mov{1},mov{2},mov{3},mov{4},mov{5},mov{6},mov{7},mov{8},mov{9}); %do this smarter
    Stim0 = reshape(white_noise_movie,nx*ny,[]);
    Stim.(ndf_names{j}) = bsxfun(@minus,Stim0,mean(Stim0(:)));
end

load('/Volumes/lab/Experiments/Array/Shared/Ruda/better decoding code/data files/2018-11-30/wninfo_4types.mat','STAS','wn_spikes','nonrep_interleave_starts','nkt','center_fits')
% STAS.(ndf_names{j}){type}{rgc}
% wn_spikes.(ndf_names{j}){type}{rgc}
% center_fits.(ndf_names{j}){type}{rgc}
% nonrep_interleave_starts{j}: use to align the wn_spikes to the white noise movie