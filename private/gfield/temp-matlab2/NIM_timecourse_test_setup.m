
% compare GLM fits of time courses to vision STA time courses


%% load stuff
%base_path='/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-06-24-0';
base_path='/Volumes/Disk2/Data/2019-06-24-0';
types_interest={'off brisk sustained';'off brisk transient';'off transient';'on brisk sustained';'on brisk transient';'on transient'};
types_interest=genvarname(types_interest);
%data_path={'/data000-map_MR/data000-map_MR';'/data001-map_MR/data001-map_MR';'/data002-map_MR/data002-map_MR';'/data003-map_MR/data003-map_MR';...
%    '/data004-map_MR/data004-map_MR';'/data005_MR/data005_MR';}; %master last
data_path={'/data000-map_MR/data000-map_MR';'/data005_MR/data005_MR';}; %master last
ndf_names={'NDF5';'NDF0'};
wnruns=[1:6];%indices of wnruns
dataruns=load_all_dataruns(base_path,data_path,ndf_names,wnruns);
[IDZ_run,IDZ_type] = mapped_cells_organize_per_datarun(dataruns,ndf_names,types_interest);

% fix ids after looking at eis
IDZ_run.(types_interest{3}){1} = [287,1051,2536,2793,3226,3437,3498,4307,4743,4831];
% IDZ_run.(types_interest{3}){2} = [287,1051,1276,2116,2161,2536,2793,3226,3437,3498,4039,4307,4743,4831,5746,5748];
% IDZ_run.(types_interest{3}){3} = [737,1051,1276,2116,2161,2536,2793,3226,3437,3498,4039,4307,4564,4743,4831,5746,5748];
% IDZ_run.(types_interest{3}){4} = [287,737,1051,1276,2116,2161,2536,2793,3226,3437,3498,4039,4307,4564,4743,4831,5746,5748];
% IDZ_run.(types_interest{3}){5} = [737,1051,1276,1696,2116,2161,2536,2793,3226,3437,3498,4307,4564,4743,4831,5746,5748];
IDZ_run.(types_interest{3}){2} = [287,737,1051,1276,1426,1696,2116,2161,2536,2793,3226,3437,3498,4039,4307,4564,4743,4831,5746,5748];

IDZ_run.(types_interest{4}){1} = [466,871,1351,3587,3661,5536,6376,7576];
IDZ_run.(types_interest{4}){2} = [331,871,1351,2041,3587,4966,5386,5536,6376,7576];
IDZ_run.(types_interest{5}){2} = [181,227,811,842,1413,1907,1966,2417,2806,3063,3737,4126,4141,4877,5581,5671,6256,6391,7216,7622];
IDZ_run.(types_interest{5}){3} = [181,227,842,1413,1907,2417,2806,3063,3737,4126,4141,4877,5581,5671,6256,6391,6916,7216,7622];
IDZ_run.(types_interest{6}){1} = [543,872,1261,1997,2568,2597,4533,4697,4846,4938];
IDZ_run.(types_interest{6}){3} = [543,1261,1997,2568,2597,3408,3991,4533,4535,4696,4697,4846,4938,5522,5822,6991];


% make time vectors
for j=1:length(wnruns)
    intervals(j)=dataruns.(ndf_names{j}).stimulus.interval;        
end
nkts = [30 30 30 30 30 30];
time_intervals=intervals/60.35;
for j=1:length(wnruns)
    time_vectors{j}=-(nkts(j)-1)*time_intervals(j):time_intervals(j):0;       
end
load(['/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/Ruda/better decoding code/data files/2019-06-24/tc_data.mat'],'bad_tcs','IDZ_run','tcs','ndf_names','time_vectors','intervals')


%% grab stimuli

nms{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
%nms{5} = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
%nms{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
%nms{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-30-1-0.48-11111-26x20-60.35_xoffset.xml';
%nms{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-30-2-0.48-11111-26x20-60.35_xoffset.xml';
nms{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/acquisition/movie-xml/BW-60-4-0.48-11111-13x10-60.35_xoffset.xml';
%trs = [4 2 1 1 1 1];
trs = [4 1];
% gather white noise movie
for j = 1:length(nms)
    tr = trs(j);
    trigs = dataruns.(ndf_names{j}).triggers;
    [this_mov,~,~,dur,refresh(j)] = get_movie(nms{j},trigs,dataruns.(ndf_names{j}).duration*60/tr);
    temp_mov = squeeze(this_mov(:,:,1,:));
    nx = size(temp_mov,2);
    ny = size(temp_mov,1);
    temp_stim = reshape(temp_mov,nx*ny,[]);
    stim_szs{j} = size(temp_mov);
    Stim.(ndf_names{j}) = bsxfun(@minus,temp_stim,mean(temp_stim(:)));
end

%% fit indpt glms
% need to pare down size around RF

nkt = 30;
for j = [1 6]
    j
    dtStim = refresh(j)/1000;
    dtSp = dtStim/10;
    usf = 10;upsamplefactor=usf;
    slen = size(Stim.(ndf_names{j}),2);
    rlen = slen*upsamplefactor;
    Stim_real = Stim.(ndf_names{j})(:,1:slen)';
    for type = 1:6
        type
        for rgc = 1:size(tcs.(types_interest{type}){j},2)
            if ~isnan(tcs.(types_interest{type}){j}(1))
%                 rgc
                
                id = IDZ_run.(types_interest{type}){j}(rgc);
                ind = get_cell_indices(dataruns.(ndf_names{j}),id);
                this_sta = dataruns.(ndf_names{j}).stas.stas{ind};
                rf = rf_from_sta(this_sta);
                
                cntr=dataruns.(ndf_names{j}).stas.fits{ind}.mean;
                lngth=max(dataruns.(ndf_names{j}).stas.fits{ind}.sd)*2;
                a=floor(cntr(1)-2*lngth);
                b=ceil(cntr(1)+2*lngth);
                c=floor(cntr(2)-2*lngth);
                d=ceil(cntr(2)+2*lngth);
                if a<1;a=1;end
                if a>stim_szs{j}(2);a=stim_szs{j}(2);end
                if b>stim_szs{j}(2);b=stim_szs{j}(2);end
                if c<1;c=1;end
                if d>stim_szs{j}(1);d=stim_szs{j}(1);end
                iixkeep = a:b; %longer (53) dimension
                iiykeep = c:d; %shorter (40) dimension

                %stim
                nky = length(iiykeep); 
                nkx = length(iixkeep);
                Stim3d = reshape(Stim_real',stim_szs{j}(1),stim_szs{j}(2),[]);
                Stim_real3d = Stim3d(iiykeep,iixkeep,:);
                Stim_use = reshape(Stim_real3d,nkx*nky,[])';
                
                %sta
                stafull = STAS.(ndf_names{j}){type}{rgc}';   
                sta3d = reshape(stafull,stim_szs{j}(1),stim_szs{j}(2),nkt);
                stasmall = sta3d(iiykeep,iixkeep,:);
                sta_real = reshape(stasmall,nkx*nky,nkt)';

                % spikes
                these_spikes = dataruns.(ndf_names{j}).spikes{ind};
                [binned_spikes,edges] = histcounts(these_spikes,dataruns.(ndf_names{j}).triggers(1):dtSp: (size(Stim.(ndf_names{j}),2)+100)*dtStim);
                sps = binned_spikes(:,1:rlen)';

                %run model
                nhbasis = 8;hpeakFinal = 0.1;
                nkbasis = 8;
                k_rank = 1; 
                ggi = makeFittingStruct_GLMbi(k_rank,dtStim,dtSp,nkt,nkbasis,sta_real,nhbasis,hpeakFinal);
                ggi.mask = [];
                ggi.sps = sps; 
                ggi.nlfun = @logexp2;
                ggi.ihw = randn(size(ggi.ihw))*1;
                % Do ML estimation of model params
                opts = {'display', 'iter'};
                [gg, fval] = MLfit_GLMbi(ggi,Stim_use,opts); % do ML (requires optimization toolbox)
                indpt_glms.(types_interest{type}){j,rgc} = gg;
                fvals.(types_interest{type}){j,rgc} = fval;
            end
        end
        save(['/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Shared/Ruda/better decoding code/data files/2019-06-24/tc_analysis_glmfits_ndf50_2019_06_24_' ndf_names{j} '.mat'],...
    'STAS','fvals','indpt_glms','type','j');
    end
end




