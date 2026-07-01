function [ dsAll_id] = DsCellFinder(DataBlock, DB, Params) 
% params=params.SkipBWMap

% JC 1/26/2016
% this function will find the DS cells in a population and provide the
% cell_ids within the BW data set

TrialTrigInterval = 10 ; % sec

% find ds cells with drifting gratings data

% load data 

if isfield(Params,'TimeBounds') ; % if you are selecting only part of dataRun
    dataRun = load_data(DataBlock(DB).DsConcatPath) ;
    dataRun = load_neurons(dataRun) ;
%     dataRun = load_ei(dataRun, 'all') ;
    dataRun.triggers  = dataRun.triggers(dataRun.triggers>=Params.TimeBounds(1) & dataRun.triggers<=Params.TimeBounds(2)) ; 
else
    dataRun = load_data(DataBlock(DB).DsPath) ;
    dataRun = load_neurons(dataRun) ;
%     dataRun = load_ei(dataRun, 'all') ;
end
%36 or 24 or 30
%dataRun.names.stimulus_path = [DataBlock(DB).DsPath(1:end-24),'/stimuli/s',DataBlock(DB).DsPath(end-1:end),'.txt'] ;

dataRun.names.stimulus_path = '/Users/amrudzite/Desktop/SC_grant/R21-10/LE/2021-07-13-1/stimuli/s01.txt' ;
dataRun = load_stim(dataRun,'user_defined_trigger_interval', TrialTrigInterval) ;

% stats
totaldur = 0 ;
cell_ids = dataRun.cell_ids ;
[NumSpikesCell, StimComb] = get_spikescellstim(dataRun,cell_ids,totaldur) ;

[mag dsindex magmax magave angle rho theta num U V spave] = dscellanalysis(NumSpikesCell, StimComb) ;

% get mag and angle for each cell at TP with highest dsi
a=1;b=2;
for rgc=1:length(dsindex{a}) ;
    [m,mi] = max([dsindex{a}(rgc),dsindex{b}(rgc)]) ; %
    mag_maxdsi(rgc) = mag{mi}(rgc) ;
    angle_maxdsi(rgc) = angle{mi}(rgc) ;

    angle_maxdsi_deg(rgc) = angle_maxdsi(rgc)*180/pi ;
    if angle_maxdsi_deg(rgc)<0 ;
        angle_maxdsi_deg(rgc) = angle_maxdsi_deg(rgc)+360 ;
    end   
end

SelectionFlag = 1 ;
while SelectionFlag==1 ; % keep going if not satisfied

    figure; clf ; 
    subplot(1,2,1)
    plot(log(dsindex{a}), log(dsindex{b}), '*');%xlim([-11 0])
    hold on
    %plot([0,1],[0,1],'k-')
%     axis tight
    title('ndf0')
%     xlabel(['TP ',num2str(num(1))])
%     ylabel(['TP ',num2str(num(2))])
    
    
    
%     figure;plot(log(dsindex{1}),log(dsindex{2}),'*');xlim([-11 0])
%     figure;plot(log(dsindex{1}),log(dsindex{3}),'*');xlim([-11 0])
%     figure;plot(log(dsindex{1}),log(dsindex{4}),'*');xlim([-11 0])
%     figure;plot(log(dsindex{2}),log(dsindex{3}),'*');
%     figure;plot(log(dsindex{2}),log(dsindex{4}),'*');
%     figure;plot(log(dsindex{3}),log(dsindex{4}),'*');

    disp('select ds set') ;
    [x, y] = ginput;
    plot(x, y);
    IN = inpolygon(log(dsindex{a}), log(dsindex{b}), x, y);
    I = find(IN == 1);
    plot(log(dsindex{a}(I)), log(dsindex{b}(I)), 'ro')

    dsAll_id = dataRun.cell_ids(I);

    % show directions 
    subplot(1,2,2)
    for rgc=1:length(I) ;
        polar([0, angle_maxdsi(I(rgc))],[0,ones(1,length(I(rgc)))],'r')
        hold on
        axis tight
    end

    magRatio = mag{1}./mag{2} ;
    
    % allow subselection of ds cells 
%     subplot(1,5,3)
%     disp('select ds subset') ;
%     plot(zeros(1,length(I)),log10(magRatio(I)),'r*') 
%     hold on
%     plot([-1,1],[0,0],'k')
%     axis tight
%     [x, y] = ginput;
%     plot(x, y);
%     IN = inpolygon(zeros(1,length(I)),log10(magRatio(I)), x, y);
%     Isub{1} = I(find(IN == 1));
%     plot(zeros(1,length(Isub{1})),log10(magRatio(Isub{1})), 'go')
% 
%     Isub{2} = setxor(I,Isub{1}) ; % those in I not in subset
% 
%     % directions (subset)
%     if ~isempty(Isub{1}) ;
%         subplot(1,5,2)
%         for rgc=1:length(Isub{1}) ;
%             polar([0, angle_maxdsi(Isub{1}(rgc))],[0,ones(1,length(Isub{1}(rgc)))],'g')
%         end
%     end

%     % direction selections
%     for subi=1:2 ; % for each subset
%         if ~isempty(Isub{subi}) ;
%             subplot(1,5,3+subi)
%             for rgc=1:length(Isub{subi}) ;
%                 polar([0, angle_maxdsi(Isub{subi}(rgc))],[0,ones(1,length(Isub{subi}(rgc)))],'k')
%                 hold on
%             end
% 
%             for drSelect=1:4 ;
%                 disp('mark direction boundaries') ;
%                 [x,y] = ginput ;
%                 Atemp = atan2d(y,x) ;
%                 Atemp(Atemp<0) = 360 + Atemp(Atemp<0) ;
%                 if max(Atemp)-min(Atemp)>acuteAngle(Atemp(1),Atemp(2))+2 ; % if its on the border
%                     Isub_dirGroup{subi}{drSelect} = Isub{subi}(angle_maxdsi_deg(Isub{subi})>max(Atemp) | angle_maxdsi_deg(Isub{subi})<min(Atemp)) ;
%                 else
%                     Isub_dirGroup{subi}{drSelect} = Isub{subi}((angle_maxdsi_deg(Isub{subi})<max(Atemp)& angle_maxdsi_deg(Isub{subi})>min(Atemp))==1) ;
%                 end
% 
%                 for rgc=Isub_dirGroup{subi}{drSelect}' ;
%                     subplot(1,5,1)
%                     plot(log(dsindex{1}(rgc)), log(dsindex{2}(rgc)), 'r*')
%                     
%                     subplot(1,5,3+subi)
%                     polar([0, angle_maxdsi(rgc)],[0,1],'r')
%                 end
% 
%                 ds_id{subi}{drSelect} = dataRun.cell_ids(Isub_dirGroup{subi}{drSelect});
%                 dsName{subi}{drSelect} = input('name cell type','s') ;
%             end
%         else
%             for drSelect=1:4 ;
%                 ds_id{subi}{drSelect} = [] ;
%                 dsName{subi}{drSelect} = [] ;
%             end    
%         end
%     end
%     
    SelectionFlag = input('enter 1 to repeat, 0 to finish') ; 
    if SelectionFlag==1 ;
        clearvars -except DataBlock DB Params dataRun SelectionFlag...
            mag dsindex magmax magave angle rho theta num U V spave NumSpikesCell StimComb...
             mag_maxdsi angle_maxdsi angle_maxdsi_deg psth_maxDivMean_Av psth_maxDivMean psth
    end
end
% 
if ~isfield(Params,'SkipBwMap')

    % correlate with ids in BW dataset
    dataRunMaster = load_data(DataBlock(DB).BwPath{1}) ;
    dataRunMaster = load_neurons(dataRunMaster) ;
    dataRunMaster = load_ei(dataRunMaster, 'all') ;

    % map using electrical images
    cell_list_map = map_ei(dataRunMaster, dataRun) ;

    for c=1:length(cell_list_map) ;
        if isempty(cell_list_map{c}) ;
            cell_list_map_mat(c) = nan ;
        else
            cell_list_map_mat(c) = cell_list_map{c} ;
        end
    end

    % cells ids in slave for each UniqueCellType set in master data
    MasterIds{1}=cell(1,4) ;
    MasterIds{2}=cell(1,4) ;
    for subi=1:2 ; % for each subset
        for drSelect=1:4 ; 
            if ~isempty(ds_id{subi}{drSelect})
                for c=1:length(ds_id{subi}{drSelect}) ;
                    ci = find(cell_list_map_mat==ds_id{subi}{drSelect}(c)) ;
                    if ~isempty(ci) ;
                        MasterIds{subi}{drSelect} = [MasterIds{subi}{drSelect},dataRunMaster.cell_ids(ci)] ;
                    end
                end
            end
        end
    end
else
    MasterIds=[] ;
end


% for igor
% ForIgor.ds_id = ds_id ;
% ForIgor.dsName = dsName ;
% ForIgor.master_id = MasterIds ;


% old code

% % get psth stats
% for c=1:length(dataRun.spikes) ; % for each cell 
%     for st=1:length(dataRun.stimulus.triggers) ; % for each grating set
%         triggs = dataRun.triggers(dataRun.triggers>=dataRun.stimulus.triggers(st)...
%             & dataRun.triggers<dataRun.stimulus.triggers(st)+TrialTrigInterval) ;
%         [psth{c}{st},psthTime] = get_psth(dataRun.spikes{c}, triggs) ; %psth
%         [powerspec_xvalues, mean_powerspec] = PowerSpectrumFinder(psth{c}{st},dataRun.sampling_rate) ;
%         psth_varF1(c,st) = mean_powerspec(2)/sum(mean_powerspec) ;
%         psth_varF2(c,st) = mean_powerspec(3)/sum(mean_powerspec) ;
%         psth_maxDivMean(c,st) = max(psth{c}{st})/mean(psth{c}{st}) ; % psth max/mean
%     end
% end
% 
% trial_list = unique(dataRun.stimulus.trial_list) ;
% for st = 1:length(trial_list) ; % for each stim type
%     tl = find(trial_list(st)==dataRun.stimulus.trial_list) ;
%     psth_varF2_mean(:,st) = nanmean(psth_varF2(:,tl),2) ;
% end
% psth_varF2_min = min(psth_varF2_mean,[],2) ; % find lowest f2 for all trials and stim
% 

% % average across all stim presentations and organize psth stats by tp, (like mag, rho, theta, ect.)
% psth_maxDivMean_Av = cell(1,2) ; % average across all stim presentations
% psth_varF1_Av = cell(1,2) ;
% psth_varF2_Av = cell(1,2) ;
% 
% trial_list = unique(dataRun.stimulus.trial_list) ;
% for st = 1:length(trial_list) ; % for each stim type
%     tl = find(trial_list(st)==dataRun.stimulus.trial_list) ;
%     if StimComb(st,2)==min(StimComb(:,2)) ; % if its the lowest temporal period
%         TpGrp = 1 ;
%     else
%         TpGrp = 2 ;
%     end
%         
%     psth_maxDivMean_Av{TpGrp} = [psth_maxDivMean_Av{TpGrp},mean(psth_maxDivMean(:,tl),2)] ; % average across all stim presentations
%     psth_varF1_Av{TpGrp} = [psth_varF1_Av{TpGrp},mean(psth_varF1(:,tl),2)] ;
%     psth_varF2_Av{TpGrp} = [psth_varF2_Av{TpGrp},mean(psth_varF2(:,tl),2)] ; 
% end


