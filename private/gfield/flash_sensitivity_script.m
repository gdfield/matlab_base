addpath(genpath('/Users/gfield/Dropbox/Code/rod-noise/ShortFlashSingles'))
addpath(genpath('/Users/gfield/Dropbox/Code/rod-noise/GeneralAnalysis'))


cd /Users/gfield/Dropbox/Code/rod-noise/ShortFlashSingles/CurrentRodNoiseAnalysis/

load PCAInfo
load CompleteConcatData
load IsoSingles
SinglesCondition = CompleteConcatData.SinglesCondition;
FailuresCondition = CompleteConcatData.FailuresCondition;
NoiseCondition = CompleteConcatData.NoiseCondition;
clear CompleteConcatData
SinglesPCs = PCAInfo.SinglesPCs;
NumDim = 5;
LocalFlag = 1;
TimingLocalFlag = 0;
%%%%%%%%%%%%%%%%%%%%%%

% Rotate PCs
% The the conversion factor between the average interpolated time and real time
for cells = 1:6
    [temp, TempMinIndex] = min(IsoSingles(cells).SinglesConditionA.AverageResponse);
    MinIndex(cells) = TempMinIndex;
end
AveMinIndex = mean(MinIndex);
[temp, InterpolatedMin] = min(SinglesCondition.AverageResponse);

% actual time [in vector] times conversion factor 
TimeShifts = floor([1] .* (InterpolatedMin/AveMinIndex));    
if (LocalFlag)
    TimeShifts = floor([400] .* (InterpolatedMin/AveMinIndex));
end
if (TimingLocalFlag)        
    TimeShifts = floor([60 30 15] .* (InterpolatedMin/AveMinIndex));    
end
NumShifts = length(TimeShifts);
[NoiseData, temp] = GetGoodEpochData([],NoiseCondition);
Verbose = 0;
for shift = 1:NumShifts
	
	clear TempPCs MeanShiftedSingle
	TempPCs(TimeShifts(shift) + 1:2000,:) = SinglesPCs(1:2000 - TimeShifts(shift),:);
    MeanShiftedSingle(TimeShifts(shift) + 1:2000) = PCAInfo.MeanSingle(1:2000 - TimeShifts(shift));

	if Verbose
        figure(1)
		for cnt = 1:NumDim
            if cnt == 1
                subplot(NumDim, 1, 1)
                plot(MeanShiftedSingle, 'r')
                title('shifted single photon response')
            else
                subplot(NumDim, 1, cnt)
                plot(TempPCs(:, cnt))
            end
		end
		pause(1)
    end
	
    %Filter time shifted PCs to avoid discontinuities
    TempPCs = ComponentFilter(TempPCs, 6);
    MeanShiftedSingle = ComponentFilter(MeanShiftedSingle', 6);
    
    % check to make sure the basic shape of the components has not been
    % changed by the filtering -- sanity check on filtering.
    if Verbose
        figure(2)
		for cnt = 1:NumDim
            if cnt == 1
                subplot(NumDim, 1, 1)
                plot(MeanShiftedSingle, 'r')
                title('shifted single photon response post filtering')
            else
                subplot(NumDim, 1, cnt)
                plot(TempPCs(:, cnt))
            end
		end
		pause(1)
    end

	RotatedPCs(shift).PCs = TempPCs;
    RotatedPCs(shift).MeanSingle = MeanShiftedSingle';
    ContNoise_ShiftedSpace = NoiseData * TempPCs;
    RotatedPCs(shift).ContNoiseCov = cov(ContNoise_ShiftedSpace);
    
    % Compute and store an explicit difference of means discriminator
    RotatedPCs(shift).DifferenceOfMeans = PCAInfo.MeanSingle - MeanShiftedSingle';
    BiasThreshold = ((PCAInfo.MeanSingle * PCAInfo.MeanSingle') - (MeanShiftedSingle' * MeanShiftedSingle)) .* 0.5;
    RotatedPCs(shift).BiasThreshold = BiasThreshold;
end

PCAInfo.RotatedPCs = RotatedPCs;
PCAInfo.TimeShifts = TimeShifts;

clear IsoSingles

%% Setup parameters of model

nl_thresh = 0; % set to one for nonlinear pooling, 0 for linear
resp_shift = 0;
%NumResponses = round(0.024 *8* 500000); 
NumResponses = 48;

NoiseParameters.PoissonFlag  = 1;
NoiseParameters.SinglesVarScale = 1.37^2;
NoiseParameters.AmpVar = 0;
NoiseParameters.ThermalRate = 0.0035;
NoiseParameters.ContNoiseScale = 1.22^2;

Verbose = 0;
FlashStrengths = [0.001 0.003 0.01 0.03 0.05 0.07 0.1 0.2 0.3 0.6 1.0];
num_FlashStrengths = length(FlashStrengths);
NumTrials = 1000;

%% Generate discriminant
TrainingNoiseParameters = NoiseParameters;
TrainingNoiseParameters.ThermalRate = 0; % set thermal rate to zero for training
% simulate single photon responses
[TrainConditionA, ~] = RotatedPCSim(PCAInfo, 1, 5000, 1, TrainingNoiseParameters, Verbose);
% simulate continuous noise
[TrainConditionB, ~] = RotatedPCSim(PCAInfo, 1, 5000, 0, TrainingNoiseParameters, Verbose);

A_traindata = TrainConditionA.EpochData.Data;
A_traindata = A_traindata - repmat(mean(A_traindata,2), [1, size(A_traindata,2)]);
B_traindata = TrainConditionB.EpochData.Data;
B_traindata = B_traindata - repmat(mean(B_traindata,2), [1, size(B_traindata,2)]);

mean_train_A = mean(A_traindata);
mean_train_B = mean(B_traindata);

% compute difference of means (DOM) discriminant
DOM_discrim = mean_train_A - mean_train_B;
norm_DOM = norm(DOM_discrim);
DOM_discrim = DOM_discrim./norm(DOM_discrim);
plot(DOM_discrim)

%% calculate NL
[TrainConditionA, ~] = RotatedPCSim(PCAInfo, 1, 5000, 1, TrainingNoiseParameters, Verbose);
[TrainConditionB, ~] = RotatedPCSim(PCAInfo, 1, 5000, 0, TrainingNoiseParameters, Verbose);
A_traindata = TrainConditionA.EpochData.Data;
A_traindata = A_traindata - repmat(mean(A_traindata,2), [1, size(A_traindata,2)]);
B_traindata = TrainConditionB.EpochData.Data;
B_traindata = B_traindata - repmat(mean(B_traindata,2), [1, size(B_traindata,2)]);

Atrain_projs = A_traindata * DOM_discrim';
Btrain_projs = B_traindata * DOM_discrim';
 
% get mean and SD of signal and noise projected through discriminant.
mean_A_projs = mean(Atrain_projs);
mean_B_projs = mean(Btrain_projs);
std_A_projs = std(Atrain_projs);
std_B_projs = std(Btrain_projs);

bin_steps = -10:0.1:20;
A_dist = 0.0001*normpdf(bin_steps, mean_A_projs,std_A_projs);
B_dist = 0.9999*normpdf(bin_steps, mean_B_projs,std_B_projs);
semilogy(bin_steps, A_dist, 'r', bin_steps, B_dist, 'k')


%%

PCorrects = zeros(1,num_FlashStrengths);
for fs = 1:num_FlashStrengths
    disp(['flash strength is ', num2str(FlashStrengths(fs))])
    temp_correct = 0;
     
    for tl = 1:NumTrials  
        
        if mod(tl, 10) == 1
            fprintf('*')
        end
        
        % sample test data
        TestConditionA = SimplePCSim(PCAInfo, 1, NumResponses, FlashStrengths(fs), NoiseParameters, Verbose);
        TestConditionB = SimplePCSim(PCAInfo, 1, NumResponses, 0, NoiseParameters, Verbose);
        
        Adata = TestConditionA.EpochData.Data;
        Adata = Adata - repmat(mean(Adata,2), [1, size(Adata,2)]);
        Bdata = TestConditionB.EpochData.Data;
        Bdata = Bdata - repmat(mean(Bdata,2), [1, size(Bdata,2)]);
              
        % project responses along discriminant
        Aproj = Adata * DOM_discrim';
        Bproj = Bdata * DOM_discrim';
        
        % apply nonlinearity if flag is true
        if nl_thresh
 
            A_A_likelihoods = normpdf(Aproj, mean_A_projs, std_A_projs) .* poisspdf(1,FlashStrengths(fs));
            A_B_likelihoods = normpdf(Aproj, mean_B_projs, std_B_projs) .* poisspdf(0,FlashStrengths(fs));
            %A_A_likelihoods = normpdf(Aproj, mean_A_projs, std_A_projs) .* poisspdf(1,0.001);
            %A_B_likelihoods = normpdf(Aproj, mean_B_projs, std_B_projs) .* poisspdf(0,0.001);
            
            B_B_likelihoods = normpdf(Bproj, mean_B_projs, std_B_projs) .* poisspdf(0,FlashStrengths(fs));
            B_A_likelihoods = normpdf(Bproj, mean_A_projs, std_A_projs) .* poisspdf(1,FlashStrengths(fs));
            %B_B_likelihoods = normpdf(Bproj, mean_B_projs, std_B_projs) .* poisspdf(0,0.001);
            %B_A_likelihoods = normpdf(Bproj, mean_A_projs, std_A_projs) .* poisspdf(1,0.001);
            
            A_photon_LR = A_A_likelihoods./A_B_likelihoods;
            B_photon_LR = B_A_likelihoods./B_B_likelihoods;

            Aproj = Aproj .* A_photon_LR;
            Bproj = Bproj .* B_photon_LR;
        end
        
        %% Pooling step
        totalA = sum(Aproj);
        totalB = sum(Bproj);

        if totalA > totalB
            temp_correct = temp_correct +1;
        elseif totalA == totalB
            tmp_rnd = rand(1);
            if tmp_rnd > 0.5
                temp_correct = temp_correct +1;
            end
        end
    end
    PCorrects(fs) = temp_correct./NumTrials;
    disp(num2str(temp_correct./NumTrials))
end

%% Plot and fit
figure(1); clf;
semilogx(FlashStrengths, PCorrects, 'r*')
hold on

coef = 0.003;
fitcoef = nlinfit(FlashStrengths, PCorrects, 'cumulative_gaussian', coef);
FitCorrect = cumulative_gaussian(fitcoef, FlashStrengths');
DetectThresh = norminv(0.82, 0, abs(fitcoef));

semilogx(FlashStrengths, FitCorrect, 'r')

pcorrect_ideal = 1 - (exp(-1*NumResponses * FlashStrengths)./2);
semilogx(FlashStrengths, pcorrect_ideal, '--k')

hold off
title_string = ['thresh = ',num2str(DetectThresh), ', rod# = ', num2str(NumResponses)];
title(title_string)
xlabel('Rh*/rod')
ylabel('Probability Correct')
axis([0.001 1 0.4 1])
axis square


print(1, '~/Desktop/rod_current_Lin_48.pdf', '-dpdf')


%%

sim_struct.FlashStrenths = FlashStrengths;
sim_struct.PCorrect = PCorrects;
sim_struct.NoiseParams = NoiseParameters;
sim_struct.NumTrials = NumTrials;
sim_struct.RF_size = NumResponses;
sim_struct.NLdiscrim = nl_thresh;
sim_struct.DetectThresh = DetectThresh;
sim_struct.fitcoef = fitcoef;
sim_struct.fit_type = 'cumulative Gaussian';

cd ~/Desktop
save Lin-48_rods sim_struct




