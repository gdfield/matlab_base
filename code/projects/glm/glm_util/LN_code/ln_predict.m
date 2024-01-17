function fittedLN = ln_predict(fittedLN, testmovie, varargin)
% input same things you would input to glm_fit
% If no STA is given, use the linear filter output from the fitted GLM

p = inputParser;
p.addParameter('testspikes', 0)
p.addParameter('predict', true)
p.addParameter('trials', 20)
p.addParameter('PSTH_smoothing', 12); % 10 ms smoothing
p.parse(varargin{:});
testspikes = p.Results.testspikes;
trials = p.Results.trials;
smoothing = p.Results.PSTH_smoothing;
clear p
monitor_refresh = 1/(fittedLN.orig_glm.t_bin*fittedLN.orig_glm.bins_per_frame);
if iscell(testspikes)
    trials = length(testspikes);
end
% new predictions
new_GS = conv(glm_gen_signal_spatial(fittedLN.orig_glm, testmovie), fittedLN.new_time, 'full');%+fittedLN.orig_glm.linearfilters.TonicDrive.Filter;
try
    new_pred = predict(fittedLN.model, new_GS(1:size(testmovie,3))');
catch
    new_pred = feval(fittedLN.model, new_GS(1:size(testmovie,3))');
end
new_pred(new_pred<=0)=0.1;

fittedLN.xval.smoothing = smoothing;
fittedLN.xval.predicted_firing_rate = new_pred;
fittedLN.xval.rasters.ln_sim = Poisson_spiking(new_pred,trials,fittedLN.orig_glm.bins_per_frame,monitor_refresh);

% need testspikes
if iscell(testspikes)
    %firing =  IDP_plot_PSTH(testspikes,1, 'color',0, 'smoothing', 1);
    %idx = 1:min(length(new_pred), length(firing));
    fittedLN.orig_glm.xvalperformance = glm_predict(fittedLN.orig_glm, testmovie,'testspikes', testspikes);
    fittedLN.xval.rasters = fittedLN.orig_glm.xvalperformance.rasters;
    fittedLN.xval.rasters.ln_sim = Poisson_spiking(new_pred,trials,fittedLN.orig_glm.bins_per_frame,monitor_refresh);
    
    % get firing rates from predicted spikes
    convolve = gausswin(smoothing)/(trials*fittedLN.orig_glm.t_bin*sum(gausswin(smoothing)));
    LN_FR = conv(sum(fittedLN.xval.rasters.ln_sim),convolve,'same');
    GLM_FR = conv(sum(fittedLN.xval.rasters.glm_sim),convolve,'same');
    REC_FR = conv(sum(fittedLN.xval.rasters.recorded),convolve,'same');

    % calculate all the metrics
    % Corr
    fittedLN.xval.orig_glm.corr = corr(REC_FR', GLM_FR');
    fittedLN.xval.corr = corr(REC_FR', LN_FR');
    % MSE
    fittedLN.xval.orig_glm.MSE = sum((REC_FR - GLM_FR).^2)/length(REC_FR);
    fittedLN.xval.MSE = sum((REC_FR - LN_FR).^2)/length(REC_FR);
    % FEV
    %%{
    firing_even = conv(sum(fittedLN.xval.rasters.recorded(2:2:end,:)),convolve,'same');
    firing_odd = conv(sum(fittedLN.xval.rasters.recorded(1:2:end,:)),convolve,'same');
    A_LN = sum((REC_FR - LN_FR).^2);
    A_GLM = sum((REC_FR - GLM_FR).^2);
    B = sum((REC_FR-mean(REC_FR)).^2);
    C = sum((firing_even-firing_odd).^2);
    D = sum((firing_even-mean(firing_even)).^2);
    fittedLN.xval.orig_glm.FEV = (1-A_GLM/B)/(1-C/D);
    fittedLN.xval.FEV = (1-A_LN/B)/(1-C/D);
    %}
    
    %%{
    fittedLN.xval.firing.sim = LN_FR;
    fittedLN.xval.firing.orig_glm = GLM_FR;
    fittedLN.xval.firing.rec = REC_FR;
    %}
    

end

end