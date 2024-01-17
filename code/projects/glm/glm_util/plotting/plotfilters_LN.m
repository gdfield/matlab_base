function plotfilters_LN(fittedLN)
fittedGLM = fittedLN.orig_glm;
dt = fittedGLM.t_bin;
tstim = fittedGLM.bins_per_frame * dt;
% Plots the GLM filters using the output from glm_fit
rows = 1;
plots = 3;
pos = [100 100 200*plots 200];

% Plot the rank 1 stimulus filters
K_time1  = fittedLN.new_time;
K_space1 = fittedLN.orig_glm.linearfilters.Stimulus.space_rk1; 

subplot(rows,plots,1)
set(gca, 'fontsize', 10);
imagesc(K_space1);
colormap gray
title('Space Filter'); axis off

LW = 2;
subplot(rows,plots,2)
set(gca, 'fontsize', 10);
frames    = fittedGLM.linearfilters.Stimulus.frame_shifts;
time_msec = 1000*tstim*frames ;
zeroline = zeros(1,length(time_msec)); 
plot(time_msec, zeroline, 'k-'); hold on
plot(time_msec,K_time1,'linewidth', LW); 
xlim([time_msec(1), time_msec(end)])
xlabel('msec'); title('Time Filter');
set(gca, 'ytick', [0]); 

LW = 2;
subplot(rows,plots,3)
set(gca, 'fontsize', 10);
plot(fittedLN.new_NL.gen_signal,fittedLN.new_NL.firing_rate,'linewidth', LW); 
xlabel('msec'); title('Nonlinearity');
set(gca, 'ytick', [0]); 


% end NBCoupling
set(gcf, 'Position', pos)

end

    



     