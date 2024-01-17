load datarun_WT 
load datarun_5M

num_frames = length(datarun_WT{1}.stas.tc_pcs(:,1));

%%
pc_dim = 1;
figure(pc_dim); clf;
plot(datarun_WT{1}.stas.tc_pcs(:,pc_dim), 'k') 
hold on
plot(datarun_WT{2}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{3}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{4}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_5M{1}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{2}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{3}.stas.tc_pcs(:,pc_dim), 'r') 

wt_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_WT)
    wt_pc = wt_pc + datarun_WT{dset}.stas.tc_pcs(:,pc_dim);
end
wt_pc = wt_pc ./ length(datarun_WT);

m5_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_5M)
    m5_pc = m5_pc + datarun_5M{dset}.stas.tc_pcs(:,pc_dim);
end
m5_pc = m5_pc ./ length(datarun_5M);

plot(wt_pc, 'k', 'LineWidth', 3) 
plot(m5_pc, 'r', 'LineWidth', 3) 
hold off
title(['PC ', num2str(pc_dim)])
axis square
axis tight
print(pc_dim, '~/Desktop/pc1.pdf', '-dpdf')


%%
pc_dim = 2;
figure(pc_dim); clf;
plot(datarun_WT{1}.stas.tc_pcs(:,pc_dim), 'k') 
hold on
plot(datarun_WT{2}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{3}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{4}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_5M{1}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{2}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{3}.stas.tc_pcs(:,pc_dim), 'r') 

wt_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_WT)
    wt_pc = wt_pc + datarun_WT{dset}.stas.tc_pcs(:,pc_dim);
end
wt_pc = wt_pc ./ length(datarun_WT);

m5_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_5M)
    m5_pc = m5_pc + datarun_5M{dset}.stas.tc_pcs(:,pc_dim);
end
m5_pc = m5_pc ./ length(datarun_5M);

plot(wt_pc, 'k', 'LineWidth', 3) 
plot(m5_pc, 'r', 'LineWidth', 3) 
hold off
title(['PC ', num2str(pc_dim)])
axis square
axis tight
print(pc_dim, '~/Desktop/pc2.pdf', '-dpdf')

%%
pc_dim = 3;
figure(pc_dim); clf;
plot(datarun_WT{1}.stas.tc_pcs(:,pc_dim), 'k') 
hold on
plot(datarun_WT{2}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{3}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{4}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_5M{1}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{2}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{3}.stas.tc_pcs(:,pc_dim), 'r') 

wt_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_WT)
    wt_pc = wt_pc + datarun_WT{dset}.stas.tc_pcs(:,pc_dim);
end
wt_pc = wt_pc ./ length(datarun_WT);

m5_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_5M)
    m5_pc = m5_pc + datarun_5M{dset}.stas.tc_pcs(:,pc_dim);
end
m5_pc = m5_pc ./ length(datarun_5M);

plot(wt_pc, 'k', 'LineWidth', 3) 
plot(m5_pc, 'r', 'LineWidth', 3) 
hold off
title(['PC ', num2str(pc_dim)])
axis square
axis tight
print(pc_dim, '~/Desktop/pc3.pdf', '-dpdf')


%%
pc_dim = 4;
figure(pc_dim); clf;
plot(datarun_WT{1}.stas.tc_pcs(:,pc_dim), 'k') 
hold on
plot(datarun_WT{2}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{3}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_WT{4}.stas.tc_pcs(:,pc_dim), 'k') 
plot(datarun_5M{1}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{2}.stas.tc_pcs(:,pc_dim), 'r') 
plot(datarun_5M{3}.stas.tc_pcs(:,pc_dim), 'r') 

wt_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_WT)
    wt_pc = wt_pc + datarun_WT{dset}.stas.tc_pcs(:,pc_dim);
end
wt_pc = wt_pc ./ length(datarun_WT);

m5_pc = zeros(num_frames, 1);
for dset = 1:length(datarun_5M)
    m5_pc = m5_pc + datarun_5M{dset}.stas.tc_pcs(:,pc_dim);
end
m5_pc = m5_pc ./ length(datarun_5M);

plot(wt_pc, 'k', 'LineWidth', 3) 
plot(m5_pc, 'r', 'LineWidth', 3) 
hold off
title(['PC ', num2str(pc_dim)])
axis square
axis tight
print(pc_dim, '~/Desktop/pc4.pdf', '-dpdf')
