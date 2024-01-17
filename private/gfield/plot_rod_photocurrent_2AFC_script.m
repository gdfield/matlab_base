cd ~/Desktop/Rod-Release-Study/rod-photocurrents

load lin-24_rods.mat
lin_pool = sim_struct;

load nl-24_rods.mat
nl_pool = sim_struct;
clear sim_struct

figure(1); clf;
semilogx(lin_pool.FlashStrenths, lin_pool.PCorrect, '-ok');
hold on
semilogx(nl_pool.FlashStrenths, nl_pool.PCorrect, '-ob');

pcorrect_ideal = 1 - (exp(-1*nl_pool.RF_size * nl_pool.FlashStrenths)./2);
semilogx(nl_pool.FlashStrenths, pcorrect_ideal, '--r')
xlabel('Rh*/rod')
ylabel('% correct')
title('2AFC detection on photocurrents')
axis tight
axis([0.001 2 0.4 1])
axis square

print(1, '~/Desktop/photocurrent_24_rods.pdf', '-dpdf')


