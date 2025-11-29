v_os_power = [0.0012 0.0166 0.053 0.0064 0.0065 0.011 0.0856 0.0345 0.052 0.1065 0.007 0.0639];

h_os_power = [0.5928 0.3803 0.4794 0.7048 0.6477 0.7572 0.4815 0.75 0.8294 0.854 0.6507 0.5885 0.7198 0.6253 0.6602 0.6284 0.6947 0.5925 0.7748 0.693 0.7336 0.6762 0.6541 0.742 0.7448];


mean_v_os_power = mean(v_os_power);
mean_h_os_power = mean(h_os_power);

v_indices = 2 * ones(1, length(v_os_power));
h_indices = ones(1, length(h_os_power));

figure(1); clf;
semilogy(h_indices, h_os_power, 'o', 'Color', [0.5 0.5 0.5])
hold on
errorbar(1, mean_h_os_power, std(h_os_power) ./ sqrt(length(h_os_power)), 'ko', 'MarkerFaceColor', 'k')

semilogy(v_indices, v_os_power, 'o', 'Color', [0.5 0.5 0.5])
errorbar(2, mean_v_os_power, std(v_os_power) ./ sqrt(length(v_os_power)),'ko', 'MarkerFaceColor', 'k')

axis([0 3 0.001 1])
axis square
ylabel('F_1 / F_0')
xticks([1 2])
xticklabels({'horizontal', 'vertical'})
title('frequency modulation')

cd ~/Desktop/
exportgraphics(gcf, 'freq_mod.pdf', 'ContentType', 'vector')
