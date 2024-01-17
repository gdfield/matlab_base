%% center surround example pair
samps = 0:0.01:20;
mu_one = 9;
sig = 1;
mu_two = 11;
surround_strength = 0.5;
sig_sur = 3;

gauss_one = normpdf(samps, mu_one, sig) - surround_strength*normpdf(samps, mu_one, sig_sur);
gauss_two = normpdf(samps, mu_two, sig) - surround_strength*normpdf(samps, mu_two, sig_sur);

gauss_one = gauss_one ./ max(gauss_one);
gauss_two = gauss_two ./ max(gauss_two);

plot(samps, gauss_one, 'k', samps, gauss_two, 'k');
hold on

%% center only example pair
samps = 0:0.01:20;
mu_one = 9;
sig = 1.5;
mu_two = 11;

gauss_one = normpdf(samps, mu_one, sig);
gauss_two = normpdf(samps, mu_two, sig);

gauss_one = gauss_one ./ max(gauss_one);
gauss_two = gauss_two ./ max(gauss_two);

plot(samps, gauss_one, 'g', samps, gauss_two, 'g');

plot(samps, zeros(length(samps)), 'k')

print(1, '~/Desktop/rf_example.pdf', '-dpdf')

hold off; clf;