a = 2;
lambda = 5;
X = 0:1:100;
num_samples = 10000;

P_x = zeros(length(X), 1);

for it = 1:length(X)

    P_x(it) = (lambda * a).^(X(it)./a) .* exp(-lambda * a) ./ factorial(X(it) ./ a);
end

% make sure the probabilities sum to 1 (normalized distribution)
P_x = P_x ./ sum(P_x);
rnd_samples = discretesample(P_x, num_samples)-1;


figure(1); clf
plot(X, P_x, 'k');
hold on
[hist_vals, temp_bins] = hist(rnd_samples, [0:1:100]);
plot(temp_bins, hist_vals ./ num_samples, 'r')

%% Poisson

mu = 5;
X = 0:1:20;
num_samples = 10000;
P_x = poisspdf(X, mu);

%% Binomial (sub-Poisson)

N = 10; %smaller numbers produce distributions with less variance
P = mu/N;

X = 0:1:20;

B_x = binopdf(X, N, P);
plot(X, P_x, 'b', X, B_x, 'r')

rnd_samples = binornd(N,P, num_samples,1);
mean(rnd_samples) % mean of random samples from binomial with N and P
var(rnd_samples) % variance of random samples from binomial with N and P

%% Neg Binomial

R = 20;
P = mu ./ (R+mu);

NB_x = nbinpdf(X, R, P);


figure(1); clf;
plot(X, P_x, 'b', X, NB_x, 'r')

rnd_samples = nbinrnd(R, P, num_samples, 1);
mean(rnd_samples)
var(rnd_samples)

P*R./(1-P)







