%% calculate information on a spike train.
% Generate a 'response' as a spike rate
spike_rate_1 = poissrnd(0.2,[100,4]);
spike_rate_2 = poissrnd(0.4, [100,4]);
spike_rate_3 = poissrnd(0.2, [100,4]);

spike_rate = [spike_rate_1, spike_rate_2, spike_rate_3];

%% explicitly compute information from a set of responses
X = [1 0; 1 0; 1 0; 1 1; 1 1; 1 1; 1 1; 1 0; 0 1; 0 1; 0 1; 0 1; 0 0; 0 0; 0 0; 0 0];

word_length = 4;
spike_rate_binned = zeros(100,3);
for w = 1:3
    start_n = 1 + ((word_length * w) - 1);
    end_n = word_length * w;
    spike_rate_binned(:,w) = sum(spike_rate(:, start_n:end_n),2);
end


[num_observations, num_vars] = size(X);
hits_vec = zeros(num_observations,1);
eventn = 0;
frequencies = [];
for obs = 1:num_observations
    if hits_vec(obs) == 0 % check that obsevation hasn't been counted
        countr = 1;
        temp_vec = X(obs,:); % grab observation
        hits_vec(obs) = 1; % mark it counted
        for sub_obs = obs+1:num_observations % loop through other obervations
            if isequal(temp_vec, X(sub_obs,:)) % find identical ones to seed
                countr = countr +1;  % count their frequency 
                hits_vec(sub_obs) = 1; % mark the event as counted
            end
            
        end
        eventn = eventn + 1; % increment the counter on observed events
        frequencies(eventn) = countr ./ num_observations; % calculate the probability
    end
   
end

H = -1.*sum(frequencies .* log2(frequencies))

%% Compute the information 'across' the response 

num_repeats = 200;
epoch_length = 2000;
add_noise_prob = 0.005;
sub_noise_prob = 0.33;
word_size = 4;

% make a fake spike train;
spike_train = poissrnd(0.25,[1,epoch_length]);
temp_inds = find(spike_train >0);
spike_train(temp_inds) = spike_train(temp_inds) -1;

spike_train_matrix = repmat(spike_train, num_repeats, 1);

% find indices in mean spike train that are non-zero;
nz_indices = find(spike_train_matrix > 0);
temp_mat = rand(length(nz_indices),1);
temp_mat(temp_mat <= sub_noise_prob) = -1;
temp_mat(temp_mat > sub_noise_prob) = 0;

% subtract some spikes for noise;
spike_train_matrix(nz_indices) = spike_train_matrix(nz_indices) + temp_mat;

% make_noise
noise_matrix = rand(num_repeats,epoch_length);

noise_matrix(noise_matrix > 1-noise_prob) =1;
noise_matrix(noise_matrix <= 1-noise_prob) =0;

spike_train_matrix = spike_train_matrix + noise_matrix;

figure(1); clf;
hold on
for rep = 1:num_repeats
    temp_spikes = find(spike_train_matrix(rep, :) > 0);
    plot(temp_spikes, rep*ones(length(temp_spikes),1), 'k.')
end
hold off

make_life_underall =     
    
    
    

