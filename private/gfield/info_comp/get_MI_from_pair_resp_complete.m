function [I,csr]=get_MI_from_pair_resp_complete(spikes,bins,word_size,starts,samples)

%inputs: spikes is a 2x1 cell, each containing the spikes of one rgc 
%           organized by bins and trials
%        bins are trial_start:bin_size:trial_stop, for eg 0:.005:10
%        word_size is in number of bins, generally between 1 and 4
%        starts is a vector of start times for each trial
%        samples is a vector of how to subsample trials, for eg [1 1/2 1/4]. 
%           note: each subsample must divide up the total number of trials evenly

%outputs: I is information in bits
%         csr (correction size ratio) is the ratio of the naive info estimate with the actual info, 
%           meant for tracking the correction size


% process: for each sample size, get MI. average across sub samples. fit
% for sample sizes and get infinity extrapolation.

num_trials=length(starts);

for s=1:length(samples)
    runs=1/samples(s);
    trials_per_run=num_trials/runs;
    %split up spikes
    for r=1:runs
        subsampled_spikes{r}{1}=spikes{1}(1+(r-1)*trials_per_run:trials_per_run*r,:);
        subsampled_spikes{r}{2}=spikes{2}(1+(r-1)*trials_per_run:trials_per_run*r,:);
    end
    for r=1:runs
        subsampled_MI(r)=get_MI_from_pair_resp(subsampled_spikes{r},bins,word_size,starts(1+(r-1)*trials_per_run:trials_per_run*r));
    end
    avg_MI=mean(subsampled_MI);
    MI(s)=avg_MI;
    subsampled_sds(s)=std(subsampled_MI);  
end

%fit 
% figure;plot(samples,MI,'o')
model=@(s0,a,b,x) s0+a./x+b./x.^2;
f=fit(samples',MI',model);
% hold on
% plot(f)
% hold on
% plot(1.1,f.s0,'g*')

I=f.s0;
csr=(MI(1)-I)/I;



end