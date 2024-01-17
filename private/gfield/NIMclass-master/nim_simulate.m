function [spksR, psth] = nim_simulate( nim, Xstims, Nreps )
% [spksR, psth] = nim_simulate( nim, Xstims, <Nreps> )
%
% Simulates spike trains of the specified model, also returns psth over repeats
%
% INPUTS:
%   nim: model structure
%   Xstims: time-embedded stimulus matrices
%   Nreps: number of repeats to simulate (default=1)
%
% OUTPUTS:
%   spksR: spike times generated from the simulation. Multiple repeats are stored in one list, with 
%          a '-1' separating each repeat.
%   psth: firing rate (in units of spks per time bin) computed from simulated spikes spksR

if ~iscell(Xstims)
	tmp = Xstims;
	clear Xstims
	Xstims{1} = tmp;
end
NT = size(Xstims{1},1);

% Calculate non-spike-history contributions to model internals (same on every repeat)
G = nim.process_stimulus( Xstims, 1:length(nim.subunits), [] );

spkhstlen = nim.spk_hist.spkhstlen;
dt = nim.stim_params.dt;
spksR = [];

if spkhstlen == 0
	% then identical rate on all trials
	pred_rate = nim.apply_spkNL( G + nim.spkNL.theta ); % apply spiking NL
	% Poisson generator for all spikes
	
	RobsR = rand(NT,Nreps) < pred_rate*ones(1,Nreps);
	
	% Generate spk time list 
	for nn=1:Nreps
		spksR = [spksR; find(RobsR(:,nn) == 1)*dt-dt/2; -1;];
	end
	
else
	% construct spike-history term
	h = zeros(1,spkhstlen); % spike-history term
	Lh = nim.spk_hist.bin_edges(end)-1;
	for n = 1:spkhstlen
		h(nim.spk_hist.bin_edges(n):(nim.spk_hist.bin_edges(n+1)-1)) = nim.spk_hist.coefs(n);
	end
	
	RobsR_buff = zeros(NT+Lh,Nreps);  % add buffer at beginning for spike history
	Gbuff = [zeros(Lh,1); G;];
	for t = (1:NT)+Lh
		Gspkhist = ones(1,Nreps)*Gbuff(t) + h * RobsR_buff(t-(1:Lh),:);
		rs = nim.apply_spkNL( Gspkhist + nim.spkNL.theta ); % apply spiking NL
		%RobsR_buff(t+spkhstlen,:) = rand(1,Nreps) < rs;
		RobsR_buff(t,:) = rand(1,Nreps) < rs;
	end
	
	% Generate spk time list 
	for nn = 1:Nreps
		spksR = [spksR; (find(RobsR_buff(Lh+1:end,nn) == 1)-0.5)*dt; -1;];
	end
end

fprintf( 'Simulated %d repeats. Average firing rate = %0.2f Hz.\n', Nreps, (length(spksR)-Nreps)/Nreps/(NT*dt) )
if Nreps == 1
	spksR = spksR(1:end-1); % take the -1 off the end if only one rep
end	
	
% Generate PSTH while at it
psth = histc( spksR, (0:NT-1)*dt )/Nreps;

