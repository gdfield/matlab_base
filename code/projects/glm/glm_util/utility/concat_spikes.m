function fitspikes = concat_spikes(spikes_cell, block_length)
% takes in spikes in cells and puts them back to back. BLOCK_LENGTH should
% be in seconds
fitspikes = [];
for block = 1:length(spikes_cell)
    fitspikes = [fitspikes; spikes_cell{block}+block_length*(block-1)];
end
end