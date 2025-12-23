function [OSIs, DSIs, mean_corr] = find_OS_RGCs(datarun, cell_spec, varargin)
% find_OS_RGCs     identify OS cells within a list
%
% usage:  [OSIs, DSI, corr_vals] = find_OS_RGCs(datarun, cell_spec, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    OSIs - OSIs for the entered cell_spec list
%               DSIs -- DSIs for the cell_spec list 
%               median_corr -- indicates how reproducible the responses are
%                           from trial to trial
%
% optional parameters, their default values, and what they specify:
%
%
% plot_tuning           true      	plot rasters for each direction and polar tuning function
%
% 
% 2025-12 gdf: Created
%

p = inputParser;

% specify list of optional parameters
p.addParameter('plot_tuning', true, @islogical)
p.addParameter('plot_tuning_os_only', false, @islogical)
p.addParameter('verbose', true, @islogical);
p.addParameter('OSI_thresh', 0.3, @isnumeric)
p.addParameter('DSI_thresh', 0.25, @isnumeric)
p.addParameter('response_cor_thresh', 0.3, @isnumeric)
p.addParameter('export_plot', false, @isnumeric)
p.addParameter('pause', false, @islogical)
p.addParameter('start_fig_num', 1, @isnumeric)
%p.addParameter('fast_plot', true, @isnumeric)


% resolve user input and default values
p.parse(varargin{:});

% check that stimulus field is there:
if isfield(datarun, 'stimulus')
    % extract spatial and temporal periods from datarun
    spatial_periods = datarun.stimulus.params.SPATIAL_PERIOD;
    temp_periods = datarun.stimulus.params.TEMPORAL_PERIOD;
    num_sps = length(spatial_periods);
    num_tps = length(temp_periods);
    num_dirs = length(datarun.stimulus.params.DIRECTION);
else
    error('datarun does not have a stimulus field; see load_stim function')
end

% get the indices for the cells of interest in cell_spec
cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);

% initialize variable to write outputs to
OSIs = zeros(num_rgcs, num_sps, num_tps);
DSIs = zeros(num_rgcs, num_sps, num_tps);
mean_corr = zeros(num_rgcs, num_sps, num_tps);

for rgc = 1:num_rgcs
    OSI = zeros(num_sps,num_tps);
    DSI = zeros(num_sps, num_tps);
    mean_cor = zeros(num_sps,num_tps);

    % get spikes from an RGC and extract tuning curves
    temp_spike_times = datarun.spikes{cell_indices(rgc)};
    
    for spat_p = 1:num_sps
        for temp_p = 1:num_tps
            [tuning_struct, temp_spike_nums, ~] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                        'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
              
            % compute Orientation-Selective Index
            [ds_max, ds_ind] = max(temp_spike_nums);
            condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs/2, 2]), 2);
            [min_val, min_ind] = min(condensed_tuning);
            
            % handles case of 8 directions
            if num_dirs == 8

                % for OS
                orth_ind = mod(min_ind+2,4);
                if orth_ind == 0
                    orth_ind = 4;
                end
                max_val = condensed_tuning(orth_ind);
                
                % for DS
                op_ind = mod(ds_ind+4, 8);
                if op_ind == 0
                    op_ind = 1;
                end
                null_val = temp_spike_nums(op_ind);
                
            % handles case of 12 directions    
            elseif num_dirs == 12
                
                % for OS
                orth_ind = mod(min_ind+3,6);
                if orth_ind == 0
                    orth_ind = 6;
                end
                max_val = condensed_tuning(orth_ind);
                
                % for DS
                op_ind = mod(ds_ind+6, 12);
                if op_ind == 0
                    op_ind = 1;
                end
                null_val = temp_spike_nums(op_ind);
                
            end
            
            % compute reliability (in the preferred
            % orientation/direction)
            [~, temp_rep_num] = size(tuning_struct);
            temp_bin_size = 0.1;
            temp_bins = 0:temp_bin_size:8;
            trial_hists = zeros(length(temp_bins)-1, temp_rep_num);
            %temp_spk_times = [];
            for rn = 1:temp_rep_num
                if isempty(tuning_struct{orth_ind, rn})
                    trial_hists(:, rn) = 0;
                else
                    temp_trial_hists = histcounts(tuning_struct{orth_ind, rn},temp_bins);                    
                    trial_hists(:,rn) = temp_trial_hists ./ norm(temp_trial_hists);
                end
            end

            mean_trial = mean(trial_hists, 2);            
            norm_mean_trial = mean_trial ./ norm(mean_trial);
            temp_cor_vals = trial_hists' * norm_mean_trial;
            temp_mean_cor_val = mean(temp_cor_vals);
            mean_corr(rgc, spat_p, temp_p) = temp_mean_cor_val;
            


            % cor_matrix = trial_hists' * trial_hists;
            % imagesc(cor_matrix)
            % diag_mask = ~eye(size(cor_matrix))
            % pause
            % size(diag_mask)
            % pause
            % median_corr(rgc, spat_p, temp_p) = median(cor_matrix(diag_mask))
            % pause

            % compute OSI and DSI for each spatial and temporal period
            OSIs(rgc, spat_p, temp_p) = (max_val - min_val) ./ (max_val + max_val);
            DSIs(rgc, spat_p, temp_p) = (ds_max - null_val) ./ (ds_max + null_val);
        end
    end

    % decide whether to plot tuning curves and rasters?
    if p.Results.plot_tuning    % is plotting on?
        if p.Results.plot_tuning_os_only    % only plot OS cells?
            % if true, determine if any cross significance thresholds
            insig_DSI = all(DSI(:) < p.Results.DSI_thresh);
            sig_OSIs = any(OSI(:) > p.Results.OSI_thresh);
            sig_cors = any(mean_cor(:) > p.Results.response_cor_thresh);
    
            if sig_OSIs && insig_DSI && sig_cors
                fig_counter = p.Results.start_fig_num;
                for spat_p = 1:num_sps
                    for temp_p = 1:num_tps
                        temp_sp = spatial_periods(spat_p);
                        temp_tp = temp_periods(temp_p);
                        fig_title = ['cell id-', num2str(datarun.cell_ids(cell_indices(rgc))), ' SP-', num2str(temp_sp), ' TP-', num2str(temp_tp)];
                        [temp_tuning, temp_spike_nums, ~] = get_direction_tuning(temp_spike_times, datarun.stimulus,...
                                                                'SP', spatial_periods(spat_p), 'TP', temp_periods(temp_p));
                        plot_direction_tuning(temp_tuning, temp_spike_nums, datarun.stimulus, 'fig_num', fig_counter, 'fig_title', fig_title)  
                        drawnow
                        if p.Results.export_plot
                            exportgraphics(gcf, [fig_title, '.pdf'], 'ContentType', 'image');
                        end
                        fig_counter = fig_counter + 1;
                    end
                end
            end
        end
    end

    % print out values for the current cells
    if p.Results.verbose
        fprintf('OSIs are %g \n', OSI)
        fprintf('DSIs are %g \n', DSI)
        fprintf('cors are %g \n', mean_cor)
    end

    if p.Results.pause
        pause
    end
end
