function [OSI,DSI] = get_selectivIDX(temp_spike_nums, num_dirs); 
% compute Orientation-Selective Index & Direction-Selective Index
[ds_max, ds_ind] = max(temp_spike_nums);
condensed_tuning = mean(reshape(temp_spike_nums, [num_dirs/2, 2]), 2);
[max_val, max_ind] = max(condensed_tuning);
            
%handles case of 8 directions
if num_dirs == 8
    
    % for OS
    orth_ind = mod(max_ind+2,4);
    if orth_ind == 0
        orth_ind = 4;
    end
    orth_val = condensed_tuning(orth_ind);
                
    % for DS
    op_ind = mod(ds_ind+4, 8);
    if op_ind == 0
          op_ind = 1;
    end
    null_val = temp_spike_nums(op_ind);
                
% handles case of 12 directions    
elseif num_dirs == 12
                
    % for OS
    orth_ind = mod(max_ind+3,6);
    if orth_ind == 0
        orth_ind = 6;
    end
    orth_val = condensed_tuning(orth_ind);
                
    % for DS
    op_ind = mod(ds_ind+6, 12);
    if op_ind == 0
       op_ind = 1;
    end
    null_val = temp_spike_nums(op_ind);
                
end
            
% compute OSI and DSI for each spatial and temporal period
 OSI = (max_val - orth_val) ./ (max_val + orth_val);
 DSI = (ds_max - null_val) ./ (ds_max + null_val);
