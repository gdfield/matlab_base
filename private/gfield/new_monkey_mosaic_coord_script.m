%% Analyze mosaics relationships (using Greg's rat data recorded on 512 array) 
datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/Chichilnisky-lab/';
savedir = '~/Desktop/';



n=1; 
ffiles{n} = '2017-11-29-0/data006-gdf/data006'; 

if exist([savedir,'datastr.mat'],'file')
    load([savedir,'datastr.mat']);
end

%% Load data into structure, get sta centers and sta fits 

opt_wn = struct('load_params',1,'load_neurons',1, 'load_ei',1, 'load_sta','all', 'verbose',1);
htype = {'ON','OFF'}; 
vision_type_ON = {'On Parasol','On Midget'};
vision_type_OFF = {'Off Parasol','Off Midget'};
real_type = {'Parasol','Midget'};
real_typename = {'Parasol','Midget'};
if (~exist('datarun','var') || ~exist('data','var') || ~isstruct(datarun{1}) || ~isstruct(data{1}))
    
    [datarun,data] = deal(cell(1,length(ffiles))); 
    for n=1:length(ffiles)
        datarun{n} = load_data([datapath,ffiles{n}], opt_wn);
        types = cell2mat(datarun{n}.cell_types);
        datarun{n}.cell_types = {types.name}; 
        marks_params.thresh = 4.5;
        datarun{n} = get_sta_summaries(datarun{n}, 'all', 'marks_params', marks_params);
        datarun{n} = get_significant_stixels(datarun{n},datarun{n}.cell_ids);
        datarun{n} = calc_rf_snrs(datarun{n}); 
        datarun{n} = rmfield(datarun{n},'globals');
        datarun{n}.stas = rmfield(datarun{n}.stas,'java_sta');

        % Get cell types 
        types = cell2mat(datarun{n}.vision.cell_types); 
        for ii=1:length(vision_type_ON)
            [~,indx,~] = intersect({types.name},vision_type_ON{ii});
            cids = double(datarun{n}.vision.cell_types{indx}.cell_ids); 
            ncells = length(cids); get_cell_indices(datarun{n},cids);  
            cinds = get_cell_indices(datarun{n}, cids);  
            for jj=1:length(cinds)
                tcind = cinds(jj); 
                fnames = fieldnames(datarun{n}.stas.fits{tcind});
                for fn=1:length(fnames)
                    data{n}.(htype{1}).(real_typename{ii}).rf_fits{jj}.(fnames{fn}) = datarun{n}.stas.fits{tcind}.(fnames{fn}); 
                end
                data{n}.(htype{1}).(real_typename{ii}).rf_coms{jj} = datarun{n}.stas.rf_coms{tcind}; 
            end
            data{n}.(htype{1}).(real_typename{ii}).rf_max_snr = datarun{n}.stas.maxsigs(cinds); 
            data{n}.(htype{1}).(real_typename{ii}).rf_med_snr = datarun{n}.stas.medsigs(cinds); 
            data{n}.(htype{1}).(real_typename{ii}).cids = cids; 
            data{n}.(htype{1}).(real_typename{ii}).ncells = ncells; 
            data{n}.(htype{1}).(real_typename{ii}).cinds = cinds; 
        end

        types = cell2mat(datarun{n}.vision.cell_types); 
        for ii=1:length(vision_type_OFF)
            [~,indx,~] = intersect({types.name},vision_type_OFF{ii});
            cids = double(datarun{n}.vision.cell_types{indx}.cell_ids); 
            ncells = length(cids); get_cell_indices(datarun{n},cids);  
            cinds = get_cell_indices(datarun{n}, cids);  
            for jj=1:length(cinds)
                tcind = cinds(jj); 
                fnames = fieldnames(datarun{n}.stas.fits{tcind});
                for fn=1:length(fnames)
                    data{n}.(htype{2}).(real_typename{ii}).rf_fits{jj}.(fnames{fn}) = datarun{n}.stas.fits{tcind}.(fnames{fn});  
                end
                data{n}.(htype{2}).(real_typename{ii}).rf_coms{jj} = datarun{n}.stas.rf_coms{tcind}; 
            end
            data{n}.(htype{2}).(real_typename{ii}).rf_max_snr = datarun{n}.stas.maxsigs(cinds); 
            data{n}.(htype{2}).(real_typename{ii}).rf_med_snr = datarun{n}.stas.medsigs(cinds); 
            data{n}.(htype{2}).(real_typename{ii}).cids = cids; 
            data{n}.(htype{2}).(real_typename{ii}).ncells = ncells; 
            data{n}.(htype{2}).(real_typename{ii}).cinds = cinds; 
        end

        data{n}.stimulus.field_width = datarun{n}.stimulus.field_width;
        data{n}.stimulus.field_height = datarun{n}.stimulus.field_height;
        data{n}.stimulus.stixel_width = datarun{n}.stimulus.stixel_width;
        data{n}.stimulus.stixel_height = datarun{n}.stimulus.stixel_height;
    end
end

clr = [1 0 0; 0 0 1]; 
%save('/Users/Suva/Documents/Projects/RetinalMosaics/datastr.mat','datarun','data','clr','-v7.3'); 



%% Plot receptive field fits and COMs from all data sets 

ctype1 = 1; % ON=1, OFF=2
ctype2 = 2; % ON=1, OFF=2
ctype_sub1 = 1; % 1=parasol, 2=midget, 
ctype_sub2 = 1; % 1=parasol, 2=midget, 

hf= figure(); clf(hf); set(hf,'position',[318   113   826   638]); 
suptitle(sprintf('%s-%s',[htype{ctype1},' ',real_type{ctype_sub1}],[htype{ctype2},' ',real_type{ctype_sub2}])); 
nc = ceil(sqrt(length(ffiles)));
nr = ceil(length(ffiles)/nc); 
for nm=1:length(ffiles)
    subplot(nr,nc,nm);
    
    on_ids = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).cids; 
    off_ids = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).cids; 
    xE = data{nm}.stimulus.field_width; % stixels
    yE = data{nm}.stimulus.field_height; % stixels     
    
    n_on = length(on_ids); 
    n_off = length(off_ids); 
    n_all = n_on+n_off; 

    for cc=1:n_on
        pos = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.mean; 
        sd = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.sd;
        angle = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.angle;
        [X,Y] = drawEllipse([pos sd angle]);
        plot(X,Y,'-','color',clr(1,:)); hold on;
        plot(pos(1),pos(2),'.','color',clr(1,:)); 
    end
    for cc=1:n_off
        pos = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc}.mean;
        sd = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc}.sd;
        angle = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc}.angle;
        [X,Y] = drawEllipse([pos sd angle]);
        plot(X,Y,'-','color',clr(2,:)); hold on;
        plot(pos(1),pos(2),'.','color',clr(2,:)); 
    end
    
    axis equal;
    set(gca,'xlim',[0 xE],'ylim',[0 yE],'fontsize',12);
    xlabel('stixels');
    ylabel('stixels');
    if nm==length(ffiles)
        [hl,ho] = legend([htype{ctype1},real_type{ctype_sub1}], [htype{ctype2},real_type{ctype_sub2}],'location','northeast');
        txt = findobj(ho,'type','text');
        set(txt(1),'color',clr(1,:),'fontsize',10);
        set(txt(2),'color',clr(2,:),'fontsize',10);
        ln = findobj(ho,'type','line'); 
        for ll=1:length(ln) 
            set(ln,'visible','off');
        end
        set(hl,'box','off');
    end
    title(sprintf('data %d',nm));
end


%% Plot receptive field fits and COMs from 1 data set 
ctype1 = 1; % ON=1, OFF=2
ctype2 = 2; % ON=1, OFF=2
ctype_sub1 = 1; % 1=brisk sustained, 2=brisk transient, 3=transient
ctype_sub2 = 1; % 1=brisk sustained, 2=brisk transient, 3=transient

nm = 1; 

on_ids = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).cids; 
off_ids = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).cids; 
xE = data{nm}.stimulus.field_width; % stixels
yE = data{nm}.stimulus.field_height; % stixels     

n_on = length(on_ids); 
n_off = length(off_ids); 
n_all = n_on+n_off; 

pos = zeros(n_all,2); 
[X,Y] = deal(cell(1,n_all)); 
for cc=1:n_on
    pos(cc,:) = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.mean; 
    sd = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.sd;
    angle = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.angle;
    [X{cc},Y{cc}] = drawEllipse([pos(cc,:) sd angle]);
end
for cc=n_on+1:n_all
    pos(cc,:) = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc-n_on}.mean;
    sd = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc-n_on}.sd;
    angle = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc-n_on}.angle;
    [X{cc},Y{cc}] = drawEllipse([pos(cc,:) sd angle]);
end

% Median cell-to-cell distance 
temp_on = pos(1:n_on,:);
temp_off = pos(1+n_on:n_all,:);
c2c_dist_on = zeros(n_on,1); 
for ii=1:n_on
    a = temp_on(ii,:); 
    b = temp_on(setxor(1:n_on,ii),:);
    DT = delaunayTriangulation(b);
    [~,c2c_dist_on(ii)] = nearestNeighbor(DT,a); 
end
c2c_dist_off = zeros(n_off,1);         
for ii=1:n_off
    a = temp_off(ii,:); 
    b = temp_off(setxor(1:n_off,ii),:);
    DT = delaunayTriangulation(b);
    [~,c2c_dist_off(ii)] = nearestNeighbor(DT,a); 
end
med_c2c_dist_on = median(c2c_dist_on);
med_c2c_dist_off = median(c2c_dist_off);
d_true_stx = min([med_c2c_dist_on med_c2c_dist_off]); 
% Illustrating ROI (in stixel units)
stxBufferX = 1.0*d_true_stx; 
stxBufferY = 1.0*d_true_stx; 
Left = min(pos(:,1));
Right = max(pos(:,1)); 
Bottom = min(pos(:,2));
Top = max(pos(:,2));
ROIx = [Left+stxBufferX Left+stxBufferX Right-stxBufferX Right-stxBufferX Left+stxBufferX]';
ROIy = [Bottom+stxBufferY Top-stxBufferY Top-stxBufferY Bottom+stxBufferY Bottom+stxBufferY]';
 

figure();
subplot(1,2,1);
for cc=1:n_on
    plot(X{cc},Y{cc},'-','color',clr(1,:)); hold on;
    plot(pos(cc,1),pos(cc,2),'.','color',clr(1,:)); 
end
plot(ROIx, ROIy, '--k','linewidth',2); 
axis equal;
set(gca,'xlim',[0 xE],'ylim',[0 yE],'fontsize',12);
xlabel('stixels');
ylabel('stixels');
title(sprintf('%s',[htype{ctype1},' ',real_type{ctype_sub1}]));
subplot(1,2,2);
for cc=n_on+1:n_all
    plot(X{cc},Y{cc},'-','color',clr(2,:)); hold on;
    plot(pos(cc,1),pos(cc,2),'.','color',clr(2,:)); 
end
plot(ROIx, ROIy, '--k','linewidth',2); 
axis equal;
set(gca,'xlim',[0 xE],'ylim',[0 yE],'fontsize',12);
xlabel('stixels');
ylabel('stixels');
title(sprintf('%s',[htype{ctype2},' ',real_type{ctype_sub2}]));
  


%% Calculate pairwise interaction energy 

micron2pix = 4; % conversion factor 
include_STA_struct = false; 
ctype1 = 1; % ON=1, OFF=2
ctype2 = 2; % ON=1, OFF=2
ctype_sub1 = 1; % 1=brisk sustained, 2=brisk transient, 3=transient
ctype_sub2 = 1; % 1=brisk sustained, 2=brisk transient, 3=transient
move = 'OFF'; % which mosaic to move: ON or OFF 
ROI_square = false; % whether to use a square ROI instead of a rectangular one 
buf = 1.5; % ROI edge buffer in units of inter-cell distance 
bin = 0.04; % step size to move mosaic (unit of inter-cell-distance)
clear accum steps_mat;
% figure; 



for nm=1:length(data) 
    
    % Load cell ids
    on_ids = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).cids; 
    on_inds = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).cinds; 
    off_ids = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).cids; 
    off_inds = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).cinds;

    % Store cell positions, and spatial RF parameters 
    n_on = length(on_ids);
    n_off = length(off_ids);
    n_all = n_on+n_off;
    [pos,rad] = deal(zeros(n_all,2)); 
    ang = zeros(n_all,1); 
    for cc=1:n_on
        pos(cc,:) = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.mean; 
        rad(cc,:) = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.sd; 
        ang(cc,:) = data{nm}.(htype{ctype1}).(real_typename{ctype_sub1}).rf_fits{cc}.angle; 
    end
    for cc=n_on+1:n_all
        pos(cc,:) = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc-n_on}.mean;
        rad(cc,:) = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc-n_on}.sd;
        ang(cc,:) = data{nm}.(htype{ctype2}).(real_typename{ctype_sub2}).rf_fits{cc-n_on}.angle;
    end
    pos_on = pos(1:n_on,:); 
    pos_off = pos(1+n_on:n_all,:); 

    % Median cell-to-cell distance 
    temp_on = pos_on;
    temp_off = pos_off;
    c2c_dist_on = zeros(n_on,1); 
    for ii=1:n_on
        a = temp_on(ii,:); 
        b = temp_on(setxor(1:n_on,ii),:);
        DT = delaunayTriangulation(b);
        [~,c2c_dist_on(ii)] = nearestNeighbor(DT,a); 
    end
    c2c_dist_off = zeros(n_off,1);         
    for ii=1:n_off
        a = temp_off(ii,:); 
        b = temp_off(setxor(1:n_off,ii),:);
        DT = delaunayTriangulation(b);
        [~,c2c_dist_off(ii)] = nearestNeighbor(DT,a); 
    end

    med_c2c_dist_on = median(c2c_dist_on);
    med_c2c_dist_off = median(c2c_dist_off);

    % Define lattice spacing in units of NN distance 
    %     Note: 
    %           The average energy map is obtained by averaging the energy 
    %           maps from different data sets (or preps). However in each 
    %           prep, the NN cell-cell distance can be different, which 
    %           means that the absolute mosaic shift distance to cover 1 NN
    %           distance can be different in different data sets. So a 
    %           constant shift applied uniformly to all data sets would 
    %           yield energy maps with X-Y axes out of scale wrt each other
    %           in units of the dimensionless inter-cell distance. To avoid
    %           this, we shift mosaics in each data by ~1.5*NN distance 
    %           corresponding to the data set.  
    
    if strcmp(move,'ON') 
        d_true_stx = med_c2c_dist_off; % in stixels 
    elseif strcmp(move,'OFF') 
        d_true_stx = med_c2c_dist_on; % in stixels 
    end
    
    % Parameter values, if eqn. from Jang, 2017 is used
    R = 1.5*d_true_stx; % interaction distance (this can be varied to change the strength of interaction)
    epsil = 0.01; % fixed constant used for setting F=0 at r<=R
    soma_d_micron = 16; % ~ micron (mouse RGCs; Ivanova, 2013)
    s = soma_d_micron/micron2pix/data{nm}.stimulus.stixel_width; 
    sf = d_true_stx; % if you need to scale with respect to the unit lattice distance 

    % Define outer bound of ROI (in stixel units)
    stxBufferX = buf*d_true_stx; 
    stxBufferY = buf*d_true_stx; 
    Left = min([min(pos_on(:,1)) min(pos_off(:,1))]);
    Right = max([max(pos_on(:,1)) max(pos_off(:,1))]);
    Bottom = min([min(pos_on(:,2)) min(pos_off(:,2))]);
    Top = max([max(pos_on(:,2)) max(pos_off(:,2))]);
    ROIx = [Left+stxBufferX Left+stxBufferX Right-stxBufferX Right-stxBufferX Left+stxBufferX]';
    ROIy = [Bottom+stxBufferY Top-stxBufferY Top-stxBufferY Bottom+stxBufferY Bottom+stxBufferY]';
    W = max(ROIx)-min(ROIx); 
    H = max(ROIy)-min(ROIy); 
    % Making the ROI a perfect square
    if ROI_square
        ROIx([1 2 5]) = ROIx([1 2 5])+(W-H)/2; 
        ROIx([3 4]) = ROIx([3 4])-(W-H)/2; 
        W = max(ROIx)-min(ROIx); 
        H = max(ROIy)-min(ROIy); 
    end 


    % Set up lattice shift parameters, ROI, and grid matrix 
    gamma = 1.2; % max mosaic shift (units of NN distance)
    x_range = sf.*[-gamma:bin:gamma]; % shift steps (units of stixel)
    y_range = x_range;
    num_steps = length(x_range);
    
    % Closest distance : 0.20 times the lattice distance or inter-cell dist
    min_sep = 0.20*d_true_stx; 

    % Calculate energy 
    [cellcount,ON_pos_temp,OFF_pos_temp] = deal([]); 
    d_ = d_true_stx; 
    energy = zeros(num_steps,num_steps);  
    ON_pos_initial = pos(1:n_on,:); 
    OFF_pos_initial = pos(1+n_on:n_all,:); 

    %figure; 
    for x_step = 1:num_steps
        for y_step = 1:num_steps 

            if strcmp(move,'ON')  
                % freeze on cell mosaic
                ON_pos_temp = [ON_pos_initial(:,1)+x_range(x_step) ON_pos_initial(:,2)+y_range(y_step)];
                % shift off cell mosaic 
                OFF_pos_temp = OFF_pos_initial;
            elseif strcmp(move,'OFF') 
                % freeze on cell mosaic
                ON_pos_temp = ON_pos_initial;
                % shift off cell mosaic 
                OFF_pos_temp = [OFF_pos_initial(:,1)+x_range(x_step) OFF_pos_initial(:,2)+y_range(y_step)];
            end

            % Get the cells inside ROI 
            in_on_updated = inpolygon(ON_pos_temp(:,1), ON_pos_temp(:,2), ROIx, ROIy); 
            in_off_updated = inpolygon(OFF_pos_temp(:,1), OFF_pos_temp(:,2), ROIx, ROIy); 

            % Get heterotypic pair-wise distance 
            pairwise_dists = ipdm(ON_pos_temp(in_on_updated,:), OFF_pos_temp(in_off_updated,:)); % row: ON, col: OFF             

            % Set energy for particles closer than min_separation to
            % constant
            pairwise_dists(pairwise_dists<min_sep) = min_sep; 
            
            % Choose only pairwise distance > soma diameter 
            %pairwise_dists = pairwise_dists((pairwise_dists>s)); 

            % Calculate energy 
            %energy_temp = epsil* ( (((R-s)./d_true_stx).^2)./(pairwise_dists-s) +  pairwise_dists );
            %energy_temp = epsil* ( (((R-s)./d_true_stx).^2)./(pairwise_dists) +  pairwise_dists );
            %energy_temp = (1./(pairwise_dists)) + pairwise_dists;
            energy_temp = (1./pairwise_dists);
            energy(x_step,y_step) = sum(energy_temp(:))./length(energy_temp(:)); 

%             cla(gca); 
%             plot(ON_pos_temp(:,1),ON_pos_temp(:,2),'o','color',clr(1,:)); hold on; 
%             plot(OFF_pos_temp(:,1),OFF_pos_temp(:,2),'o','color',clr(2,:)); hold on; 
%             plot(ON_pos_temp(in_on_updated,1),ON_pos_temp(in_on_updated,2),'^m');
%             plot(OFF_pos_temp(in_off_updated,1),OFF_pos_temp(in_off_updated,2),'^m');
%             plot(ROIx,ROIy,'--k'); axis equal; 
%             set(gca,'xlim',[min(bound(:,1))-s max(bound(:,1))+s],'ylim',[min(bound(:,2))-s max(bound(:,2))+s]);
%             pause(0.01);

        end
    end
    fprintf('Data %d out of %d \n',nm,length(ffiles)); 


    % Smooth energy surface 
    energy_smooth = imgaussfilt(energy, 0.25*d_true_stx);  

    % Gather raw energy maps and smoothed energy maps
    accum.energy_map(:,:,nm) = energy; 
    accum.energy_smooth_map(:,:,nm) = energy_smooth;

    % Gather COMs and ROIs
    accum.ON_pos_initial{nm} = ON_pos_initial; 
    accum.OFF_pos_initial{nm} = OFF_pos_initial; 
    accum.Left{nm} = Left; 
    accum.Right{nm} = Right;
    accum.Bottom{nm} = Bottom;
    accum.Top{nm} = Top;
    accum.ROIx{nm} = ROIx; 
    accum.ROIy{nm} = ROIy; 
    accum.d_true_stx{nm} = d_true_stx; 
    accum.x_range{nm} = x_range; 
end
fprintf(sprintf('Mosais: %s, %s \n',[htype{ctype1},' ',real_typename{ctype_sub1}],[htype{ctype2},' ',real_typename{ctype_sub2}]));


% Gather matrix of X,Y,rho,theta corresponding to each shift (units of unitless NN distance)
sc_x_range = -gamma:bin:gamma;
sc_y_range = sc_x_range;
global_sc_steps_mat = zeros(length(sc_x_range),length(sc_x_range),4); 
for x_step = 1:num_steps
    for y_step = 1:num_steps 
        global_sc_steps_mat(x_step,y_step,1) = sc_x_range(x_step); % Z dimension: X,Y,r,theta
        global_sc_steps_mat(x_step,y_step,2) = sc_x_range(y_step); 
        [th,rd] = cart2pol(sc_x_range(x_step),sc_x_range(y_step)); 
        global_sc_steps_mat(x_step,y_step,3) = th;
        global_sc_steps_mat(x_step,y_step,4) = rd;
    end
end



% Figure of energy maps
hf= figure(12); clf(hf); %set(hf,'position',[-1536 385 1413 310]); 
suptitle(sprintf('%s-%s',[htype{ctype1},' ',real_type{ctype_sub1}],[htype{ctype2},' ',real_type{ctype_sub2}])); 
nc = ceil(sqrt(length(ffiles)));
nr = ceil(length(ffiles)/nc); 
for nm=1:length(ffiles)
    subplot(nr,nc,nm);
    localenergy = accum.energy_map(:,:,nm); 
    if robust_std(localenergy(:))~=0        
        localenergy = (localenergy-mean(localenergy(:)))./robust_std(localenergy(:)); 
    end
    localenergy = imgaussfilt(localenergy, 0.5*d_true_stx); % gaussian filter with sigma = 1/2*inter-cell dist
    h = pcolor(global_sc_steps_mat(:,:,1),global_sc_steps_mat(:,:,2), localenergy); hold on; 
    plot([sc_x_range(1) sc_x_range(end)],[0 0],'-w',[0 0],[sc_y_range(1) sc_y_range(end)],'-w');
    colorbar; 
    set(h,'EdgeColor','None'); axis square; 
    set(gca,'fontsize',12);
end
ylabel('Mosaic shift (unit of d)');



%% Zscore energy. Plot energy vs radial shift 

% Average Z scored energy maps
energy_map = zeros(size(accum.energy_smooth_map,1), size(accum.energy_smooth_map,2)); 
for nm=1:size(accum.energy_smooth_map,3) 
    temp = accum.energy_smooth_map(:,:,nm); 
    if robust_std(temp(:))~=0
        temp = (temp-mean(temp(:)))./robust_std(temp(:)); 
    end
    energy_map = energy_map + temp; 
end
mean_energy_map = energy_map./size(accum.energy_smooth_map,3); 
mean_energy_map_zsc = (mean_energy_map-mean(mean_energy_map(:)))./robust_std(mean_energy_map(:));


% set axes units for image 
x_un = sc_x_range;
y_un = fliplr(sc_y_range);
tickmarks = [-1 0 1];
[~,lcx] = intersect(x_un, tickmarks); 
vlx = sc_x_range(lcx); 
[~,lcy] = intersect(y_un, tickmarks); 
if ~isequal(lcy,sort(lcy,'ascend'))
    [lcy,invrt] = sort(lcy,'ascend'); 
    vly = tickmarks(invrt); 
end
    
    

% Plot averaged energy maps 
FS = 20; 
figure(13); 
imagesc(mean_energy_map_zsc); hold on; 
plot([1 length(sc_x_range)],repmat((length(sc_x_range)+1)/2,1,2),'-w','linewidth',2);
plot(repmat((length(sc_y_range)+1)/2,1,2),[1 length(sc_y_range)],'-w','linewidth',2);
caxis([-1 1].*max(mean_energy_map_zsc(:))); 
hc = colorbar; ylabel(hc,'Z-score'); axis square;
set(hc,'fontsize',FS);
xlabel('Mosaic shift (unit of inter-cell dist)');
title(sprintf('Averaged over %d data sets. \n %s, %s mosaics',length(ffiles),[htype{ctype1},' ',real_type{ctype_sub1}],[htype{ctype2},' ',real_type{ctype_sub2}]));
set(gca,'fontsize',FS);
set(gca,'xtick',lcx,'xticklabel',vlx,'ytick',lcy,'yticklabel',vly); 



% Get radial energy dependence  
r_bin_edge = linspace(0,max(sc_x_range),21);
r_bin_cent = (r_bin_edge(1:end-1) + r_bin_edge(2:end))./2;
r_max = r_bin_edge(end);
[lin_energy, lin_energy_std] = deal(zeros(1,length(r_bin_cent))); 
for rb=1:length(r_bin_edge)-1
    linind = find(global_sc_steps_mat(:,:,4)>r_bin_edge(rb) & global_sc_steps_mat(:,:,4)<=r_bin_edge(rb+1));
    lin_energy(rb) = median(mean_energy_map_zsc(linind));
    lin_energy_std(rb) = std(mean_energy_map_zsc(linind)); 
end
rads = global_sc_steps_mat(:,:,4); 
[sorted_rads, sorted_inds] = sort(rads(:),'ascend');
sorted_energies = mean_energy_map_zsc(sorted_inds); 


% Plot radial energy 
figure(14); cla; 
errorbar(r_bin_cent, lin_energy,lin_energy_std./2,'o-b','linewidth',2); hold on; 
plot([0 r_bin_cent(end)],[0 0],'--k','linewidth',2); 
xlabel('Mosaic shift (unit of inter-cell dist)');
ylabel('Z-score');
title('Energy profile');
set(gca,'fontsize',20);


% Get distribution of radial energy and plot 
figure(15); cla; 
plot([0 ceil(sc_x_range(end)*sqrt(2))],[0 0],'--k','linewidth',2); hold on; 
plot(sorted_rads, sorted_energies,'.b'); hold on; 
xlabel('Radial dist (unit of inter-cell dist)');
ylabel('Z-scored energy');
set(gca,'fontsize',20);



% Plot radial scatter and mean together 
clx = [0    0.4470    0.7410]; % from Matlab's parula color map
figure(); cla; 
% scatter(sorted_rads, sorted_energies,15,...
%     'MarkerEdgeColor','none','MarkerFaceColor','k','MarkerFaceAlpha',0.3); hold on; 
errorbar(r_bin_cent, lin_energy,lin_energy_std./2,'o-b','linewidth',2); hold on; 
plot([0 r_bin_cent(end)],[0 0],'--k','linewidth',2);
ylabel('Z-scored energy');
xlabel('Radial distance'); 
set(gca,'xlim',[0 1.2],'fontsize',25);

