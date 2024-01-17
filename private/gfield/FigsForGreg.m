% FIGURES
load('/Volumes/lab/temp/cx57_psam_analysis/KoSpatialTuningAnalyzerConcatPlusMappingDB12')

lineWidthRange = 3 ;
figure % psths of each cell in control, +psem, wash
cells=1 ;
while cells <=l_cells ; % for each cell
    clf
    for dset=1:3 ;
        for sp=1:l_sp ; % for each spatial period 
            for cntrst=1:l_cntrst ; % for each contrast
                subplot(l_sp,l_cntrst,l_cntrst*(sp-1)+cntrst) ;
                plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cells,:),'color',Color_list{dset},'linewidth',lineWidthRange*psth_deltaWash{cntrst}{sp}(cells)+.1)
                axis tight
                hold on
            end
        end
    end
    title(num2str(cells))
    nxtval = input('next (cell num, 0=back, return=forward)') ;
    if isempty(nxtval) ;
        cells=cells+1 ;
    elseif nxtval == 0 ;
        cells=cells-1 ;
    elseif nxtval>0 ;
        cells=nxtval ;
    end
end


% psth [control, +PSEM, wash] organized by unique cell types for a single spatial freq and contrast
sp = 7 ; % spatial frequency default
cntrst = 2 ; % contrast default
for uc = 1:length(UniqueCellTypes) ; % for each cell type
    figure
    set(gcf,'name',UniqueCellTypes{uc})
    lc = length(cell_i{uc}) ; % number of cell of that type
    for cells=1:lc ; % for each cell of this type
        subplot(ceil(lc/3),3,cells) 
        for dset=1:3 ; 
            plot(PsthTime,psth{dset}{1}{cntrst}{sp}(cell_i{uc}(cells),:),'color',Color_list{dset})
            hold on
        end
        title(num2str(cell_i{uc}(cells)))
    end
end
    
figure % spatial tuning plots normalized with peak comparisons
for cntrst=1:l_cntrst ; % for each contrast
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst)
        plot(SpatialPeriod, psth_var_norm_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4) 
        hold on
    end
    for dset=1:3 ;
        subplot(2,l_cntrst,cntrst+l_cntrst)
        plot(SpatialPeriod, psth_var_mean{dset}{1}{cntrst},Color_list{dset},'linewidth',4) 
        hold on
    end
end
