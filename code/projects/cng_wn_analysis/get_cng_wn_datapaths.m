function datapath = get_cng_wn_datapaths(cohort, ndf_val)


if ndf_val ~= 0 && ndf_val ~= 2
    error('ndf_val not recognized, must be 0 or 2')
end            


switch cohort
    case 'WT'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-12-16-0/data008/data008';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-05-13-0/data006/data006';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-04-14-0/data005/data005';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data006/data006';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-12-16-0/data004/data004';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-05-13-0/data004/data004';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-04-14-0/data003/data003';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-22-0/data003/data003';
        end
        
    case '1M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-18-0/data005/data005';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-19-0/data005/data005';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-04-21-0/data006/data006';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-11-19-0/data006/data006';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-18-0/data002/data002';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-19-0/data002/data002';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-04-21-0/data003/data003';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-11-19-0/data003/data003';
        end

    case '2M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-26-0/data005/data005';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-15-0/data005/data005';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-30-0/data005/data005';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-03-15-0/data004/data004';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-07-26-0/data002/data002';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-08-15-0/data002/data002';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-30-0/data002/data002';
            %datapath{4} = no 4th ndf 2 data data
        end        
      
    case '3M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-16-0/data005/data005';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-11-07-0/data005/data005';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-09-17-0/data006/data006';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-01-07-0/data005/data005';
            datapath{5} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-03-14-0/data004/data004';
            datapath{6} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-12-0/data004/data004';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-16-0/data002/data002';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-11-07-0/data002/data002';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-09-17-0/data003/data003';
            datapath{4} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-01-07-0/data002/data002';
            datapath{5} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-03-14-0/data001/data001';
            datapath{6} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-04-12-0/data002/data002';
        end        
        
    case '4M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-12-13-0/data005/data005';
            %datapath{2} ='/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-03-16-0/data005/data005'; %% very few ON cells
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-05-10-0/data006/data006';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-07-13-0/data005/data005';
            %datapath{5} = /Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-11-10-0/data006/data006';%poor recording
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-12-13-0/data002/data002';
            %datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-03-16-0/data002/data002';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-05-10-0/data004/data004';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-07-13-0/data003/data003';
            %datapath{5} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-11-10-0/data003/data003';
        end
        
    case '5M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-10-0/data005/data005';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-11-0/data005/data005';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-18-0/data005/data005';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-10-0/data002/data002';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-11-0/data002/data002';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2019-10-18-0/data002/data002';
        end
 
    case '6M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-07-06-0/data004/data004';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-07-06-0/data002/data002';
        end             
        
    case '7M'
        if ndf_val == 0
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-05-19-0/data005/data005';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-01-22-0/data005/data005';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-05-28-0/data005/data005';
        elseif ndf_val == 2
            datapath{1} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-05-19-0/data003/data003';
            datapath{2} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2021-01-22-0/data002/data002';
            datapath{3} = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2020-05-28-0/data003/data003';
        end        
end

