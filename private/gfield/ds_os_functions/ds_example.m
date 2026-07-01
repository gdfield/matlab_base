addpath '/Users/amrudzite/Documents/MATLAB/ds_analysis/dscellanalysis/dscellanalysis';

DataBlock.DsPath='/Users/amrudzite/Desktop/SC_grant/R21-10/LE/2021-07-13-1/data001/data001';
DataBlock.BwPath='/Users/amrudzite/Desktop/SC_grant/R21-10/LE/2021-07-13-1/data002/data002';
Params=[];
DB=1;
[dsall_id] = DsCellFinder2(DataBlock,DB,Params);


