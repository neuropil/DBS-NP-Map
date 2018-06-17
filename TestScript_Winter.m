

% TEST script for executing the functions to generate DBS lead alignment
% with the firing rates extracted from microelectrode recordings

close all

ele_nii = 'SamplePostOP-Rlead.nii.gz';
mr_nii = 'SamplePostOpBrain.nii';
neuroDAT = 'SampleFiringRate.csv';
coronalSlice = 260;
sagittalSlice = [];
axialSlice = [];
dataColor = 0;
dataSize = 0;
border = 1;

DBS_Align_SpikeParam_05(ele_nii, mr_nii ,...
    coronalSlice, sagittalSlice , axialSlice , dataColor, dataSize , border , neuroDAT)

%%

dataColor = 1;
border = 0;

DBS_Align_SpikeParam_05(ele_nii, mr_nii ,...
    coronalSlice, sagittalSlice , axialSlice , dataColor, dataSize , border , neuroDAT)


%%

dataColor = 1;
border = 1;
dataSize = 1;

DBS_Align_SpikeParam_05(ele_nii, mr_nii ,...
    coronalSlice, sagittalSlice , axialSlice , dataColor, dataSize , border , neuroDAT)