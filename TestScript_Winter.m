

% TEST script for executing the functions to generate DBS lead alignment
% with the firing rates extracted from microelectrode recordings

close all

ele_nii = 'c260_NATele.nii.gz';
mr_nii = 'c260_brain.nii';
neuroDAT = 'neurodata.csv';
coronalSlice = 260;
sagittalSlice = [];
axialSlice = [];
dataColor = 0;
dataSize = 0;
border = 1;

DisplayEphys2DBS(ele_nii, mr_nii ,...
    coronalSlice, sagittalSlice , axialSlice , dataColor, dataSize , border , neuroDAT)

%%

dataColor = 1;
border = 0;

DisplayEphys2DBS(ele_nii, mr_nii ,...
    coronalSlice, sagittalSlice , axialSlice , dataColor, dataSize , border , neuroDAT)


%%

dataColor = 1;
border = 1;
dataSize = 1;

DisplayEphys2DBS(ele_nii, mr_nii ,...
    coronalSlice, sagittalSlice , axialSlice , dataColor, dataSize , border , neuroDAT)

%%







