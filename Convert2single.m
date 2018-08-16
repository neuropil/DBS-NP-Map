function [mriOUT] = Convert2single(MRI)
% Convert2single 
% 
% Purpose:
%   Converts the data type of MRI to 'single'
%
% Inputs (required):
%   MRI = 3D matrix represented the MRI data
% 
% Outputs 
%   mriOUT = a 3D matrix of similar dimensions to input matrix, however
%   values are stored as 'single'
%
% Example:
% *Using NIFTITools to read .nii
% >> mri_nii = load_nii('MRI.nii')
% >> mriSingle = Convert2single(mri_nii.img)
%
%
% Last edit 8/14/2018


mriOUT = zeros(size(MRI.img),'single');

for mi = 1:size(mriOUT,3)
    
    mriIm = single(MRI.img(:,:,mi));
    mriIm = mriIm/(max(max(mriIm)));
    mriOUT(:,:,mi) = mriIm;
    
end



end
