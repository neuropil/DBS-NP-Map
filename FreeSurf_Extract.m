function [ blobPoints , blobBounds ] = FreeSurf_Extract( FS_blobDims )
% FreeSurf_Extract 
% 
% Purpose:
%   Identify faces and vertices for points that represent the 3D dimensions
%   of the Freesurfer region-of-interest
%
% Inputs (required):
%
%   FS_blobDims = Output from 'Extract3DObjectFS'
% 
% Outputs 
%   blobBounds = Output from Matlab 'boundary'
%   blobPoints = matrix with X,Y,Z coordinates for z-levels with binary mask elements
%
% Example:
% 
% *Using NIFTITools to read .nii
% >> 

blobPoly = FS_blobDims.blobDims;

[dbsINDS] = cellfun(@(x) ~isempty(x), blobPoly);

% Create a circle at top Centroid and Bottom Centroid and connect with
% Boundary

blobBounds = boundary(cell2mat(blobPoly(dbsINDS)),0.5); 
blobPoints = cell2mat(blobPoly(dbsINDS)); 


end


