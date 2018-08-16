function [] = Plot3D_EleBoundary( eleStruct )
% PLOT3D_ELECBOUNDARY
% 
% Purpose:
%   Converts the data type of MRI to 'single'
%
% Inputs (required):
%   eleStruct = output struct from 'ExtractDBSPolygon'
% 
% Outputs 
%    A plot is generated as an output
%
% Example:
% *Using NIFTITools to read .nii
% >> ele_nii = load_nii('RstnElectrodeTrace.nii');
% >> eleDiamMM = 1.27;
% >> numPtsCircle = 80;
% >> extrPolyOutput = ExtractDBSPolygon(ele_nii, eleDiamMM, numPtsCircle);
% >> Plot3D_EleBoundary( extrPolyOutput )
%
%
% Last edit 8/14/2018

zINDS = cellfun(@(x) ~isempty(x), eleStruct.meanCirSMint(:,1));
circlesAll = eleStruct.meanCirSMint(zINDS);
circlesAll = cellfun(@(x) x(:,1:3), circlesAll, 'UniformOutput', false);

circleContainer = [];

for i = 1:length(circlesAll)
   
    tmp = circlesAll{i};
    
    circleContainer = [circleContainer ; tmp];

end

boundaryObj = boundary(circleContainer, 0.5);
hold on

trisurf(boundaryObj,circleContainer(:,1),circleContainer(:,2),circleContainer(:,3),'Facecolor','black','FaceAlpha',0.3,'Edgecolor','none')

xlim([0 512])
ylim([0 512])

end

