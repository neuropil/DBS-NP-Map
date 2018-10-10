function [Xc, Yc, Zc,...
    idC, Xa, Ya,...
    Za, idA, Xp,...
    Yp, Zp, idP,...
    Xm, Ym, Zm,...
    idM, Xl, Yl,...
    Zl, idL, featureOut] =...
    DeriveXYZ_NEUROverlay(neurData , extrPolyOutput , sliceTh)
% DERIVEXYZ_NEUROVERLAY
% 
% Purpose:
%   Generates the SIZE and COLOR parameters for the X,Y,Z coordinates for
%   visual display - used in Matlab SCATTER
%
% Inputs (required):
%   'neuroDat' = Matlab Table array with two variables: depth from target,
%              orientation in microelectrode holder (e.g., a = anterior, c = center).
%   'extrPolyOutput' = output struct from 'ExtractDBSPolygon'
%   'sliceTh' = slice thickness in mm
% 
% Outputs 
%   For each possible recording trajectory in the microelectrode holder:
%   anteior, posterior, center, lateral, and medial, the XYZ coordinates
%   and ephys parameter are stored in vector for each output.
%   NOTE: if the trajectory is not present, the output will be an empty
%   vector
%
% Example:
% *Using NIFTITools to read .nii
% >> ele_nii = load_nii('RstnElectrodeTrace.nii');
% >> eleDiamMM = 1.27;
% >> numPtsCircle = 80;
% >> extrPolyOutput = ExtractDBSPolygon(ele_nii, eleDiamMM, numPtsCircle);
% >> neuroDat = readtable('ephysDATA.csv');
% >> solidCol = 1;
% >> sizeMag = 0;
% >> sliceThick = 2; Can be obtained from nii struct NII.hdr.dime.pixdim(4)
% >> Plot3D_EleBoundary( extrPolyOutput );
% >> [Xc, Yc, Zc, idC,...
%     Xa, Ya, Za, idA, Xp,...
%     Yp, Zp, idP, Xm, Ym,...
%     Zm, idM, Xl, Yl, Zl,...
%     idL, featureOut] = DeriveXYZ_NEUROverlay( neuroDat ,...
%     extrPolyOutput , sliceThick);

% Last edit 8/14/2018

lastRecDep = neurData{height(neurData),4} - 100;
rawXYZvals = extrPolyOutput.centroidsInt;

oldX = rawXYZvals(:,1);
oldY = rawXYZvals(:,2);
oldZ = rawXYZvals(:,3);
% oldZs = rawXYZvals(:,4);

xVals = linspace(min(oldX),max(oldX),numel(oldX)*5);
yVals = interp1(oldX,oldY,xVals,'linear','extrap');
zValsMRI = interp1(oldX,oldZ,xVals,'linear','extrap');
% zValsEphys = interp1(oldX,oldZs,xVals,'linear','extrap');

cen = 1;
ant = 1;
post = 1;
med = 1;
lat = 1;

lenData = height(neurData);

Xc = nan(lenData,1);
Yc = nan(lenData,1);
Zc = nan(lenData,1);
idC = nan(lenData,1);

Xa = nan(lenData,1);
Ya = nan(lenData,1);
Za = nan(lenData,1);
idA = nan(lenData,1);

Xp = nan(lenData,1);
Yp = nan(lenData,1);
Zp = nan(lenData,1);
idP = nan(lenData,1);

Xm = nan(lenData,1);
Ym = nan(lenData,1);
Zm = nan(lenData,1);
idM = nan(lenData,1);

Xl = nan(lenData,1);
Yl = nan(lenData,1);
Zl = nan(lenData,1);
idL = nan(lenData,1);

featureOut = nan(lenData,1);

for ii = 1:height(neurData)
    
    tempCell = neurData(ii,:);
    trackID = tempCell.trackAnat{1};
    featureOut(ii) = tempCell.param;
    
    switch trackID
        case 'c'
            
            Zc(cen) = (((tempCell.depth - lastRecDep)/1000)/sliceTh) + min(zValsMRI); %
            zInd = knnsearch(zValsMRI',Zc(cen));
            
            % Find Z index in Zvals
            Xc(cen) = xVals(zInd);
            Yc(cen) = yVals(zInd);
            
            idC(cen) = ii;
            cen = cen + 1;
        case 'a'
            
            Za(ant) = (((tempCell.depth - lastRecDep)/1000)/sliceTh) + min(zValsMRI);
            zInd = knnsearch(zValsMRI',Za(ant));
            
            Xa(ant) = xVals(zInd) + (2 * sliceTh);
            Ya(ant) = yVals(zInd);
            
            idA(ant) = ii;
            ant = ant + 1;
        case 'p'
            
            Zp(post) = (((tempCell.depth - lastRecDep)/1000)/sliceTh) + min(zValsMRI);
            zInd = knnsearch(zValsMRI',Zp(post));
            
            Xp(post) = xVals(zInd) - (2 * sliceTh);
            Yp(post) = yVals(zInd);
            
            idP(post) = ii;
            post = post + 1;
        case 'm'
            
            Zm(med) = (((tempCell.depth - lastRecDep)/1000)/sliceTh) + min(zValsMRI);
            zInd = knnsearch(zValsMRI',Zm(med));
            
            Xm(med) = xVals(zInd);
            Ym(med) = yVals(zInd) - (2 * sliceTh);
            
            idM(med) = ii;
            med = med + 1;
        case 'l'
            
            Zl(lat) = (((tempCell.depth - lastRecDep)/1000)/sliceTh) + min(zValsMRI);
            zInd = knnsearch(zValsMRI',Zl(lat));

            Xl(lat) = xVals(zInd);
            Yl(lat) = yVals(zInd) + (2 * sliceTh);
            
            idL(lat) = ii;
            lat = lat + 1;
    end
end














end

