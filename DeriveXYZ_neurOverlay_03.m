function [Xc, Yc, Zc,...
    idC, Xa, Ya,...
    Za, idA, Xp,...
    Yp, Zp, idP,...
    Xm, Ym, Zm,...
    idM, Xl, Yl,...
    Zl, idL, featureOut] =...
    DeriveXYZ_neurOverlay_03( neurData , output_args , sliceTh, lastRecDep)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

rawXYZvals = output_args.centroidsInt;

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

