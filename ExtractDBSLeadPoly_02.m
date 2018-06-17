function [ output_args ] = ExtractDBSLeadPoly_02( ele_nii  , solidType , actELEDIM , numPts)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% OLD extractDBSLeadPoly2

% Requires NIFTI tools = URL*************
%  Part of this file is copied and modified from:
%  http://www.mathworks.com/matlabcentral/fileexchange/1878-mri-analyze-tools
%  
%  NIFTI data format can be found on: http://nifti.nimh.nih.gov
%  
%  - Jimmy Shen (jimmy@rotman-baycrest.on.ca)

%% Electrode

% Get mean circle size
% Compute solid cylinder

niiE = load_nii(ele_nii);

sliceTHICK = niiE.original.hdr.dime.pixdim(4);

testIMS = niiE.img;

imInd = false(size(testIMS,3),1);
blobDims = cell(size(testIMS,3),1);  % xyz
circleDims = cell(size(testIMS,3),1);% xyz
centroidS = nan(size(testIMS,3),4); % xyz
diameterS = nan(size(testIMS,3),1);
actEledims = cell(size(testIMS,3),1);% xyz
actEleDia = nan(size(testIMS,3),1);

for iii = 1:size(testIMS,3)
   
    tmpIM = testIMS(:,:,iii);
%     tmpIMflv = flip(tmpIM,1);
%     tmpIMflh = flip(tmpIMflv,2);
    
    tmpIM8 = uint8(tmpIM);
    tmpIM8(tmpIM8 ~= 0) = 255;
    
    if max(max(tmpIM8)) ~= 255
        continue
    else
        imInd(iii) = 1;
        
        tmpCoords = bwboundaries(tmpIM8);
        xyCords = tmpCoords{1};
        lenCords = size(xyCords,1);
        
        blobDims{iii}(1:lenCords,1:2) = xyCords;
        blobDims{iii}(1:lenCords,3) = iii;
        
        regDat = regionprops(imbinarize(tmpIM8),'Centroid');
        
        if size(regDat,1) > 1
            regMean = zeros(size(regDat,1),2);
            for ri = 1:size(regDat,1)
                regMean(ri,1) = regDat(ri).Centroid(1);
                regMean(ri,2) = regDat(ri).Centroid(2);
            end
            
            regMeans = mean(regMean);
            regDat = struct;
            regDat.Centroid(1) = regMeans(1);
            regDat.Centroid(2) = regMeans(2);
        end

        centroidS(iii,1:2) = regDat.Centroid;
        centroidS(iii,3) = iii;
        centroidS(iii,4) = iii*sliceTHICK;
        
        [ eqDia ] = getREgionDims(xyCords(:,2), xyCords(:,1));
        
        [xunit, yunit] = SpecifyCircle_01(regDat.Centroid(1),regDat.Centroid(2),eqDia - 2 , numPts);
        [aExunit, aEyunit] = SpecifyCircle_01(regDat.Centroid(1),regDat.Centroid(2), actELEDIM , numPts);
        
        diameterS(iii,1) = eqDia - 2;
        actEleDia(iii,1) = actELEDIM;
        
        lenEle = size(aExunit,2);  
        
        actEledims{iii}(1:lenEle,1:2) = transpose([aExunit ; aEyunit]);
        actEledims{iii}(1:lenEle,3) = iii;

        lenCir = size(xunit,2);

        circleDims{iii}(1:lenCir,1:2) = transpose([xunit ; yunit]);
        circleDims{iii}(1:lenCir,3) = iii;
        
        
    end

end

%% Compute smooth electrode
orgCenX = centroidS(:,1);
nanIND = ~isnan(orgCenX);
orgCenX = orgCenX(nanIND);
orgCenY = centroidS(:,2);
orgCenY = orgCenY(nanIND);
orgCenZ = centroidS(nanIND,3);
orgCenZs = centroidS(nanIND,4);

smCenY = smooth(orgCenX,orgCenY,0.7,'rloess');
smCenZ = smooth(orgCenX,orgCenZ,0.7,'rloess');
smCenZs = smooth(orgCenX,orgCenZs,0.7,'rloess');

centroidSM = nan(size(centroidS,1),4);
centroidSM(nanIND,1) = orgCenX;
centroidSM(nanIND,2) = smCenY;
centroidSM(nanIND,3) = smCenZ;
centroidSM(nanIND,4) = smCenZs;


%% Compute interpolation



intSMCenX = linspace(min(orgCenX),max(orgCenX),numel(orgCenX)*3);

if sum(ismember(intSMCenX,orgCenX)) ~= 0
    intSMCenX(ismember(intSMCenX,orgCenX)) = intSMCenX(ismember(intSMCenX,orgCenX)) + 0.05;
end

if sum(ismember(intSMCenX,smCenZ)) ~= 0
    intSMCenX(ismember(intSMCenX,smCenZ)) = intSMCenX(ismember(intSMCenX,smCenZ)) + 0.05;
end

if length(smCenZ) ~= length(unique(smCenZ))
    tbl = tabulate(smCenZ);
    toAdj = tbl(tbl(:,2) > 1,1);
    for ai = 1:length(toAdj)
        
        indi = ismember(smCenZ,toAdj(ai));
        
        randAdd = rand(sum(indi),1);
        
        smCenZ(indi) = smCenZ(indi) + randAdd;
    end
    
end

if length(orgCenX) ~= length(unique(orgCenX))
    tbl = tabulate(orgCenX);
    toAdj = tbl(tbl(:,2) > 1,1);
    for ai = 1:length(toAdj)
        
        indi = ismember(orgCenX,toAdj(ai));
        
        randAdd = rand(sum(indi),1);
        
        orgCenX(indi) = orgCenX(indi) + randAdd;
    end
    
end

intSMCenY = interp1(orgCenX,smCenY,intSMCenX,'linear','extrap');
intSMCenZ = interp1(orgCenX,smCenZ,intSMCenX,'linear','extrap');
intSMCenZs = interp1(orgCenX,smCenZs,intSMCenX,'linear','extrap');

centroidsInt = transpose([intSMCenX ; intSMCenY ; intSMCenZ ; intSMCenZs]);

meanAEDia = nanmean(actEleDia);
meanCirSMint = cell(size(centroidsInt,1),1);% xyz
for ccII = 1:size(centroidsInt,1)
    
    [xunitINTsm, yunitINTsm] = SpecifyCircle_01(centroidsInt(ccII,1),centroidsInt(ccII,2),meanAEDia, numPts); % original centroid , actual electrode diameter
    
    % original centroid , original electrode diameter
    lenCirSMint = size(xunitINTsm,2);
    meanCirSMint{ccII}(1:lenCirSMint,1:2) = transpose([xunitINTsm ; yunitINTsm]);
    meanCirSMint{ccII}(1:lenCirSMint,3:4) = repmat(centroidsInt(ccII,[3 , 4]),lenCirSMint,1);
    
end

%% Compute Centered circle

centroidsCen = bsxfun(@minus, centroidsInt, min(centroidsInt));
centroidsCen(:,3) = 0:1:numel(centroidsCen(:,1))-1;

meanAEDia = nanmean(actEleDia);
meanCirCENSM = cell(size(centroidsCen,1),1);% xyz
for ccIII = 1:size(centroidsCen,1)
    
    [xunitINTCSM, yunitINTCSM] = SpecifyCircle_01(centroidsCen(ccIII,1),centroidsCen(ccIII,2),meanAEDia, numPts); % original centroid , actual electrode diameter
    
    % original centroid , original electrode diameter
    lenCirCENSM = size(xunitINTCSM,2);
    meanCirCENSM{ccIII}(1:lenCirCENSM,1:2) = transpose([xunitINTCSM ; yunitINTCSM]);
    meanCirCENSM{ccIII}(1:lenCirCENSM,3:4) = repmat(centroidsCen(ccIII,[3 , 4]),lenCirCENSM,1);
    
end




%% Derive mean circle

meanDia = nanmean(diameterS);
meanCir = cell(size(testIMS,3),1);% xyz

meanAECir = cell(size(testIMS,3),1);% xyz
meanSMCir = cell(size(testIMS,3),1);
for ci = 1:size(testIMS,3)
    
    if isnan(diameterS(ci,1))
        continue
    else
        
        [xunitMAE, yunitMAE] = SpecifyCircle_01(centroidS(ci,1),centroidS(ci,2),meanAEDia, numPts); % original centroid , actual electrode diameter
        [xunitM, yunitM] = SpecifyCircle_01(centroidS(ci,1),centroidS(ci,2),meanDia, numPts); % original centroid , original electrode diameter
        [xunitSM, yunitSM] = SpecifyCircle_01(centroidSM(ci,1),centroidSM(ci,2),meanAEDia, numPts); % smooth centroid , actual electrode diameter
        
        % original centroid , original electrode diameter
        lenCir = size(xunitM,2); 
        meanCir{ci}(1:lenCir,1:2) = transpose([xunitM ; yunitM]); 
        meanCir{ci}(1:lenCir,3) = ci;
        
        % original centroid , actual electrode diameter
        lenCirAE = size(xunitMAE,2);
        meanAECir{ci}(1:lenCirAE,1:2) = transpose([xunitMAE ; yunitMAE]);
        meanAECir{ci}(1:lenCirAE,3) = ci;
       
        % smooth centroid , actual electrode diameter
        lenCirSM = size(xunitSM,2);
        meanSMCir{ci}(1:lenCirSM,1:2) = transpose([xunitSM ; yunitSM]);
        meanSMCir{ci}(1:lenCirSM,3) = ci;
    end

end

%% Set up Solid

if ismember(solidType,[0 1 2])
    
    cirInds = find(imInd);
    rows = repmat(transpose(1:size(tmpIM,1)),[1,size(tmpIM,2),length(cirInds)]);
    cols = repmat(1:size(tmpIM,1),[size(tmpIM,2),1,length(cirInds)]);
    zVals = zeros(size(rows));
    
    for i = 1:length(cirInds)
        
        zVals(:,:,i) = cirInds(i);
        
    end
    
    cdata = zeros(size(rows));
    for i = 1:length(cirInds)
        
        switch solidType
            case 0 % original centroid , original electrode diameter
                tmpCir = meanCir{cirInds(i)};
            case 1 % original centroid , actual electrode diameter
                tmpCir = meanAECir{cirInds(i)};
            case 2 % smooth centroid , actual electrode diameter
                tmpCir = meanSMCir{cirInds(i)};
        end
        cdata(round(tmpCir(:,1)),round(tmpCir(:,2)),i) = 1;
        
    end
    
else
    
    rows = repmat(transpose(1:size(tmpIM,1)),[1,size(tmpIM,2),length(centroidsInt)]);
    cols = repmat(1:size(tmpIM,1),[size(tmpIM,2),1,length(centroidsInt)]);
    zVals = zeros(size(rows));
    
    for i = 1:length(centroidsInt)
        
        zVals(:,:,i) = i;
        
    end
    
    cdata = zeros(size(rows));
    for i = 1:length(centroidsInt)
        
        tmpCir = meanCirSMint{i};

        cdata(round(tmpCir(:,1)),round(tmpCir(:,2)),i) = 1;
        
    end
    
    
    
end


%% DATA
output_args = struct;
output_args.imInd = imInd;
output_args.blobDims = blobDims;
output_args.circleDims = circleDims;
output_args.centroidS = centroidS;
output_args.centroidsInt = centroidsInt;
output_args.centroidSM = centroidSM;
output_args.centroidsCen = centroidsCen;
output_args.diameterS = diameterS;
output_args.actEledims = actEledims;
output_args.actEleDia = actEleDia;
output_args.meanDia = meanDia;
output_args.meanCir = meanCir;
output_args.meanAEDia = meanAEDia;
output_args.meanAECir = meanAECir;
output_args.meanSMCir = meanSMCir;
output_args.meanCirSMint = meanCirSMint;
output_args.meanCirCENSM = meanCirCENSM;
output_args.rows = rows;
output_args.cols = cols;
output_args.zVals = zVals;
output_args.cdata = cdata;



end



function [xunit, yunit] = SpecifyCircle_01(x,y,r,numPts)

if nargin == 3
    numPts = 100;
end

steps = round(numPts/2);

th = 0:pi/steps:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;

end



function [ eQdiameter ] = getREgionDims(x_vec, y_vec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


minX = min(x_vec);
maxX = max(x_vec);

Xdist = maxX - minX;

minY = min(y_vec);
maxY = max(y_vec);

Ydist = maxY - minY;

eQdiameter = mean([Xdist , Ydist]);


end
















