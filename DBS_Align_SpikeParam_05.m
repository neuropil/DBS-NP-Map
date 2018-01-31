function [] = DBS_Align_SpikeParam_05(ele_nii, mr_nii ,...
    sliceCnum, sliceSnum , sliceAnum , solidCol , sizeMag , bubBorder, neuroDATcsv)

% PLOT GRAY BEST FIT LINE
% PLOT DBS CONTACTS
% DBS_Align_SpikeParam_01
% OLD Plot_MW_DBSLinearBubble_v1

figure;
[ output_args ] = ExtractDBSLeadPoly_02( ele_nii  , 3 , 1.3 , 80);

mriLoad = load_nii(mr_nii); % NIFTI TOOLS

[brainIm] = Process_MRI_01(mriLoad);

[Xsl , Ysl , Zsl] = size(brainIm);

% sliceNum = 279;
% sliceO = 'S';
hold on
if isnan(sliceCnum)
    sliceCu = round(Xsl/2);
else
    sliceCu = sliceCnum;
end
c = slice(brainIm,sliceCu,[],[]);
squzDataC = squeeze(double(brainIm(:,sliceCu,:)));
set(c, 'alphadata', squzDataC, 'facealpha','interp');alim([0 0.5]);

if isnan(sliceSnum)
    sliceSu = round(Ysl/2);
else
    sliceSu = sliceSnum;
end
s = slice(brainIm,[],sliceSu,[]);
squzDataS = squeeze(double(brainIm(sliceSu,:,:)));
set(s, 'alphadata', squzDataS, 'facealpha','interp');alim([0 0.5]);

if isnan(sliceAnum)
    sliceAu = round(Zsl/2);
else
    sliceAu = sliceAnum;
end
a = slice(brainIm,[],[],sliceAu);
squzDataA = squeeze(double(brainIm(:,:,sliceAu)));
set(a, 'alphadata', squzDataA, 'facealpha','interp');alim([0 0.5]);

shading('interp')
colormap('bone')

hold on

Zticks = get(gca,'ZTick');
sliceThick = mriLoad.hdr.dime.pixdim(4);
ZticksSlice = Zticks*sliceThick;
z2cell = num2cell(ZticksSlice);
z2str = cellfun(@(x) num2str(x), z2cell, 'UniformOutput', false);
zticklabels(z2str)

hold on

[~] = Plot3D_ElecBound_01( output_args );

% load('mkAACAData.mat','allNeurons')


%% Get Neuron Data
% frRates = cellfun(@(x) x.FR, allNeurons);
% depthSS = cellfun(@(x) x.Depth.Actual, allNeurons);
% eleN = cellfun(@(x) str2double(x.CaseInfo.electrode), allNeurons);
% eleI1 = cellfun(@(x) x.trackIDS, allNeurons, 'UniformOutput',false);
% 
% cutThr = mean(frRates) + (std(frRates)*3);
% cutInd = frRates >= cutThr;
% 
% 
% allNeurons = allNeurons(~cutInd);

neuroDat = readtable(neuroDATcsv);

%%
[Xc, Yc, Zc, idC,...
    Xa, Ya, Za, idA, Xp,...
    Yp, Zp, idP, Xm, Ym,...
    Zm, idM, Xl, Yl, Zl,...
    idL, featureOut] = DeriveXYZ_neurOverlay_03( neuroDat ,...
    output_args , 2, 100);

%% Compute Feature coordinates

cNan = ~isnan(Xc);
Xc = Xc(cNan);
Yc = Yc(cNan);
Zc = Zc(cNan);
allc = idC(cNan);

aNan = ~isnan(Xa);
Xa = Xa(aNan);
Ya = Ya(aNan);
Za = Za(aNan);
alla = idA(aNan);

pNan = ~isnan(Xp);
Xp = Xp(pNan);
Yp = Yp(pNan);
Zp = Zp(pNan);
allp = idP(pNan);

mNan = ~isnan(Xm);
Xm = Xm(mNan);
Ym = Ym(mNan);
Zm = Zm(mNan);
allm = idM(mNan);

lNan = ~isnan(Xl);
Xl = Xl(lNan);
Yl = Yl(lNan);
Zl = Zl(lNan);
alll = idL(lNan);

%%

if ~solidCol
    
    colormapC = [transpose(linspace(226/255,74/255,length(Xc))) ,...
        transpose(linspace(0/255,142/255,length(Xc))) ,...
        transpose(linspace(12/255,203/255,length(Xc)))];
    
    colormapA = [transpose(linspace(226/255,74/255,length(Xa))) ,...
        transpose(linspace(0/255,142/255,length(Xa))) ,...
        transpose(linspace(12/255,203/255,length(Xa)))];
    
    colormapP = [transpose(linspace(226/255,74/255,length(Xp))) ,...
        transpose(linspace(0/255,142/255,length(Xp))) ,...
        transpose(linspace(12/255,203/255,length(Xp)))];
    
    colormapM = [transpose(linspace(226/255,74/255,length(Xm))) ,...
        transpose(linspace(0/255,142/255,length(Xm))) ,...
        transpose(linspace(12/255,203/255,length(Xm)))];
    
    colormapL = [transpose(linspace(226/255,74/255,length(Xl))) ,...
        transpose(linspace(0/255,142/255,length(Xl))) ,...
        transpose(linspace(12/255,203/255,length(Xl)))];
    
else
    colormapC = repmat([1 0 1],length(Xc),1);
    colormapA = repmat([1 1 0],length(Xa),1);
    colormapP = repmat([1 1 1],length(Xp),1);
    colormapM = repmat([1 1 1],length(Xm),1);
    colormapL = repmat([0 1 0],length(Xl),1);
end

featureOut = featureOut(~isnan(featureOut));

% Firing Rate
allFRSizeSpace = round(linspace(min(featureOut),max(featureOut),5));
[~,~,binFRall] = histcounts(featureOut,5);

allFRSize = zeros(length(featureOut),1);
for fralli = 1:length(featureOut)
    allFRSize(fralli) = allFRSizeSpace(binFRall(fralli));
end

% Figure Setup

[~,sortOrCZ] = sort(Zc);
[~,sortOrAZ] = sort(Za);
[~,sortOrPZ] = sort(Zp);
[~,sortOrMZ] = sort(Zm);
[~,sortOrLZ] = sort(Zl);

sXc = Xc(sortOrCZ);
sYc = Yc(sortOrCZ);
sZc = Zc(sortOrCZ);

sXa = Xa(sortOrAZ);
sYa = Ya(sortOrAZ);
sZa = Za(sortOrAZ);

sXp = Xp(sortOrPZ);
sYp = Yp(sortOrPZ);
sZp = Zp(sortOrPZ);

sXm = Xm(sortOrMZ);
sYm = Ym(sortOrMZ);
sZm = Zm(sortOrMZ);

sXl = Xl(sortOrLZ);
sYl = Yl(sortOrLZ);
sZl = Zl(sortOrLZ);

cSizeFR = featureOut(allc(sortOrCZ));
aSizeFR = featureOut(alla(sortOrAZ));
pSizeFR = featureOut(allp(sortOrPZ));
mSizeFR = featureOut(allm(sortOrMZ));
lSizeFR = featureOut(alll(sortOrLZ));

[~,colOrderC] = sort(cSizeFR);
[~,colOrderA] = sort(aSizeFR);
[~,colOrderP] = sort(pSizeFR);
[~,colOrderM] = sort(mSizeFR);
[~,colOrderL] = sort(lSizeFR);

colormapC2 = colormapC(colOrderC,:);
colormapA2 = colormapA(colOrderA,:);
colormapP2 = colormapP(colOrderP,:);
colormapM2 = colormapM(colOrderM,:);
colormapL2 = colormapL(colOrderL,:);
        
if sizeMag
    
    cSizeFR2 = cSizeFR;
    aSizeFR2 = aSizeFR;
    pSizeFR2 = pSizeFR;
    mSizeFR2 = mSizeFR;
    lSizeFR2 = lSizeFR;

else
    
    cSizeFR2 = 60;
    aSizeFR2 = 60;
    pSizeFR2 = 60;
    mSizeFR2 = 60;
    lSizeFR2 = 60;

end

hold on
sC = scatter3(sXc , sYc , sZc, cSizeFR2, colormapC2, 'filled');
sA = scatter3(sXa , sYa , sZa, aSizeFR2, colormapA2, 'filled');
sP = scatter3(sXp , sYp , sZp, pSizeFR2, colormapP2, 'filled');
sM = scatter3(sXm , sYm , sZm, mSizeFR2, colormapM2, 'filled');
sL = scatter3(sXl , sYl , sZl, lSizeFR2, colormapL2, 'filled');

grid off
set(gca,'XTickLabel',[])
set(gca,'YTickLabel',[])
set(gca,'ZTickLabel',[])
set(gca,'XTick',[])
set(gca,'YTick',[])
set(gca,'ZTick',[])
set(gca,'View', [-116.0000  22.0000])
set(gca,'Color','none')

if bubBorder
    
    sC.MarkerEdgeColor = 'k';
    sA.MarkerEdgeColor = 'k';
    sP.MarkerEdgeColor = 'k';
    sM.MarkerEdgeColor = 'k';
    sL.MarkerEdgeColor = 'k';
    
end



end

function [ mriOUT] = Process_MRI_01(MRI)
%Process_MRI converts values within MRI to 'single'
% Input is NIFTI Tools structure .img input

mriOUT = zeros(size(MRI.img),'single');

for mi = 1:size(mriOUT,3)
    
    mriIm = single(MRI.img(:,:,mi));
    mriIm = mriIm/(max(max(mriIm)));
    mriOUT(:,:,mi) = mriIm;
    
end



end














