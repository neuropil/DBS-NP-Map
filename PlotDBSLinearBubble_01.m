function [plotDATA] = PlotDBSLinearBubble_01(ele_nii, mr_nii, caseID , zoomIN, sliceThick,...
    sliceNum, sliceO, ColorFR)

% PLOT GRAY BEST FIT LINE
% PLOT DBS CONTACTS
% EXCLUDE ALL BUT PATIENT
% FIND actual centroid and NOT BEST FIT

figure;
[ output_args ] = ExtractDBSLeadPoly_02( ele_nii  , 3 , 1.3 , 80);

nanInd = ~isnan(output_args.centroidSM(:,3));
Zvals2 = output_args.centroidSM(nanInd,3);

mriLoad = load_nii(mr_nii); % NIFTI TOOLS

[brainIm] = Process_MRI_01(mriLoad);

switch sliceO
    case 'C'
        s = slice(brainIm,sliceNum,[],[]);
        squzData = squeeze(double(brainIm(:,sliceNum,:)));
        set(gca,'View',[-270.7000 -0.4000])

    case 'S'
        s = slice(brainIm,[],sliceNum,[]);
        squzData = squeeze(double(brainIm(sliceNum,:,:)));
        set(gca,'View',[180 0])

    case 'A'
        s = slice(brainIm,[],[],sliceNum);
        squzData = squeeze(double(brainIm(:,:,sliceNum)));
        set(gca,'View',[180 90])

end
shading('interp')
colormap('gray')

set(s, 'alphadata', squzData, 'facealpha','interp');alim([0 0.2]);

hold on

Zticks = get(gca,'ZTick');
ZticksSlice = Zticks*sliceThick;
z2cell = num2cell(ZticksSlice);
z2str = cellfun(@(x) num2str(x), z2cell, 'UniformOutput', false);
zticklabels(z2str)

hold on

[~] = Plot3D_ElecBound_01( output_args );

load('MW_Fig1Data.mat','allNeurons')


%% Get Neuron Data
frRates = cellfun(@(x) x.FR, allNeurons); %#ok<NODEF>
cutThr = mean(frRates) + (std(frRates)*3);
cutInd = frRates >= cutThr;
allNeurons = allNeurons(~cutInd);
%% Get Patient cells
caseNumS = cellfun(@(x) x.caseNum, allNeurons);
caseInd = caseNumS == caseID;
allNeurons = allNeurons(caseInd);
feature = 1;
%%
[Xc, Yc, Zc, idC,...
    Xa, Ya, Za, idA, Xp,...
    Yp, Zp, idP, Xm, Ym,...
    Zm, idM, Xl, Yl, Zl,...
    idL, featureOut] = DeriveXYZ_neurOverlay_02( allNeurons,...
    feature,...
    output_args, 'Thalamus',sliceThick);

%% Compute Feature coordinates
if feature == 3
    idC = idC(~isnan(idC));
    idA = idA(~isnan(idA));
    idM = idM(~isnan(idM));
    idP = idP(~isnan(idP));
    idL = idL(~isnan(idL));
end

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

if ColorFR
    
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
    
    colormapA = repmat([1 1 0],length(Xc),1);
    
    colormapP = repmat([1 1 1],length(Xc),1);
    
    colormapM = repmat([1 1 1],length(Xc),1);
    
    colormapL = repmat([0 1 0],length(Xc),1);
    
    
end


%

featureOut = featureOut(~isnan(featureOut));

switch feature
    
    case 1
        % Firing Rate
        allFRSizeSpace = round(linspace(min(featureOut),max(featureOut),5));
        [~,~,binFRall] = histcounts(featureOut,5);
        
        allFRSize = zeros(length(featureOut),1);
        for fralli = 1:length(featureOut)
            allFRSize(fralli) = allFRSizeSpace(binFRall(fralli));
        end
        
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

switch feature
    case 1
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
end

%

switch feature
    
    case 1
        
        % CHANGE COLOR PER ELECTRODE
        hold on
        sC = scatter3(sXc , sYc , sZc, 60, colormapC2, 'filled');
        sA = scatter3(sXa , sYa , sZa, 60, colormapA2, 'filled');
        sP = scatter3(sXp , sYp , sZp, 60, colormapP2, 'filled');
        sM = scatter3(sXm , sYm , sZm, 60, colormapM2, 'filled');
        sL = scatter3(sXl , sYl , sZl, 60, colormapL2, 'filled');
        
        if zoomIN
            
            xlim([230 290])
            ylim([260 310])
            zlim([round(min(Zvals2))-1 15])
            
        end

        sC.MarkerEdgeColor = 'k';
        sA.MarkerEdgeColor = 'k';
        sP.MarkerEdgeColor = 'k';
        sM.MarkerEdgeColor = 'k';
        sL.MarkerEdgeColor = 'k';
        
        %         title('Firing Rate')
        set(gca,'ZLim',[1 35])
        set(gca,'XLim',[70 420])
        set(gca,'YLim',[120 400])
        grid off
        set(gca,'XTickLabel',[])
        set(gca,'YTickLabel',[])
        set(gca,'ZTickLabel',[])
        set(gca,'XTick',[])
        set(gca,'YTick',[])
        set(gca,'ZTick',[])
        %         set(gca,'View', [-116.0000  22.0000])
        set(gca,'Color','none')
        
        if ~zoomIN
            
            set(gca,'XColor','none')
            set(gca,'YColor','none')
            set(gca,'ZColor','none')
            
        end
        
        
end

plotDATA.scatCen = [sXc , sYc , sZc];
plotDATA.scatAnt = [sXa , sYa , sZa];
plotDATA.scatPost = [sXp , sYp , sZp];
plotDATA.scatMed = [sXm , sYm , sZm];
plotDATA.scatLat = [sXl , sYl , sZl];

plotDATA.scatCenC = colormapC2;
plotDATA.scatAntC = colormapA2;
plotDATA.scatPostC = colormapP2;
plotDATA.scatMedC = colormapM2;
plotDATA.scatLatC = colormapL2;


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








