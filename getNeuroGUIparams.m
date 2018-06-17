function [x,y,z,sz,col] = getNeuroGUIparams(neuroDat, output_args,solidCol,sizeMag)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

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

x.C = sXc;
x.A = sXa;
x.P = sXp;
x.M = sXm;
x.L = sXl;
y.C = sYc;
y.A = sYa;
y.P = sYp;
y.M = sYm;
y.L = sYl;
z.C = sZc;
z.A = sZa;
z.P = sZp;
z.M = sZm;
z.L = sZl;
sz.C = cSizeFR2;
sz.A = aSizeFR2;
sz.P = pSizeFR2;
sz.M = mSizeFR2;
sz.L = lSizeFR2;
col.C = colormapC2;
col.A = colormapA2;
col.P = colormapP2;
col.M = colormapM2;
col.L = colormapL2;


end

