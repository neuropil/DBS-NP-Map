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
