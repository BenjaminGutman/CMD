% for each blurry and noisy image, this function performs down sampeling in
% pixel space
% options used:
%           options.LRfactor - denoted as f in the artical
%           options.PSF.number
%     
%
%
%
function LRImages = Decimation(OrigImages, options)

LRF = options.LRfactor;
LRImages = cell(LRF^2, options.PSF.number);


for i=1:options.PSF.number
    inx = 1;
    for j = 1:LRF
        for k = 1:LRF
            LRImages{inx, i} = OrigImages{i}(j:LRF:end, k:LRF:end);
            inx = inx + 1;
        end
    end
end