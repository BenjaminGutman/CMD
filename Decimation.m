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
LRImages = zeros(floor(options.gridSz/LRF), floor(options.gridSz/LRF),...
    LRF^2, options.PSF.number);

lastpxl = floor(options.gridSz/LRF)*LRF;            % ensures that decimates images are of same size

for i=1:options.PSF.number
    inx = 1;
    for j = 1:LRF
        for k = 1:LRF
            LRImages(:, :, inx, i) = OrigImages{i}(j:LRF:lastpxl, k:LRF:lastpxl);
            inx = inx + 1;
        end
    end
end