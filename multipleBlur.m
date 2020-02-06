% performs convolutional blur of ground truth adn the different PSFs and
% adds poiison noise.
% options used: 
%           options.PSF.number
%           options.ET - exposure time equivalent
%
%
%
function corrupted = multipleBlur(GT, PSF, options)

corrupted = cell(1, options.PSF.number);
for ik = 1:options.PSF.number
    corrupted{ik} = single(poissrnd(options.ET*conv2(GT, PSF{ik}, 'same')));
end

end