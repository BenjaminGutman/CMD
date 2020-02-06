% fucntion for creating 2D PSF with double error function profile
%
% input: 
%       a and b - 1x2 vectors each for the parameters of x and y axes
%       respectivly, 
%       PSFsz - size of the PSF image
%       f, g, ax1, ax2 - optional input to plot the images of the PSF (f, ax1) and ints
%       profile (g, ax2)
%
function PSF = generatePSF(a, b, PSFsz, f, g, ax1, ax2)


x = -(PSFsz-1)/2:(PSFsz-1)/2;
fun = @(x, a) erf((a(1)+x)/a(2)) + erf((a(1)-x)/a(2));
PSF = fun(x, b).*fun(x, a)';
PSF = PSF/max(PSF(:));

%% optional ploting
if nargin==3
    return
end
figure(f);
subplot(ax1);
imagesc(PSF)
title(num2str(sum(PSF(:))));
axis equal
%  colormap

figure(g);
subplot(ax2);
hold on; grid on; grid minor
y = fun(x, b);
plot(x, y/max(y));
y = fun(x, a);
plot(x, y/max(y));

end