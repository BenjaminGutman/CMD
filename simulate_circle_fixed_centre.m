% creates image with a lorenzian profiled circle
% input: r0 - radius of circle
%        sigma - circle width
%        gridSz - of image
function Final = simulate_circle_fixed_centre(r0, sigma, gridSz)

if r0*2+5*sigma>gridSz
    warning('radius larger than image size')
end
t = linspace(-10,10,gridSz);
[x,y] = meshgrid(t);
A = sqrt((x).^2 + (y).^2);

A = A.*5;
r0 = r0/gridSz*100;
Final = 1/(pi)*sigma./((A - r0).^2 + sigma^2);
Final = Final - max(Final, [], 'all')/100;
Final(Final<0) = 0;
Final = Final/max(Final(:));


end