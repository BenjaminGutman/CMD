% M is a matrix of mx3, where is number of circles required (typically 2)
% and the columns of M are: r0, sigma, Amplitude of the circles
% GT witout noise and background
function GT = simulate_GT(options)

GT = zeros(size(options.gridSz));
for i = 1:size(options.CirclesParams, 1)
    GT = GT + options.CirclesParams(i, 3) * ...
        simulate_circle_fixed_centre(options.CirclesParams(i, 1),...
        options.CirclesParams(i, 2), options.gridSz);
end
    
end