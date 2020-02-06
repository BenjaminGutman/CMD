% given I(q) with only two circles, this function calculates the distance
% between the outer peak and the valey and devides by the ditrances of the
% extrimum in 1 space
% 
% input: known_separation_between_x_axis - a rough estimation of the
% distance between peaks
%
%
function score = separation_criterion(q, I, known_separation_between_x_axis)

% handles if I is not normilized
if max(I)>1
    I = I./max(I);
end

% identifies peaks
[pks,locs] = findpeaks(I, q, 'MinPeakDistance', known_separation_between_x_axis,...
    'MinPeakHeight', 0.3);
[~, inx] = sort(pks, 'descend');

% in case a single peak is found
if length(pks)<2
    score = nan;
    return
end

% lower peak
SecondPeakX = locs(inx(2));
SecondPeakY = pks(inx(2));

% valey
ValeyInx = find(q > locs(inx(1)) & q < locs(inx(2)));
ValeyX = q(ValeyInx);
ValeyY = I(ValeyInx);

% valey minima
[ValeySpotY, ValeySpotInx] = min(ValeyY);
ValeySpotX = ValeyX(ValeySpotInx);

% delta
score = (SecondPeakY - ValeySpotY) / (SecondPeakX - ValeySpotX);

end