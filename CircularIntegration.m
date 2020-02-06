% suppose A is a matrix with circular symmetry,  
% this function performs uses the center to evaluate the intensity at each
% distance from the center
%
% input: A - imgae
%        center - of image, optional
%        LR_factor - if A is decimated, the LR_factor may be used to start
%        the integration at some initial radius, optional
%
% output: I - integral at each q
%         q - the used radii
%         STD - stardant deviation of intensity for each q
%
% options: Initial_q - niminal q to integrate
%          StepSize - integration quanta
%
%
%
%

function [I, q, STD] = CircularIntegration(A, center, LR_factor)

%% parameter
Initial_q = 2;                  % of which the algoruthm starts integrating.
StepSize = .5;                  % integration quanta.

%% haldel bags
% handles in case resolution enhancment factor is not specified 
if ~exist('LR_factor', 'var')
    LR_factor = 1;
end

% handles in case center is not provided 
if ~exist('center', 'var')
    [center(1), center(2)] = FindCentre(A);
end

% handles the case where A is not valid
if any([sum(abs(A(:)))==0 isnan(sum(abs(A(:))))])
    q = 1:10;
    I = nan(size(q));
    STD = I;
    return
end

%% algorithm
% pixel distance from center
[x,y] = meshgrid(1:size(A, 1), 1:size(A, 2));
r = sqrt((x-center(1)).^2+(y-center(2)).^2);
R = round(max(center));

I = zeros(1,R);
STD = I;
q = I;
i = 1;
inner = Initial_q/LR_factor;
outer = inner + StepSize;

while outer<R
    inx = find(r<outer & r>=inner);
    I(i) = mean(A(inx), 'omitnan');
    STD(i) = std(A(inx), 'omitnan');
    q(i) = inner;
    i = i + 1;
    inner = outer;
    outer = outer + StepSize;

end


end

