% performs cross correlation of each row and column with itself, 
% assuming mirror symetry, the center at each axes is the maximum of the
% histogram of the cross correlation.
% 
%
%
function [rowsOfMaxes, colsOfMaxes] = FindCentre(A)

[m, n] = size(A);

% colomn
I = zeros(1, m);
for i = 1 : m
    vec = A(i, :);
    croscorr = xcorr(vec, vec(end:-1:1));         % xcorr of colomn with itself
    [~, I(i)] = max(croscorr);                    % find maximum pixel
end

[N,edges] = histcounts(I, 0:(2*m-1));
edges = edges(2:end);
N = N(edges>5 & edges<(2*m-6));                   % removal of eadges
[~, maxN_Inx] = max(N);
rowsOfMaxes = I(abs(I - maxN_Inx - 5)<10);        % averaging around the maximum value
rowsOfMaxes = mean(rowsOfMaxes)/2 + 0.5;          % adding half a pixel due to bias

% row
J = zeros(1, n);
for j=1:n
    vec = A(:, j);
    croscorr = xcorr(vec, vec(end:-1:1));        % xcorr of row with itself
    [~, J(j)] = max(croscorr);                   % find maximum pixel
end

[N,edges] = histcounts(J, 0:(2*n-1));
edges = edges(2:end);
N = N(edges>5 & edges<(2*n-6));                  % removal of eadges
[~, maxN_Inx] = max(N);
colsOfMaxes = J(abs(J - maxN_Inx - 5)<10);       % averaging around the maximum value
colsOfMaxes = mean(colsOfMaxes)/2 + 0.5;         % adding half a pixel due to bias

end
