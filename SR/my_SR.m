% this function and all its' reffered function were taken from
% http://www.ee.ucsc.edu/~milanfar from which we used the SR section
%
function HR = my_SR(LR_vec, resFactor, Hpsf, props)

D = RegisterImageSeq(LR_vec);

% Round translation to nearest neighbor
D = round(D.*resFactor);

% Shift all images so D is bounded from 0-resFactor
% Dr = floor(D/resFactor);

ind = sign(D);

if not((all(ind(:, 1) >= 0) | all(ind(:, 1) <= 0)) & (all(ind(:, 2) >= 0) | all(ind(:, 2) <= 0)))
    [x, y] = find(ind<0);
    for i=1:length(x)
        if y(i) == 2
            LR_vec(1:end-1, :, x(i)) = LR_vec(2:end, :, x(i));
        else
            LR_vec(:, 1:end-1, x(i)) = LR_vec(:, 2:end, x(i));
        end
    end
    D = RegisterImageSeq(LR_vec);
    D = round(D.*resFactor);
end
D = mod(D, resFactor) + resFactor;

if nargin<4
    props.alpha = 0.7;
    props.beta = 0.01;
    props.lambda = 0.01;
    props.P = 6;
    props.maxIter = 100;
end


HR = FastRobustSR(LR_vec(3:end-2,3:end-2,:), D, resFactor, Hpsf, props);
end
