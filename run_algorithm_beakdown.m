tic
%% options
options = struct;
options.gridSz          = 160;                   % size of image in pixels 
options.BKG.A           = 10;                      % BKG amplifier
options.BKG.r           = -1.7;             % BKG power law
options.CirclesParams   = [30, 0.7, 1;35, 0.7, 1];      % [radius , width, Amplifier ;...
options.PSF.params      = [0.1 1; 1.6 1.5; 3.1 2; 4.6 2.5];       % [width, erf incline ;...
options.PSF.sz          = 19;                        % in pixels
options.ET              = 1.2;
options.LRfactor        = 3;

options.PSF.number = size(options.PSF.params, 1)^2;


score = zeros(1, 10, 10);
scoreGT = zeros(1, 10);
I = cell(size(score));
q = I;
IGT = cell(size(scoreGT));
qGT = IGT;
inner = 5;
outer = 1:10;
for i = 1:length(inner)
    for j = 1:length(outer)
        options.CirclesParams(1, 1) = inner(i);
        options.CirclesParams(2, 1) = outer(j) + inner(i);
        options.BeamStopSz = min(options.CirclesParams(:, 1));   % smallest circle radius - 5 pixel
        %% GT creation
        GT = simulate_GT(options);
        GT = GT + BKG_noise(options);
        GT = GT./max(GT, [], 'all');
        [IGT{i, j}, qGT{i, j}, ~] = CircularIntegration(GT);
        % plot(q, vec)
        scoreGT(i, j) = separation_criterion(qGT{i, j}, ...
            IGT{i, j}./max(IGT{i, j}), outer(j)/2);
        %% PSF + blur
        PSF = MultiplePSFGenerator(options);
        for MC = 1:2

            corrupted = multipleBlur(GT, PSF, options);
            % LRImages = Decimation(corrupted, options);
            %% CMD
            I0 = conv2(corrupted{1}, zeros(size(PSF{1})));
            flux = zeros(options.PSF.number);
            for ik = 1:options.PSF.number
                flux(ik) = sum(corrupted{ik}, 'all');
            end
            S = 0.2 ;
            lambda = 0.05;
            sigmas = normpdf(flux, mean(flux)*S, mean(flux)/3);
            sigmas = sigmas./max(sigmas);
            temp = CMD(corrupted, PSF, ...
                I0, eps, 200, 0, lambda,1./sigmas);
            temp = gather(temp); temp(temp<0)=0;
            [I{i, j, MC}, q{i, j, MC}, ~] = CircularIntegration(temp);
            score(i, j, MC) = separation_criterion(q{i, j , MC}, ...
                I{i, j, MC}./max(I{i, j, MC}), outer(j)/2);
        end
        sprintf('%.0f %d', toc, j)
    end
end