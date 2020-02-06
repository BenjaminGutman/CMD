
%% options
options = struct;
options.gridSz          = 160;                            % size of image in pixels 
options.BKG.A           = 10;                             % BKG amplifier
options.BKG.r           = -1.7;                           % BKG power law
options.CirclesParams   = [30, 0.7, 1;35, 0.5, 0.5];      % [radius , width, Amplifier ;...
options.PSF.params      = [0.1 1; 1.6 1.5; 3.1 2; 4.6 2.5]; % [width, erf incline ;...
options.PSF.sz          = 19;                             % in pixels
options.ET              = 0.6;                            % exposure timr
options.LRfactor        = 3;                              % denoted as f in the artical
options.separation_init = 4;                              % used in 'separation_criterion' function

options.PSF.number = size(options.PSF.params, 1)^2;
options.BeamStopSz = min(options.CirclesParams(:, 1))-3;  % smallest circle radius - 5 pixel

%% GT creation
GT = simulate_GT(options);                                % simulates ground truch circles
GT = GT + BKG_noise(options);                             % adds background noise
GT = GT./max(GT, [], 'all');                              % normalizing
[IGT, qGT, ~] = CircularIntegration(GT);                  % circular int of GT
% plot(q, vec)
scoreGT = separation_criterion(qGT, IGT, options.separation_init); % delta of GT

%% PSF + blur
PSF = MultiplePSFGenerator(options);                      % generation of PSFs
corrupted = multipleBlur(GT, PSF, options);               % bluring and addition of poisson noise
LRImages = Decimation(corrupted, options);                % decimation for low rasolution image

%% CMD
I0 = conv2(corrupted{1}, zeros(size(PSF{1})));            %initial guess
flux = zeros(options.PSF.number);                         % flux at each corrupted image for sigma
for ik = 1:options.PSF.number
    flux(ik) = sum(corrupted{ik}, 'all');
end
S = 0.2;                                                  % sigmas mean parameter
lambda = 0.05;   
sigmas = normpdf(flux, mean(flux)*S, mean(flux)/3);
sigmas = sigmas./max(sigmas);                             % normalization of sigma
Rec = CMD(corrupted, PSF, I0, eps, 200, 0, lambda,1./sigmas); % main algorithm
Rec = gather(Rec); Rec(Rec<0)=0;                          % reconstruction final adhustment
[I, q, ~] = CircularIntegration(Rec);                     % circular int. of reconstruction
score = separation_criterion(q, I, options.separation_init); % delta of reconstruction

