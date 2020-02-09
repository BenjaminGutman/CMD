
addpath 'SR';
addpath 'SR\registration';
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

options.flag.SR_alone = 1;
options.flag.CMD_alone = 1;
options.flag.srSAXS = 1;
options.flag.plot = 1;
options.flag.enableGPU = 0;

options.PSF.number = size(options.PSF.params, 1)^2;
options.BeamStopSz = min(options.CirclesParams(:, 1))-3;  % smallest circle radius - 5 pixel

%% GT creation
GT = simulate_GT(options);                                % simulates ground truch circles
GT = GT + BKG_noise(options);                             % adds background noise
GT = GT./max(GT, [], 'all');                              % normalizing
[IGT, qGT, ~] = CircularIntegration(GT);                  % circular int of GT
deltaGT = separation_criterion(qGT, IGT, options.separation_init); % delta of GT

%% PSF + blur
PSF = MultiplePSFGenerator(options);                      % generation of PSFs
corrupted = multipleBlur(GT, PSF, options);               % bluring and addition of poisson noise
LRImages = Decimation(corrupted, options);                % decimation for low rasolution image

%% SR
% SR as first step of srSAXS
if options.flag.srSAXS            
    HR_for_CMD = cell(1, options.LRfactor);
    f = waitbar(0, 'SR in progress');
    
    for i = 1:options.PSF.number
        f = waitbar(i/options.PSF.number, f ,'SR in progress');
        HR_for_CMD{i} = my_SR(squeeze(LRImages(:, :, :, i)), options.LRfactor, 1);
    end
    close(f)
    
end

% SR alone
if any([options.flag.SR_alone  options.flag.srSAXS])  
    f = waitbar(0, 'SR in progress');
    HR_alone = cell(1, options.LRfactor);
    I_SR_alone = HR_alone;
    q_SR_alone = HR_alone;
    deltaSRalone = zeros(1, options.LRfactor);
    
    for i = 1:options.PSF.number
        f = waitbar(i/options.PSF.number, f ,'SR in progress');
        HR_alone{i} = my_SR(squeeze(LRImages(:, :, :, i)), ...
            options.LRfactor, PSF{i});
        [I_SR_alone{i}, q_SR_alone{i}, ~] = ...
            CircularIntegration(HR_alone{i});                     % circular int. of reconstruction
        deltaSRalone(i) = separation_criterion(q_SR_alone{i}, ...
            I_SR_alone{i}, options.separation_init/options.LRfactor);           % delta of reconstruction
    end
    close(f)
    
end

%% CMD

% CMD alone
if options.flag.CMD_alone
    % CMD params
    S = 0.3;                                                 % sigmas mean parameter
    lambda = 0.05;   

    I0 = conv2(corrupted{1}, zeros(size(PSF{1})));            % initial guess
    flux = zeros(options.PSF.number);                         % flux at each corrupted image for sigma
    for ik = 1:options.PSF.number
        flux(ik) = sum(corrupted{ik}, 'all');
    end

    sigmas = normpdf(flux, mean(flux)*S, mean(flux)/3);
    sigmas = sigmas./max(sigmas);                             % normalization of sigma

    CMD_alone_Rec = CMD(corrupted, PSF, I0, eps, 200, 0, lambda,1./sigmas, options.flag.enableGPU); % main algorithm

	CMD_alone_Rec(CMD_alone_Rec<0)=0;
    if options.flag.enableGPU
        CMD_alone_Rec = gather(CMD_alone_Rec);                           % reconstruction final adjustment
    end
    [I_CMD_alone_Rec, q_CMD_alone_Rec, ~] = CircularIntegration(CMD_alone_Rec);                     % circular int. of reconstruction
    deltaCMDalone = separation_criterion(q_CMD_alone_Rec, ...
        I_CMD_alone_Rec, options.separation_init);            % delta of reconstruction
end


% CMD as second step of srSAXS
if options.flag.srSAXS
    S = 0.3;                                                 % sigmas mean parameter
    lambda = 0.05;   

    I0 = conv2(HR_for_CMD{1}, zeros(size(PSF{1})));           % initial guess
    flux = zeros(options.PSF.number);                         % flux at each corrupted image for sigma
    for ik = 1:options.PSF.number
        flux(ik) = sum(HR_for_CMD{ik}, 'all');
    end

    sigmas = normpdf(flux, mean(flux)*S, mean(flux)/3);
    sigmas = sigmas./max(sigmas);                             % normalization of sigma

    srSAXS_Rec = CMD(HR_for_CMD, PSF, I0, eps, 200, 0, lambda,1./sigmas, options.flag.enableGPU); % main algorithm
	srSAXS_Rec(srSAXS_Rec<0)=0;
    if options.flag.enableGPU
        srSAXS_Rec = gather(srSAXS_Rec);                           % reconstruction final adjustment
    end
    [I_srSAXS, q_srSAXS, ~] = CircularIntegration(srSAXS_Rec);                % circular int. of reconstruction
    deltasrSAXS = separation_criterion(q_srSAXS, I_srSAXS, options.separation_init); % delta of reconstruction
end
%% plot 
if options.flag.plot
    figure('WindowStyle', 'docked', 'NumberTitle', 'off', 'Name', 'results')
    hold on; grid on; grid minor
    plot(qGT, IGT./max(IGT), 'DisplayName', sprintf('GT  %.2f', deltaGT))
        if options.flag.CMD_alone
            plot(q_CMD_alone_Rec, I_CMD_alone_Rec./max(I_CMD_alone_Rec), ...
                'DisplayName', sprintf('CMD reco  %.2f', deltaCMDalone))
        end
        if options.flag.srSAXS
            plot(q_srSAXS, I_srSAXS./max(I_srSAXS), ...
                'DisplayName', sprintf('srSAXS reco  %.2f', deltasrSAXS))
        end
        if options.flag.SR_alone
            plot(q_SR_alone{1}, I_SR_alone{1}./max(I_SR_alone{1}), ...
                'DisplayName', sprintf('SR reco  %.2f', deltaSRalone(1)))
        end
    l = legend;
    l.Title.String = 'numbers for delta critiria';
end