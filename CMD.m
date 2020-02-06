% Minimizes the sum over ||D(conv(Image,Kernel_i)-Measured_i)||_F^2 over all
% kernels and measured images to find the original image, in operator
% notation: Solve ||PI-M||^2
% Equivalent to solving P* D^2 P I=P* D^2 M 

function I = CMD(Measured, Kernels, I, tol, MaxIter, MinIter, Lambda, sigmas)

    % normilize PSF
    for i=1:length(Measured)
        Kernels{i} = gpuArray(Kernels{i}./sum(Kernels{i}(:))); 
    end
    
    % start by computing P*M
    adjPM = AdjointOperator(Measured, Kernels, sigmas);  %%adding the images together 

    % Solve P*PI=P*M using conjugate gradients
    temp = SquareOperator(I, Kernels, Lambda, sigmas);
    r = adjPM - temp;
    p = r;
    k = 0;
    err = inf;
    % f = waitbar(0, 'CMD in progress');

    % CG iterations
    while (err>tol) && (k<MaxIter) || (k<MinIter)
        PadjPp = SquareOperator(p, Kernels, Lambda, sigmas);
        alpha = sum2(r.^2)/sum2(p.*PadjPp);
        I = I + alpha*p;
        r_next = r - alpha*PadjPp;
        beta = sum2(r_next.^2)/sum2(r.^2);
        p = r_next + beta*p;
        r = r_next;
        k = k + 1;
        err = sum2(abs(r));
    %     waitbar(k/MaxIter, f, 'CMD in progress');
    end
    
% restore original size
croping = (size(I, 1) - size(Measured{1}, 1))/2;
I = I(croping + 1:end - croping, croping + 1:end - croping);

end

% Compute P X
function Outputs = DirectOperator(Image, Kernels)
    N = length(Kernels);
    Outputs = cell(N,1);
    for i = 1:N
        Outputs{i} = conv2(Image, Kernels{i}, 'valid');

    end
    
end

% Compute P* D^2 X
function Output = AdjointOperator(Images, Kernels, sigmas)
    N = length(Kernels);
    Output = gpuArray(conv2(Images{1}, zeros(size(Kernels{1}), 'single'), 'full'));
    for i = 1:N
        Output = Output + conv2(Images{i}, rot90((Kernels{i}),2), 'full')/(2*sigmas(i)^2);
    end
end

% Compute P* D^2 PX
function Output = SquareOperator(I, Kernels, Lambda, sigmas)
    Output = AdjointOperator(DirectOperator(I, Kernels), Kernels, ...
        sigmas) + Lambda*I;
end
    

function s = sum2(x)
    s = x*x';
    s = sqrt(trace(s));
end

