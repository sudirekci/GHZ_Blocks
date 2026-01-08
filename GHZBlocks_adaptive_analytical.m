

bmse_threshold = 0.00005;

deltaPhi = 0.1;


bits = [16 8 4 2 1];
noBlocks = [1 1 0 0 2];

noAtoms = sum(bits.*noBlocks)

dPhi = 12*deltaPhi/2500;
phiArray = -6*deltaPhi:dPhi:6*deltaPhi;
priorFunc = exp(-(phiArray.^2)/(2*deltaPhi^2))/sqrt(2*pi*deltaPhi^2)*dPhi;
phiPriorFunc = phiArray.*exp(-(phiArray.^2)/(2*deltaPhi^2))/sqrt(2*pi*deltaPhi^2)*dPhi;

evaluate_derivatives = true;

initial_phases = rand([2^sum(noBlocks)-1, 1])*2*pi-pi;

[bmse, systemPhases] = optimize_adaptive(bits, noBlocks, deltaPhi, ...
       bmse_threshold, evaluate_derivatives, initial_phases, phiArray,...
       priorFunc, phiPriorFunc);


function [bmse, systemPhases] = optimize_adaptive(bits, noBlocks, deltaPhi, ...
    bmse_threshold, evaluate_derivatives, systemPhases, phiArray, priorFunc, ...
    phiPriorFunc)

    epsilon_adam = 10^(-8.);
    beta1 = 0.9;
    beta2 = 0.999;

    stepSize = 0.008; 

    bitsRepeated = repelem(bits, noBlocks).';

    noBranches = 2^sum(length(bitsRepeated));

    m = zeros([noBranches-1, 1]);
    v = zeros([noBranches-1, 1]);

    bmse_1 = -1;
    bmse_2 = -3;

    cont = true;
    increase_count = 0;

    while abs((bmse_2-bmse_1)/bmse_1)> bmse_threshold && cont

        bmse_1 = bmse_2;

        [systemPhases, derivatives, bmse_2] = forward_pass(systemPhases, ...
            bitsRepeated, noBlocks, phiArray, priorFunc, phiPriorFunc,deltaPhi);

        disp(bmse_2);

        if bmse_2 > bmse_1

            increase_count = increase_count + 1;

        end

        [systemPhases, m, v, ~] = take_step(systemPhases, derivatives, ...
            m, v, beta1, beta2, epsilon_adam, stepSize);

        cont = evaluate_derivatives;

    end

    bmse = bmse_2;

end



function [systemPhases, derivatives, bmse] = ...
    forward_pass(systemPhases, bitsRepeated, noBlocks, ...
    phiArray, priorFunc, phiPriorFunc,deltaPhi)

    noBranches = 2^(sum(noBlocks));

    derivatives = zeros([noBranches-1, 1]); 

    no_digits = sum(noBlocks);

    denoms = 2.^(no_digits-1:-1:0);
    nums = (2.^(1:no_digits-1));

    ind = ones(no_digits,1);

    bmse = deltaPhi^2;

    for ii = 0:noBranches-1

        div = floor(ii./denoms);
    
        tan_cot_indicator = (mod(div,2)).';

        ind(2:end) = div(1:end-1) + nums;

        phases = systemPhases(ind);

        tans = tan((bitsRepeated.*phiArray-phases)/2);
        cots = 1./tans;

        branchProbs = prod(sin((bitsRepeated.*phiArray-phases)/2).^2.*...
            (tan_cot_indicator)+ cos((bitsRepeated.*phiArray-phases)/2).^2.*...
            (1-tan_cot_indicator), 1);

        est_num = sum(branchProbs.*phiPriorFunc);
        est_denom = sum(branchProbs.*priorFunc);

        opt_estimator = est_num/est_denom;

        bmse = bmse - est_num^2/est_denom;
    
        derivatives(ind) = derivatives(ind) - sum(-2*opt_estimator.*...
            ((1-tan_cot_indicator).*(-tans)+(tan_cot_indicator).*cots).*...
            branchProbs.*phiPriorFunc + opt_estimator^2.*...
            ((1-tan_cot_indicator).*(-tans) + (tan_cot_indicator).*...
            (cots)).*priorFunc.*branchProbs, 2);

    end

end



function [systemPhases, m, v, derivatives] = take_step(systemPhases, ...
    derivatives, m, v, beta1, beta2, epsilon_adam, stepSize)

    %disp(derivatives);

    m = beta1*m + (1 - beta1)*derivatives;
    v = beta2*v + (1 - beta2)*derivatives.^2;

    systemPhases = systemPhases - stepSize*(m./(1 - beta1))./(sqrt(v./(1 - beta2)) + epsilon_adam);

end
