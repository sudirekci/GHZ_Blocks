
delta_phi_array = 10.^(-1.9:0.18:log10(3.5));
bmse_array = zeros(1, length(delta_phi_array));

tic

for i=1:length(delta_phi_array)

    deltaPhi = delta_phi_array(i);
    
    bits = [4,2,1];
    noBlocks = [2,5,8];
    
    noSamples = 1e6;
    
    noPhi = 1500;
    
    dPhi = 12*deltaPhi/(noPhi-1);
    phi = -6*deltaPhi:dPhi:6*deltaPhi;
    
    priorFunc = 1/sqrt(2*pi*deltaPhi^2).*exp(-phi.^2./(2*deltaPhi^2));
    phiPriorFunc = phi.*priorFunc;
    
    noAtoms = sum(bits.*noBlocks);
    
    bitsRepeated = repelem(bits, noBlocks).';
    
    systemPhases = zeros(length(bitsRepeated),1);
    
    ind = 1;
    
    for jj=1:length(noBlocks)
    
        systemPhases(ind:ind+noBlocks(jj)-1,:) = (0:noBlocks(jj)-1)*pi/noBlocks(jj);
    
        ind = ind+noBlocks(jj);
    
    end
    
    probsPlus = sin((bitsRepeated.*phi - systemPhases)/2).^2;
    probsMinus = cos((bitsRepeated.*phi - systemPhases)/2).^2;
    
    
    [bmse, ~] = forward_pass(noSamples, bitsRepeated, systemPhases, deltaPhi, ...
        probsPlus, probsMinus, priorFunc, phiPriorFunc);

    bmse_array(i) = bmse;

end

toc

figure
loglog(delta_phi_array, sqrt(bmse_array)./delta_phi_array, 'o', ...
    'MarkerFaceColor', 'blue')


function [bmse, systemPhases] = forward_pass(noSamples, bitsRepeated, ...
    systemPhases, deltaPhi, probsPlus, probsMinus, priorFunc, phiPriorFunc)

    noBranches = 2^(length(bitsRepeated));

    samples = normrnd(0,deltaPhi,[noSamples, 1]);
   
    bmse = zeros([noSamples, 1]);
    
    for ii = 1:noSamples

        r = binornd(1, sin((bitsRepeated.*samples(ii) - systemPhases)/2).^2);

        opt_est = prod(r.*probsPlus + (1-r).*probsMinus,1);

        opt_est = sum(opt_est.*phiPriorFunc)/sum(opt_est.*priorFunc);

        bmse(ii) = (samples(ii)-opt_est)^2;

    end

    bmse = sum(bmse)/noSamples;

end
