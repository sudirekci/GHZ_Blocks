

bits = [4,2,1].';
n0 = 5;

n1 = floor(n0/2);
n2 = n0-n1;

deltaPhi = 0.7;

noSamples = 1e5;

dPhi = 12*deltaPhi/(2000);
phi = -6*deltaPhi:dPhi:6*deltaPhi;

priorFunc = (exp(-(phi.^2)/(2*deltaPhi^2))).';

phiPriorFunc = (phi.*exp(-(phi.^2)/(2*deltaPhi^2))).';

noAtoms = sum(bits.*n0)

probsY = zeros([n2+1, length(bits), length(phi)]);
probsX = zeros([n1+1, length(bits), length(phi)]);

for j=0:n1

    probsX(j+1,:,:) = nchoosek(n1, j)*(sin((bits.*phi)/2).^2).^(j).*(cos((bits.*phi)/2).^2).^(n1-j);

end

for j=0:n2

    probsY(j+1,:,:) = nchoosek(n2, j)*(sin((bits.*phi+pi/2)/2).^2).^(j).*(cos((bits.*phi+pi/2)/2).^2).^(n2-j);

end

tic

[bmse_list] = forward_pass(noSamples, bits, deltaPhi, n1, n2,...
                    probsX, probsY, priorFunc, phiPriorFunc);

toc

bmse = mean(bmse_list)


function [bmse] = forward_pass(noSamples, bits, ...
                        deltaPhi, n1, n2, probsX, probsY, ...
                            priorFunc, phiPriorFunc)

    noBranches = ((n1+1)*(n2+1))^length(bits);

    samples = normrnd(0,deltaPhi,[noSamples, 1]);

    bmse = zeros([noSamples, 1]);
    
    for ii = 1:noSamples

        rx = binornd(n1, sin((bits.*samples(ii))/2).^2)+1;
        ry = binornd(n2, sin((bits.*samples(ii) + pi/2)/2).^2)+1;

        opt_est = 1;
        
        for j = 1:length(bits)

            opt_est = opt_est.*probsX(rx(j), j, :).*probsY(ry(j), j, :);

        end

        opt_est = squeeze(opt_est);

        opt_est = sum(opt_est.*phiPriorFunc, "all")/sum(opt_est.*priorFunc, "all");

        bmse(ii) = (samples(ii)-opt_est)^2;

    end

end
