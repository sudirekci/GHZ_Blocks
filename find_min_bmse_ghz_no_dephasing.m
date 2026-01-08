

function [min_bmse, m1, m2] = find_min_bmse_ghz_no_dephasing(coeffs, delta_phi, m1, m2)

    no_atoms = length(coeffs)-1;

    rho = coeffs'*coeffs;

%     exp_factors = exp(-(0:no_atoms).^2*delta_phi^2/2);
% 
%     exponential_matrix = toeplitz(exp_factors);
% 
%     rho_bar = rho.*exponential_matrix;
% 
%     c = (0:-1:-no_atoms).*exp_factors;
%     r = (0:no_atoms).*exp_factors;
%     exponential_matrix = toeplitz(c, r);
% 
%     rho_bar_prime = 1i*delta_phi^2*rho.*exponential_matrix;
    
    rho_bar = rho.*m1;
    rho_bar_prime = rho.*m2;

    [V,D] = eig(rho_bar);
    D = diag(D);

    YMatrix = 2*V'*rho_bar_prime*V;

    for ii=1:no_atoms+1

        for jj=1:no_atoms+1

            YMatrix(ii, jj) = YMatrix(ii, jj)/(D(ii)+D(jj));

        end

    end

%     LOpt = V*YMatrix*V';
% 
%     min_bmse = delta_phi^2 - trace(rho_bar*LOpt*LOpt);
%     
%     LOpt = V*YMatrix*V';

    min_bmse = delta_phi^2 - sum(D.*diag(YMatrix*YMatrix));

end
