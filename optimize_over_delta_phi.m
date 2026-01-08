

NoAtoms = 26;

delta_phi_list = 10.^(-2.:0.03:log10(3.2));

%delta_phi_list = 10.^(-1.9:0.18:log10(3.5));

bmse_min_list = zeros(length(delta_phi_list),1);

not_found = true;
k = -1;

while not_found
    
    k = k + 1;
    
    if 2^k > NoAtoms
        
        not_found = false;
        
    end
    
end

K = k-1;

partition_list = zeros(10000, K+1);

[num_partitions, partition_list] = find_partitions(NoAtoms, K, partition_list, 1);

partition_list = partition_list(1:num_partitions, :);

min_bmse_partitions = zeros(length(delta_phi_list), K+1);


for dd=1:length(delta_phi_list)

    delta_phi = delta_phi_list(dd);

    exp_factors = exp(-(0:NoAtoms).^2*delta_phi^2/2);

    m1 = toeplitz(exp_factors);

    c = (0:-1:-NoAtoms).*exp_factors;
    r = (0:NoAtoms).*exp_factors;
    exponential_matrix = toeplitz(c, r);

    m2 = 1i*delta_phi^2*exponential_matrix;

    bmse_min = 1000;
    min_ind = 1;
    
    
    for jj = 1:num_partitions
    
        coeffs = find_number_base_coeffs2(NoAtoms, partition_list(jj, :));
         
        bmse = find_min_bmse_ghz_no_dephasing(coeffs, delta_phi, m1, m2);

        if bmse < bmse_min
    
            bmse_min = bmse;
            min_ind = jj;
    
        end
    
    end
    
    %min
    delta_phi
    partition_list(min_ind, :)
    bmse_min
    % bmse_min/delta_phi^2

    bmse_min_list(dd) = bmse_min;
    min_bmse_partitions(dd, :) = partition_list(min_ind, :);

end

bmse_min_list = real(bmse_min_list)

close all
figure
loglog(delta_phi_list, sqrt(bmse_min_list).'./delta_phi_list)
xlabel('\delta\phi')
ylabel('\Delta\phi/\delta\phi')



