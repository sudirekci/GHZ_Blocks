
deltaPhi = 0.7;

no_atoms_list = [140];

bmse_min_list = zeros(1, length(no_atoms_list));
min_bmse_partitions = zeros(length(no_atoms_list), 10);

tic

for dd=1:length(no_atoms_list)

    NoAtoms = no_atoms_list(dd);
    
    exp_factors = exp(-(0:NoAtoms).^2*deltaPhi^2/2);

    m1 = toeplitz(exp_factors);

    c = (0:-1:-NoAtoms).*exp_factors;
    r = (0:NoAtoms).*exp_factors;
    exponential_matrix = toeplitz(c, r);

    m2 = 1i*deltaPhi^2*exponential_matrix;

    not_found = true;
    k = -1;
    
    while not_found
        
        k = k + 1;
        
        if 2^k > NoAtoms
            
            not_found = false;
            
        end
        
    end
    
    K = k-1;
    
    partition_list = zeros(50000, K+1);
   
    
    [num_partitions, partition_list] = find_partitions(NoAtoms, K, partition_list, 1);
    
    num_partitions
    
    partition_list = partition_list(1:num_partitions, :);

    % partition_list = zeros(1,7)
    % partition_list(1,:) = [3,3,3,3,3,3,3]
    % num_partitions=1;
    
    
    bmse_min = 1000;
    min_ind = 1;
    
    
    for jj = 1:num_partitions
        
        %tic
    
        coeffs = find_number_base_coeffs2(NoAtoms, partition_list(jj, :));
        
        %toc
        
        %tic
         
        [bmse, ~, ~] = find_min_bmse_ghz_no_dephasing(coeffs, deltaPhi, m1, m2);
        
        %toc
        %partition_list(jj, :)
        %bmse
    
        if abs(bmse) < abs(bmse_min) && ~isnan(real(bmse))
    
            bmse_min = bmse;
            min_ind = jj;
            %disp(bmse_min);
    
        end
        
        if mod(jj, 10000) == 0
            disp(jj/num_partitions*100);
        end
    
    end
    
    %min
    NoAtoms
    bmse_min
    partition_list(min_ind, :)
    
    bmse_min_list(dd) = bmse_min;
    min_bmse_partitions(dd, end-(K):end) = partition_list(min_ind, :);

end

toc

bmse_min_list


