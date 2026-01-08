
function [coeff_list] = find_number_base_coeffs2(NoAtoms, partition)

    pwrs_of_2 = 2.^(length(partition)-1:-1:0);
    
    syms x
    
    expr = prod(expand((1.+x.^pwrs_of_2).^partition));

    coeff_list = sym2poly(expr);

    % coeff_list = [1];
    % 
    % for i=1:length(partition)
    % 
    %     a = sym2poly((1+x^(pwrs_of_2(i)))^partition(i));
    % 
    %     coeff_list = conv(coeff_list, a);
    % 
    % end

    coeff_list = sqrt(1.0*coeff_list)/sqrt(2^(sum(partition)));

    % coeff_list = zeros(sum(partition),NoAtoms);
    % 
    % ind = 0;
    % 
    % for i=1:length(partition)
    % 
    %     if partition(i) == 0
    %         continue;
    %     end
    % 
    %     coeff_list(ind+1:ind+partition(i), pwrs_of_2(i)+1) = 1;
    %     coeff_list(ind+1:ind+partition(i), 1) = 1;
    % 
    %     ind = ind + partition(i);
    % 
    % end
    % 
    % coeff_list = ifft(prod(fft(coeff_list,NoAtoms,2),1));
    % 
    % coeff_list = sqrt(coeff_list)/sqrt(2^(sum(partition)));


    %coeff_list = coeff_list/norm(coeff_list);

end