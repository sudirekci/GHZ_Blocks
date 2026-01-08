function [num_partitions, partition_list] = find_partitions(NoAtoms, K, partition_list, last_ind)
    
    if NoAtoms == 0
        
        partition_list(last_ind, end-K:end) = zeros(1, K+1);
        num_partitions = 1;
        
    end

    if K>0
        
        m = idivide(int16(NoAtoms), int16(2^K));
        
        last_ind_i = last_ind;

            for ii = 0:m

                [num, partition_list] = find_partitions(NoAtoms-ii*2^K, K-1, partition_list, last_ind);
                
                partition_list(last_ind:last_ind+num-1, end-K) = ii;
                
                last_ind = last_ind + num;

            end
            
            num_partitions = last_ind - last_ind_i;
            
    else
        
        partition_list(last_ind, end) = NoAtoms;
        num_partitions = 1;
        
    end

end