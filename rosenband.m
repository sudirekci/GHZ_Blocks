
close all

noSamples = 1e3;

delta_phi = 1.4;
speeds = [4, 2, 1];

N_arr = 40:20:300;
error_var = zeros([1, length(N_arr)]);

for i = 1:length(N_arr)

    no_atoms = [40, 40, N_arr(i)];

    error_var(i) = rosenband_sim(delta_phi, speeds, no_atoms, noSamples);

end

n_arr = N_arr+40/4+40/2;

figure
hold on
plot(n_arr, error_var, LineWidth=3)
plot(n_arr, 1./n_arr, LineWidth=3)
set(gca, 'XScale', 'log')
set(gca, 'YScale', 'log')
ax = gca;
ax.FontSize = 20; 
ax.FontName = "times";
xlabel("Number of Atoms")
ylabel("MSE (\Delta\phi)^2")


function [error_variance] = rosenband_sim(delta_phi, speeds, no_atoms, noSamples)
    
    samples = rand([noSamples, 1])*4*pi-pi*2;
    
    %samples = normrnd(0, delta_phi, [noSamples, 1]);
    
    for l = 1:length(speeds)

        %num_atoms = N*speeds(l);
        num_atoms = no_atoms(l);

        if num_atoms == 0
            continue
        end
    
        r1 = binornd(round(num_atoms), (1+sin(samples./speeds(l)))/2)./round(num_atoms);
        %r2 = binornd(round(num_atoms/2), (1+sin(samples./speeds(l)))/2)./round(num_atoms/2);
        
        beta = asin(2*r1-1);

        %beta = angle(r1-1/2+1j*(r2-1/2));
    
        samples = samples - beta*speeds(l);
    
    end
    
    error_variance = sum(samples.^2)/noSamples;

end