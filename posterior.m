function P = posterior(phi1, phi2, model, data)

P = ones(length(phi2), length(phi1));
wb_in = waitbar(0, 'Posterior computation.. ');
for i = 1:1:length(phi2)
    for j = 1:1:length(phi1)
        for k = 1:1:length(data)
            
            P(i, j) = P(i, j) * exp( -0.5 * ((data(k).vt ...
                    - vt(data(k).dv, phi1(j), phi2(i), model)) ...
                    / data(k).sigma)^2);
        end
    end
    waitbar(i/length(phi2), wb_in)
end
close(wb_in)