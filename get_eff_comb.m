function [E,V] = get_eff_comb(m, C, NumPorts, col)
% To draw the scatter plot of the efficient combination of the portfolio 
    % m - mean (expected return)
    % C - Covariances (risk)
    % NumPorts - number of random portfolio
    % col - number of assets
    
    % initialization
    w = zeros(NumPorts, col);
    E = ones(NumPorts, 1);
    V = zeros(NumPorts, 1);
    
    for i = 1: NumPorts
        w(i,:) = diff([0; sort(rand(col-1,1)); 1]);  %  portfolio weights 
        % Expected Return and Risk
        E(i) = w(i,:) * m;
        V(i) = w(i,:) * C * w(i,:)'; 
    end
    
end
