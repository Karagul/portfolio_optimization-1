function [] = frontcon(m, C, NumPorts)
    % To draw the efficient portfolio frontier
    % m - mean (expected return)
    % C - Covariances (risk)
    % NumPorts - number of random portfolio
    p = Portfolio;
    p = setAssetMoments(p, m, C);
    p = setDefaultConstraints(p);
    figure; 
    plotFrontier(p, NumPorts);
end