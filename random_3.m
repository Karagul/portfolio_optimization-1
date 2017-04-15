function [r] = random_3(N)
% Randomly select N portfolios
% N -> Number of portfolios
    r = round(1 + rand(N, 1)*29);
    for i=1:N
        if size(ftse100_30(i))<500
            r = random_3();
        end
    end
end