function [w_greedy, S_greedy] = greedySelection(f, R, Co, N)
% INPUTS: y <- FTSE 100 index return. 
%         R <- Return matrix
%         Co<- target number of combinations
%         N <- Total number of assets
% OUTPUT: w_greedy <- weight by the algorithm
%         S_greedy <- set of greedy combinations of the assets

% initialize
    S = {};
    w = {};
    k = 1;
    err = 10000;
    while k <= Co
        for j = 1:N
            if(~ismember(j, S))
                w_co = computeWeight([S,  R(j)]);
                w_co_vector = [w; w_co];
                index(j) = w_co;
            end
        end
        
        Ret_err = f * w_co_vectors;
        
%         % find the minimum
%         for i = 1: size(w_co_vector)            
%             Ret_err = f * w_co_vectors;
%             if(Ret_err < err)
%                 err = Ret_err;
%                 new_w =
        min_s = min(Ret_err);
        w = [w w_];
    end
end

