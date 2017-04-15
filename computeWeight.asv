% function w = computeWeight(R)
% %     w_co = inv(R' * R + (tau * Eret)) * (R' * I - rand * e)
%     
%     % Using CVX to get the optimal weight
%     [T, N]  = size(R);
%     rho     = 0.03;
%     tau     = 1;
%     mu      = mean(R);       % Ensuring its column vector
% 
%     cvx_begin quiet
%       variable w(N);
%         minimize( norm( rho * ones(T,1) - R*w) + tau * norm(w,1) );
%         subject to
%             w' * ones(N,1) == 1;
%             w' * mu        == rho;
%             w               > 0;
%     cvx_end
% end
% 
% function w = computeWeight(f, R)
%     [T, N]  = size(R);
%     % Using CVX to get the optimal weight    
%     cvx_begin quiet
%       variable w(N);
%         minimize( norm(f - R*w) );
%         subject to
%             w' * ones(N,1) == 1;
%             w              >= 0;
%     cvx_end
% end

function w = computeWeight(f, R)
    [T, N]  = size(R);
    ERet = mean(R);
    % Using CVX to get the optimal weight    
    cvx_begin quiet
      variable w(N);
        minimize( ERet * w );
        subject to
            w' * ones(N,1) == 1;
            w              >= 0;
    cvx_end
end