% cd H:\MATLAB\cvx-w64\cvx\examples
% quickstart























% m = 20; n = 10; p = 4;
% A = randn(m,n); b = randn(m,1);
% C = randn(p,n); d = randn(p,1); e = rand;
% cvx_begin
%     variable x(n)
%     minimize( norm( A * x - b, 2 ) )
%     subject to
%         C * x == d
%         norm( x, Inf ) <= e
% cvx_end
%  figure(1), clf, bar(x); grid on
% 
% % clc,clear 
% % T = 150; N = 50;
% % R = randn(T, N);
% % rho = 0.02;
% % tau = 1;
% % mu = rand(N,1);
% % cvx_begin quiet
% % variable w(N)
% %     minimize( norm(rho*ones(T,1)-R*w) + tau*norm(w,1) )
% %     subject to
% %         w'*ones(N,1) == 1
% %         w'*mu == rho
% %         w > 0
% % cvx_end
% % figure(1), clf, bar(w); grid on