clear; clc;  close all;
% Load financial data and organize it
ftse_matrix     = load('ftse.mat');
ftse_data       = ftse_matrix.data;
ftse_companies  = ftse_matrix.companies;
NumPorts = 30; 

% Select 3 random portfolio - index 1 is UKX with is the index FTSE
rand_3 = round(1 + rand(3, 1)*29);
stock1 = getfield(ftse_data, ftse_companies{rand_3(1)});
stock2 = getfield(ftse_data, ftse_companies{rand_3(2)});
stock3 = getfield(ftse_data, ftse_companies{rand_3(3)});

% Get the  AdjClose;
st1 = stock1.AdjClose;
st2 = stock2.AdjClose;
st3 = stock3.AdjClose;

% size and half size
min_size = getMin(st1, st2, st3);
half = round(min_size/2);

% Reverse the data
st1 = st1(end:-1:1);
st2 = st2(end:-1:1);
st3 = st3(end:-1:1);

% Build matrix of the 3 random stock
ftse = [st1(1:min_size) st2(1:min_size) st3(1:min_size)];

% Daily Return
ftse_return = price2ret(ftse);
[n, m] = size(ftse_return);

% Split into training and testing data
ftse_test  = ftse_return(1:half, :);
ftse_train = ftse_return(1+half:end, :);

% Expected Returns and Covariances
ERet_tr = mean(ftse_train)';
ECov_tr = cov(ftse_train);

%%
% efficient combinations
[E,V] = get_eff_comb(ERet_tr, ECov_tr, NumPorts, length(ERet_tr)); 
%figure;clf; scatter(E,V, 5, 'b'); title({'Efficient Combination Scatter ';' '}); xlabel('V'); ylabel('E'); grid;


%%   Using CVX to get the optimal weight - Markowitz
[T, N]  = size(ftse_train);
rho     = max(E);
R       = ftse_train;
tau     = 1;
mu      = ERet_tr(:);       % Ensuring its column vector

cvx_begin quiet
variable w(N);
    minimize( norm( rho * ones(T,1) - R*w));% + tau * norm(w,1) );
    subject to
        w' * ones(N,1) == 1;
        w' * mu        == rho;
        w              >= 0;
cvx_end
figure; clf, bar(w([3 1 2])); grid on;
title('Markowitz Optimal Portfolio Weight');
xlabel('Stocks');
ylabel('Weights');

%% % 1/N Naive portfolio Weights
naive_w = ones(N, 1) * (1/N);

figure; clf, bar(naive_w); grid on;
title('1/N Naive Portfolio Weight');
xlabel('Stocks');
ylabel('Weights');
%% Using the Estimation to design efficient frontier


       
% Q1(b) Efficient portfolio
p = get_frontcon(ERet_tr, ECov_tr);
figure; plotFrontier(p);  title('Efficient Frontier ');
figure; plotFrontier(p);  title('Efficient Frontier and Choosen Portfolios'); hold on;
fr_w = estimateFrontier(p);
[fr_prsk, fr_pret] = estimatePortMoments(p, fr_w);
scatter(fr_prsk, fr_pret, 'filled'); hold off;
text(fr_prsk(5), fr_pret(5),'\leftarrow My Chosen portfolio')
legend('Efficient Frontier','Optimal Portfolios');
hold off;

%% EFFICIENT PORTFOLIOS
figure; plotFrontier(p);  title('Efficient Frontier on Training Data'); hold on;

% Markowitz
mkwz_w = estimateFrontier(p);
[prsk, pret] = estimatePortMoments(p, mkwz_w);
mkwz_w = (max(mkwz_w'))';
[mkwz_prsk, mkwz_pret] = estimatePortMoments(p, w);
scatter(mkwz_prsk, mkwz_pret, 'm','filled');

% % 1/N Naive
% mkwz_w = estimateFrontier(p);
[naive_prsk, naive_pret] = estimatePortMoments(p, naive_w);
scatter(naive_prsk, naive_pret, 'g','filled');
legend('Efficient Frontier', 'Markowitz Portfolios', '1/N Naive Portfolio');

%text(prsk(5),pret(5),'\leftarrow Markowitz portfolio')
% scatter(sqrt(prsk), pret);



%% TESTING
% Expected Return and Risk of test data
ERet_ts = mean(ftse_test)';
ECov_ts = cov(ftse_test);

% 
% Efficient portfolio
pp = get_frontcon(ERet_ts, ECov_ts);
figure; plotFrontier(pp);  title('Efficient Frontier of test data'); hold on;
ts_w = estimateFrontier(pp);
[prsk_ts, pret_ts] = estimatePortMoments(pp, ts_w);
% scatter(fr_prsk, fr_pret, 'filled'); hold off;
% text(fr_prsk(5), fr_pret(5),'\leftarrow My Chosen portfolio')
% legend('Efficient Frontier','Optimal Portfolios');
% hold off;
%% Constructing the Portfolio Risk and Return
% Markowitz
mkwz_pret_ts  = mkwz_w' * ERet_ts;
mkwz_prsk_ts  = sqrt(mkwz_w' * ECov_ts * w);
hold on; scatter(mkwz_prsk_ts, mkwz_pret_ts, 'm', 'filled');


% 1/N portfolio
naive_pret_ts  = naive_w' * ERet_ts;
naive_prsk_ts  = sqrt(naive_w' * ECov_ts * naive_w);
hold on; scatter(naive_prsk_ts, naive_pret_ts, 'g', 'filled');

legend('Efficient Frontier','Markowitz Portfolio', '1/N Naive portfolio');


%% SHARPE RATIO
% Sharpe ratio

% figure; plotFrontier(p);  title('Efficient Frontier with Maximum Sharpe Ratio');
% % Markowitz
% mkwz_sp = setInitPort(p, 0);
% mkwz_sw = estimateMaxSharpeRatio(mkwz_sp);
% [mkwz_srsk, mkwz_sret] = estimatePortMoments(mkwz_sp, mkwz_sw);
% hold on; scatter(mkwz_srsk, mkwz_sret, 'filled');
% text(mkwz_srsk, mkwz_sret,'\leftarrow Markowitz Sharpe Ratio');
% 
% % Naive
% naive_sp = setInitPort(p, 0);
% naive_sw = estimateMaxSharpeRatio(naive_sp);
% [naive_srsk, naive_sret] = estimatePortMoments(naive_sp, naive_sw);
% hold on; scatter(naive_srsk, naive_sret,  'filled');
% text(naive_srsk, naive_sret,'\leftarrow Naive Sharpe Ratio');
% 
% legend('Efficient Frontier', 'Markowitz Sharpe', '1/N Naive Sharpe'); hold off; 
% 

% mkwz_w = estimateFrontier(p);
% 



%Maximum Portfolio sharpe ratio
sp = setInitPort(pp, 0);
sw = estimateMaxSharpeRatio(sp);
[srsk, sret] = estimatePortMoments(sp, sw);
% [prsk_ts, pret_ts] = estimatePortMoments(sp, sw);

% Stock risk and returns
% Markowitz
[mkwz_srsk, mkwz_sret] = estimatePortMoments(pp, w);
% 1/N Naive
[naive_srsk, naive_sret] = estimatePortMoments(pp, naive_w);


% Markowitz and Naive Sharpe ratio
mkwz_psr  = (mkwz_sret) ./ mkwz_srsk;
naive_psr = (naive_sret)./naive_srsk;
% Portfolio sharpe ratio
psr  = (pret_ts)./ prsk_ts;

figure;
subplot(2,1,1); plot(prsk_ts, pret_ts); hold on;
scatter(mkwz_prsk_ts, mkwz_pret_ts, 'm', 'filled');
scatter(naive_prsk_ts, naive_pret_ts, 'g', 'filled');
title('\bf Efficient Frontier');
xlabel('Portfolio Risk');
ylabel('Portfolio Return');
legend('Efficient Frontier', 'Markowitz', '1/N Naive'); grid; hold off; 


subplot(2,1,2); plot(prsk_ts, psr); hold on;
scatter(srsk, mkwz_psr, 'm', 'filled');
scatter(srsk, naive_psr, 'g', 'filled');
title('\bf Sharpe Ratio');
xlabel('Portfolio Risk');
ylabel('Sharpe Ratio');
legend('Sharpe Ratio', 'Markowitz', '1/N Naive'); grid; hold off; 
hold off;

% 
% 
% %% % 1/N Naive portfolio Weights
% w_equal = ones(N, 1) * (1/N);
% 
% [prsk_equal, pret_equal] = estimatePortMoments(p, w_equal);
% scatter(prsk_equal, pret_equal);
% text(prsk_equal,pret_equal,'\leftarrow Markowitz portfolio')
% % scatter(sqrt(prsk), pret);
% 
% sp_equal = setInitPort(p, 0);
% swgt_equal = estimateMaxSharpeRatio(sp_equal);
% [srsk_equal, sret_equal] = estimatePortMoments(sp_equal, swgt_equal);
% scatter(srsk_equal, sret_equal);
% text(srsk_equal ,sret_equal,'\leftarrow Markowitz Sharpe Ratio');
% hold off; 
% 
% 
% psratio_equal = (pret_equal) ./ prsk_equal;
% ssratio_equal = (sret_equal)./ srsk_equal;
% figure;
% subplot(2,1,1); plot(prsk_equal, pret_equal); hold on;
% scatter(srsk_equal, sret_equal, 'g', 'filled');
% title('\bf Efficient Frontier');
% xlabel('Portfolio Risk');
% ylabel('Portfolio Return');
% hold off;
% 
% subplot(2,1,2); plot(prsk_equal, psratio_equal); hold on;
% scatter(srsk_equal, ssratio_equal, 'g', 'filled');
% title('\bf Sharpe Ratio');
% xlabel('Portfolio Risk');
% ylabel('Portfolio Return');
% hold off;
% 
% 
% 
% 
% %% Constructing the return on testing data
% 
% %%
% % Sharpe Ratio 
% 
% sharpeRatio1 = R_test1/C_test1
% sharpeRatio2 = R_test2/C_test2
% 
% 
% % Certainity equivalent return
% ceq1 = R_test1 - (1/2)*C_test1^2
% ceq2 = R_test2 - (1/2)*C_test2^2
% 
% % Turn over
% 
