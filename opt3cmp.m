clear; clc;  close all;
% Load financial data and organize it
ftse_matrix     = load('ftse.mat');
ftse_data       = ftse_matrix.data;
ftse_companies  = ftse_matrix.companies;
NumPorts        = 30; 

% Select & Build the matrix of the 30 stocks
row_size = 753;
series30 = [];
for i = 1:NumPorts
    stock = getfield(ftse_data, ftse_companies{i});
    ftse30(:,i) = stock.Close(1:row_size);
    series30 = [series30 ftse30(:,i) ];
end

% Select the ftse index
ukx = getfield(ftse_data, ftse_companies{31});
ftse100 = ukx.Close(1:row_size);

%% Time series plots
% FTSE30
ts30 = fints(ftse_data.UKX.Date(1:row_size), series30);

%FTSE100
ts100 = fints(ftse_data.UKX.Date(1:row_size), ftse_data.UKX.Close(1:row_size));
figure; plot(ts100); hold on; plot(ts30); hold off;
xlabel('\bf Date', 'FontSize', 20); ylabel('\bf Market Close Price', 'FontSize', 20); title('\bf Financial Time Series of the Market Index (FTSE100 and 30 stocks)', 'FontSize', 25);
leg = ['FTSE100 - UKX', ftse_companies(1:30)];
legend(leg);


%%

% Reorder the data % Reverse the data according to earliest date
ftse30 = ftse30(end:-1:1, :);
ftse100 = ftse100(end:-1:1, :);
series30 = series30(end:-1:1, :);

% Calculate Daily Return
ftse30_return = price2ret(ftse30);
ftse100_return = price2ret(ftse100);

% Split into Training and Testing Data
half = round(row_size/2);
ftse100_test  = ftse100_return(1:half, :);
ftse100_train = ftse100_return(1+half:end, :);
ftse30_test  = ftse30_return(1:half, :);
ftse30_train = ftse30_return(1+half:end, :);

% Expected Returns and Covariances
ERet30 = mean(ftse30_train)';
ECov30 = cov(ftse30_train);
ERet100 = mean(ftse100_return)';
ECov100 = cov(ftse100_return);

p = get_frontcon(ERet30, ECov30);
fr_w = estimateFrontier(p);
figure;clf; plotFrontier(p, 30);  title('Efficient Frontier');

%% splitted Time series plots
% FTSE30
ts30_test = fints(ftse_data.UKX.Date(half:row_size), series30(half:row_size,:));
ts30_train = fints(ftse_data.UKX.Date(1:half), series30(1:half,:));

%FTSE100
ts100_train = fints(ftse_data.UKX.Date(1:half), ftse_data.UKX.Close(1:half));
ts100_test = fints(ftse_data.UKX.Date(half:row_size), ftse_data.UKX.Close(half:row_size));

figure; plot(ts100_train); hold on; plot(ts30_train); hold off;
xlabel('\bf Date', 'FontSize', 20); ylabel('\bf Market Close Price', 'FontSize', 20); title('\bf Financial Time Series of Training half(FTSE100 and 30 stocks)', 'FontSize', 25);
leg = ['FTSE100 - UKX', ftse_companies(1:30)];
legend(leg);

figure; plot(ts100_test); hold on; plot(ts30_test); hold off;
xlabel('\bf Date', 'FontSize', 20); ylabel('\bf Market Close Price', 'FontSize', 20); title('\bf Financial Time Series of Testing half(FTSE100 and 30 stocks)', 'FontSize', 25);
legend(leg);

%% Index Tracking - Greedy selection

y   = ftse100_train;
R   = ftse30_train;
My_greedy_weight  = zeros(NumPorts,1);
for i = 1:6
    v = 1/(2^i);
    
    % set the values of some spare portfolio
    for j = 1:NumPorts
        if(My_greedy_weight(j) >= v) 
            %wt(j) = v;
        end
    end
    
    % greedy selection
    for k=1:NumPorts
        if(My_greedy_weight(k) == 0)
            My_greedy_weight(k) = v;
            track(k) = norm(y - R*My_greedy_weight);
            My_greedy_weight(k) = 0;
        else
            track(k) = norm(y - R*My_greedy_weight);
        end
    end
    
    % pick the minimum
    mini = inf;
    t_index = 1;
    for p = 1:NumPorts
        if(track(p) < mini)
            mini = track(p);
            t_index = p;
        end
    end  
    My_greedy_weight(t_index) = v;
end
figure; clf, bar((My_greedy_weight)); grid on; %xlabel(ftse_companies(:))
title('Greedy Selected Weight and rest of zero-weights');
%set(gca,'XTickLabel',ftse_companies(1:30));
xlabel('Stocks');
ylabel('Weights');

%% Niranjan Greedy Algorithm
f   = ftse100_train;
R   = ftse30_train;
S = [];
Niranjan_greedy_weight = [];
Co = 6; 
min_err = 10000;
for k = 1: Co
    for j = 1:NumPorts
        if(~ismember(j, S))
            % Using CVX to get the optimal weight                              
            R_co = ([ R(:,S), R(:,j) ]);      % Combine new Stock to build R matrix          
            w_co = computeWeight(f, R_co);    % compute weight vector
            Ret_err = norm(f -  (R_co * w_co));
            if(Ret_err < min_err)   % Choose the stock with minimum error of approxiamtion
                opt_w = w_co;
                stock = j;
            end
        end
    end
    Niranjan_greedy_weight =  opt_w;
    S = [S, stock];
end

%build the weight vector
Niranjan_greedy_weight_30 = zeros(NumPorts,1);
for i=1:length(S)
    Niranjan_greedy_weight_30(S(i)) = Niranjan_greedy_weight(i);
end

figure; clf, bar(Niranjan_greedy_weight); grid on;
title('Sparse Portfolio Weight by Greedy Algorithm');
Niranjan_greedy_selected_stocks = ftse_companies(S);
set(gca,'XTickLabel',Niranjan_greedy_selected_stocks);
ylabel('Weights');
xlabel('Stocks');
%%
% Q1(b) Efficient portfolio
p30 = get_frontcon(ERet30, ECov30);
figure; plotFrontier(p30, NumPorts);  title('Efficient Frontier (30 stocks)');

p100 = get_frontcon(ERet100, ECov100);
figure; plotFrontier(p100, 100);  title('Efficient Frontier FTSE100 index'); 
[port_w100] = estimateFrontier(p100);
[prsk, pret] = estimatePortMoments(p100, port_w100);
hold on; scatter(prsk, pret); hold off;


hold off
% Portfolio weight
[port_w, port_buy, port_sell] = estimateFrontier(p30);
w1 = port_w(:,1);



%%   Using CVX to get the optimal weight
% rho     = 0.2786;

[T, N]  = size(ftse30_train);
R       = ftse30_train;
% tau     = [0.000, 0.055, 0.085, 0.105, 0.125, 0.145, 0.165, 0.185, 0.195, 0.205];    % Threshold = 1e-5
tau     = [0.03, 0.05, 0.07, 0.09, 0.11, 0.13, 0.15, 0.17, 0.19, 0.3];    % Threshold = 1e-5
% tau     = [0.04, 0.80, 1.00, 1.02, 1.05, 1.08, 1.12, 1.20, 1.30, 1.50]; % Threshold = 1e-10
%  tau     = [0.52, 0.82, 0.92, 0.96, 1, 1.02, 1.06, 1.12, 1.22, 1.32]; % Threshold = 1e-10
mu      = ERet30(:);       % Ensuring its column vector
y       = ftse100_train;
% legd = zeros(length(tau));
for i=1:length(tau)
    cvx_begin quiet
    variable w2(N);
        minimize( norm(y - R*w2) + (tau(i) * norm(w2,1)) );
%         subject to 
%             w2' * ones(N,1) == 1;
%             w2 >= 0;
    cvx_end
    
%     Non-zero coefficient that are not witched off by the regularizer
    nnzw(i) = sum(w2>1e-5);
    sparse_w(:,i) = w2;
    list_legend(i) = string(sprintf(['tau =',num2str(tau(i))]));
end
nnzw


figure;
bar(sparse_w); grid on; hold on;
legend(list_legend);
title(' Portfolio Weights (Tuning Regularizer (0.52 - 1.32)');
xlabel('Stocks');
ylabel('Weights');

figure;
plot(tau, nnzw); grid; hold;
scatter(tau(8), nnzw(8), 'r', 'filled');
text(tau(8), nnzw(8),'\leftarrow 6 nonzero weights')
title(' Regularizer Tuning on Training Set (Threshold = 1e-5)');
xlabel('Regularizer');
ylabel('Number of Non-zero Weights');


%% Lab 4  Using weight that gives 6 stocks to compare.
[T, N]  = size(ftse30_train);
R       = ftse30_train;
gamma   =  0.17;    % Threshold = 1e-5
y       = ftse100_train;
cvx_begin quiet
variable w2(N);
    minimize( norm(y - R*w2) + (gamma * norm(w2,1)) );
cvx_end
%     Non-zero coefficient that are not witched off by the regularizer
nnzw = sum(w2>1e-5);
regularized_sparse_weight = w2;

figure; clf, bar(regularized_sparse_weight); grid on;
title({'Sparse Index Tracking Weight on Training set','Using L1 regularization (Gamma = 0.3 and Threshold =1e-5) '});
xlabel('\bfStocks');
ylabel('\bfWeights');

%% Comparisons
Wgreedy = My_greedy_weight;
WgreedyN = Niranjan_greedy_weight_30;
Wsparse = regularized_sparse_weight;

% Evaluation on Training data
Rgreedy = Wgreedy' * ftse30_train';
RNgreedy = WgreedyN' * ftse30_train';
Rsparse = Wsparse' * ftse30_train';

% Tracking Error 
greedy_Tracking_error = norm(ftse100_test - Rgreedy);
greedy_Tracking_errorN = norm(ftse100_test - RNgreedy);
sparse_Tracking_error = norm(ftse100_test - Rsparse);

disp('Training Tracking errors:');
disp({'My Greedy = ',greedy_Tracking_error});
disp({'Niranjan Greedy = ',greedy_Tracking_errorN});
disp({'Sparse = ',sparse_Tracking_error});

% Evaluation on Test data
Rgreedy = Wgreedy' * ftse30_test';
RNgreedy = WgreedyN' * ftse30_test';
Rsparse = Wsparse' * ftse30_test';

% Tracking Error 
greedy_Tracking_error = norm(ftse100_test - Rgreedy);
greedy_Tracking_errorN = norm(ftse100_test - RNgreedy);
sparse_Tracking_error = norm(ftse100_test - Rsparse);
disp('Tesing Tracking errors:');
disp({'My Greedy = ',greedy_Tracking_error});
disp({'Niranjan Greedy = ',greedy_Tracking_errorN});
disp({'Sparse = ',sparse_Tracking_error});


%% Sharpe
% Expected Returns and Covariances
ERet30_ts = mean(ftse30_test)';
ECov30_ts = cov(ftse30_test);
ERet100 = mean(ftse100_return)';
ECov100 = cov(ftse100_return);

p_ts = get_frontcon(ERet30_ts, ECov30_ts);
ts_w = estimateFrontier(p_ts);
[prsk_ts, pret_ts] = estimatePortMoments(p_ts, ts_w);
figure;clf; plotFrontier(p_ts, 30);  title('Efficient Frontier');

%Maximum Portfolio sharpe ratio
sp = setInitPort(p_ts, 0);
sw = estimateMaxSharpeRatio(sp);
[srsk, sret] = estimatePortMoments(sp, sw);
% [prsk_ts, pret_ts] = estimatePortMoments(sp, sw);

% Stock risk and returns
% Greedy
[mkwz_srsk, mkwz_sret] = estimatePortMoments(p_ts, WgreedyN);
% 1/N Sparse
[naive_srsk, naive_sret] = estimatePortMoments(p_ts, Wsparse);


% Markowitz and Naive Sharpe ratio
mkwz_psr  = (mkwz_sret) ./ mkwz_srsk;
naive_psr = (naive_sret)./naive_srsk;
% Portfolio sharpe ratio
psr  = (pret_ts)./ prsk_ts;

% figure;
% subplot(2,1,1); plot(prsk_ts, pret_ts); hold on;
% scatter(mkwz_prsk_ts, mkwz_pret_ts, 'm', 'filled');
% scatter(naive_prsk_ts, naive_pret_ts, 'g', 'filled');
% title('\bf Efficient Frontier');
% xlabel('Portfolio Risk');
% ylabel('Portfolio Return');
% legend('Efficient Frontier', 'Greedy', 'Sparse'); grid; hold off; 


figure(20); plot(prsk_ts, psr); hold on;
scatter(srsk, mkwz_psr, 'm', 'filled');
scatter(srsk, naive_psr, 'g', 'filled');
title('\bf Sharpe Ratio');
xlabel('Portfolio Risk');
ylabel('Sharpe Ratio');
legend('Sharpe Ratio', 'Greedy', 'Sparse'); grid; hold off; 
hold off;

% 
% % My greedy algorithm method
% Return_My_greedy1 = My_greedy_weight' * ftse30_train';
% Return_My_greedy1a = My_greedy_weight' * ERet30;
% Risk_My_greedy  = sqrt(My_greedy_weight' * ECov30 * My_greedy_weight);
% 
% % Niranjan greedy algorithm
% Return_My_greedy2 = Niranjan_greedy_weight_30' * ftse30_return';
% Return_My_greedy2a = Niranjan_greedy_weight_30' * ERet30;
% 
% % Sparse Index Tracking algorithm
% Return_Sparse1 = regularized_sparse_weight' * ftse30_return';
% Return_Sparse1a = regularized_sparse_weight' * ERet30;
% 
% 
% %% Graphs
% %FTSE 100
% % p100 = get_frontcon(ERet100, ECov100);
% % w100 = estimateFrontier(p100);
% % [prsk100, pret100] = estimatePortMoments(p100, w100);
% % plot(prsk100, pret100); hold on;
% 
% % FTSE 30
% p30 = get_frontcon(ERet30, ECov30);
% w30 = estimateFrontier(p30);
% [prsk30, pret30] = estimatePortMoments(p30, w30);
% plot(prsk30, pret30); hold on;
% 
% % p = get_frontcon(ERet30, ECov30);
% % plotFrontier(p); grid;
% 
% % Greed 6
% [mg_prsk, mg_pret] = estimatePortMoments(p30, My_greedy_weight);
% scatter(mg_prsk, mg_pret, 'y','filled');
% % Niranjan
% [ng_prsk, ng_pret] = estimatePortMoments(p30, Niranjan_greedy_weight_30);
% scatter(ng_prsk, ng_pret, 'r','filled');
% % Sparse
% [sp_prsk, sp_pret] = estimatePortMoments(p30, regularized_sparse_weight);
% scatter(sp_prsk, sp_pret, 'g','filled');
% % 
% % p6g = get_frontcon(Return_My_greedy1a, Risk_My_greedy);
% % w6g = estimateFrontier(p6g);
% % [prsk6g, pret6g] = estimatePortMoments(p6g, w6g);
% % scatter(sqrt(Risk_My_greedy), Return_My_greedy1a); hold on;
% %scatter(fr_prsk, fr_pret, 'filled'); hold off;
% %text(fr_prsk(5), fr_pret(5),'\leftarrow My Chosen portfolio')
% 
% 
% title('\bf Efficient Frontier');
% xlabel('Portfolio Risk');
% ylabel('Portfolio Return');
% legend('FTSE30', 'My Greedy 6','Greedy 6', 'Sparse 6'); grid; hold off; 
% 
% 
% % Compare by Sharpe ratio

% Tracking
plot(ftse100_test(1:100))
hold on;
plot(Rsparse);
% hold on;
% plot(Rgreedy);
% plot(RNgreedy);
xlabel('\bf Date', 'FontSize', 12);
ylabel('\bf Returns', 'FontSize', 12); 
title('\bf Index Tracking by Greedy Selection', 'FontSize', 12);
legend('FTSE100 returns', 'Greedy Algorithm'); grid on;
