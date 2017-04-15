clear; clc;  close all;
% Load financial data and organize it
ftse_matrix     = load('ftse.mat');
ftse_data       = ftse_matrix.data;
ftse_companies  = ftse_matrix.companies;
NumPorts = 30; 

% Select 3 random portfolio
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
[N, M] = size(ftse_return);

% Split into training and testing data
ftse_test  = ftse_return(1:half, :);
ftse_train = ftse_return(1+half:end, :);


% Expected Returns and Covariances
ERet_tr = mean(ftse_train)';
ECov_tr = cov(ftse_train);

%% Using the Estimation to design efficient frontier

% Q1(a) Random portfolios
[E,V] = get_eff_comb(ERet_tr, ECov_tr, NumPorts, length(ERet_tr)); 
%figure;clf; scatter(E,V, 5, 'b'); title({'Efficient Combination Scatter ';' '}); xlabel('V'); ylabel('E'); grid;
       
% Q1(b) Efficient portfolio
p = get_frontcon(ERet_tr, ECov_tr);
figure; plotFrontier(p, NumPorts);  title('Efficient Frontier on Training set'); 

%%
%{
 Get weight and expected V and M (Niranja way)
[V1, M1, PWts1] = NaiveMV(ERet_tr, ECov_tr, NumPorts);
figure, clf,
subplot(1,1,1); plot(V1,M1, 'm', 'linewidth', 2); grid;
title('Mean Variance Portfolio NaiveMV ', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);

% or this way
%}
%%
% Portfolio weight
[port_w, port_buy, port_sell] = estimateFrontier(p);
w1 = port_w(:,1);

%%
% Using CVX method
ERet_tr = ERet_tr(:);       % makes sure it is a column vector
NAssets = length(ERet_tr);  % get number of assets
V0 = zeros(NAssets, 1);     % vector of lower bounds on weights
V1 = ones(1, NAssets);      % row vector of ones
    
cvx_begin quiet
    variable wMax(NAssets);
    maximize(ERet_tr' * wMax)
    subject to
        wMax' * ones(NAssets,1) == V1;
        wMax >= V0;
cvx_end

cvx_begin quiet
        variable wMin(NAssets);
        minimize(wMin' * ECov_tr * wMin);
        subject to
            wMin' * ones(NAssets,1) == V1;
            wMin >= V0;
cvx_end
    
%%
% 1/N portfolio Weights
w_equal = ones(M, 1) * (1/M);

%%
% Expected Return and Covariance of test set
ERet_ts = mean(ftse_test)';
ECov_ts = cov(ftse_test);

%%
% Estimate expected returns on test set
Return1 = w1' * ERet_ts;
Risk1 = w1' * ECov_ts * w1; 
% figure; plot(Risk1, Return1,'--o'); hold on;

% option CVX
Return2 = wMax' * ERet_ts;
Risk2 = (wMin' * ECov_ts * wMin);
% plot(Risk2, Return2,'--o');

% option 3: Equal weights
Return3 = w_equal' * ERet_ts;
Risk3 = w_equal' * ECov_ts * w_equal; 
% plot(Risk3, Return3,'--o');

