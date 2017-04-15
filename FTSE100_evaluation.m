clear; clc;  close all;
% Load financial data and organize it
ftse_matrix     = load('ftse.mat');
ftse_data       = ftse_matrix.data;
ftse_companies  = ftse_matrix.companies;

% Select 3 random portfolio
rand_3 = round(1 + rand(3, 1)*29);

% Create their time series
stock1 = getfield(ftse_data, ftse_companies{rand_3(1)});
ts1 = fints(stock1.Date, stock1.AdjClose);
% ts1 = fillts(ts1, 'n', ftse_data.UKX.Date);    % Fill the missing data

stock2 = getfield(ftse_data, ftse_companies{rand_3(2)});
ts2 = fints(stock2.Date, stock2.AdjClose);
% ts2 = fillts(ts2, 'n', ftse_data.UKX.Date);

stock3 = getfield(ftse_data, ftse_companies{rand_3(3)});
ts3 = fints(stock3.Date, stock3.AdjClose);
% ts3 = fillts(ts3, 'n', ftse_data.UKX.Date);

% Get their retuens
ts_return1 = (ts1 - lagts(ts1, 1));% ./ lagts(ts1, 1);
ts_return1 = fillts(ts_return1, 'zero', {'17-Feb-2014';'17-Feb-2017'});
ts_return2 = (ts2 - lagts(ts2, 1));% ./ lagts(ts2, 1);
ts_return2 = fillts(ts_return2, 'zero', {'17-Feb-2014';'17-Feb-2017'});
ts_return3 = (ts3 - lagts(ts3, 1));% ./ lagts(ts3, 1);
ts_return3 = fillts(ts_return3, 'zero', {'17-Feb-2014';'17-Feb-2017'});


min_size = getMin(ts_return1,ts_return2,ts_return3);

ftse_return = [fts2mat(ts_return1(2:min_size)) fts2mat(ts_return2(2:min_size)) fts2mat(ts_return3(2:min_size))];


% Split into training and test data
% size and half size
[N, M] = size(ftse_return);
half = round(min_size/2);

% half x 3 matrix
ftse_train = ftse_return(1:half,:);
ftse_test  = ftse_return(1+half:end,:);


% Expected Returns and Covariances
[ExpReturn, ExpCovariance] = ewstats(ftse_train, 1);
ExpReturn = mean(ftse_train);
ExpCovariance = cov(ftse_train);

%% Using the Estimation to design efficient frontier
NumPorts = 100;
ExpReturn = ExpReturn';
 
% Q1(a) Random portfolios
[E,V] = get_eff_comb(ExpReturn, ExpCovariance, NumPorts, length(ExpReturn)); 
figure;clf; scatter(E,V, 5, 'b'); title({'Efficient Combination Scatter ';' '}); xlabel('V'); ylabel('E'); grid;
       
% Q1(b) Efficient portfolio
p = get_frontcon(ExpReturn, ExpCovariance);
figure;clf; plotFrontier(p, NumPorts);  title('Efficient Frontier');
 
%% Get weight and Calculating Maximum Return and Minimu Variance

% Get the weight using built-in function
[port_w, port_buy, port_sell] = estimateFrontier(p);
w = port_w(:,1);
Return1 = w' * ftse_test';
Risk1 = w' * ExpCovariance * w; 

% Weights 1/N portfolio
equal_w = ones(M, 1) * (1/M);

% 

% Using CVX
ExpReturn = ExpReturn(:); % makes sure it is a column vector
NAssets = length(ExpReturn); % get number of assets
% vector of lower bounds on weights
V0 = zeros(NAssets, 1);
% row vector of ones
V1 = ones(1, NAssets);
    
cvx_begin quiet
    variable MaxReturnWeights(NAssets);
    maximize(ExpReturn' * MaxReturnWeights)
    subject to
        MaxReturnWeights' * ones(NAssets,1) == V1;
        MaxReturnWeights >= V0;
cvx_end

% Estimate expected return
MaxReturn = MaxReturnWeights' * ftse_test';
MaxR = w' * ftse_test';

cvx_begin quiet
        variable MinVarWeights(NAssets);
        minimize(MinVarWeights' * ExpCovariance * MinVarWeights);
        subject to
            MinVarWeights' * ones(NAssets,1) == V1;
            MinVarWeights >= V0;
cvx_end
    
MinVarReturn = MinVarWeights' * ftse_test';
MinVarStd = sqrt(MinVarWeights' * ExpCovariance * MinVarWeights);
%% Compare to 1/N portfolio






















% for i=1:3
%     stock = getfield(ftse_data, ftse_companies{rand_3(i)});
%     ts(:,i) = stock(1:754,7);
% end
% ts(:,i) = fints(stock.Date(1:754), stock.AdjClose(1:754));
% % Estimate expected return and covariances of the time series
% ts1 = fints(stock.Date(1:754), stocks(:,1));

