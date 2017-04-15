%Computational Finance: Assignment 1: Portfolio Optmization
clc; close all; clear
%% Q1.  % expected returns and corresponding covariances
m = [0.10  0.20  0.15]';
C = [0.005  -0.010  0.004; 
    -0.010   0.040 -0.002; 
     0.004  -0.002  0.023];
NumPorts = 100;


% Q1(a) Random portfolios
scatter_plot(m, C, NumPorts, length(m));

% Q1(b) Efficient portfolio
frontcon(m, C, NumPorts);

%Three two-asset portfolios
% pair wise 1, 2
C_12 = C(1:2,1:2);
m_12 = m(1:2);
scatter_plot(m_12, C_12, NumPorts, length(m_12));   %scatter plot
frontcon(m_12, C_12, NumPorts);       % Efficient frontier

% pair wise 1, 3
m_13 = m([1;3]);
C_13 = C([1,3; 7,9]);
scatter_plot(m_13, C_13, NumPorts, length(m_13));   %scatter plot
frontcon(m_12, C_13, NumPorts);       % Efficient frontier

% pair wise  2,3
m_23 = m(2:3);
C_23 = C(2:3,2:3);
scatter_plot(m_23, C_23, NumPorts, length(m_23));   %scatter plot
frontcon(m_23, C_23, NumPorts);       % Efficient frontier

% Q1(c) Why linprog is used in fig 6.11 [2]

% Q1(d) CVX optmization convex 
% CVX
    %run H:\MATLAB\cvx-w64\cvx\cvx_setup.m 
   % run H:\MATLAB\cvx-w64\cvx\cvx_startup.m


[V1, M1, PWts1] = CVX_NaiveMV(m, C, 25);
[V2, M2, PWts2] = NaiveMV(m, C, 25);

figure, clf,
plot(V1,M1, 'm', 'linewidth', 2); grid;
title('Mean Variance Portfolio NaiveMV ', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);
figure
plot(V2, M2, 'r', 'linewidth', 2); grid;
title('Mean Variance Portfolio CVX', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);








