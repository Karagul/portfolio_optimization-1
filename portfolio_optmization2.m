%{ Computational Finance: Assignment 1: Portfolio Optmization
clc; close all; clear
%% Q1.  % expected returns and corresponding covariances
m = [0.10  0.20  0.15]';
C = [0.005  -0.010  0.004; 
    -0.010   0.040 -0.002; 
     0.004  -0.002  0.023];
NumPorts = 100;

figure;
% Q1(a) Random portfolios
[E,V] = get_eff_comb(m, C, NumPorts, length(m)); 
subplot(4,2,1); scatter(V,E, 5, 'b'); title({'Efficient Combination Scatter ';' '}); xlabel('V'); ylabel('E'); grid;
       
% Q1(b) Efficient portfolio
p = get_frontcon(m, C);
subplot(4,2,2); plotFrontier(p, NumPorts);  title('Efficient Frontier');
    
%Three two-asset portfolios

% pair wise 1, 2
C_12 = C(1:2,1:2);
m_12 = m(1:2);
[E12,V12] = get_eff_comb(m_12, C_12, NumPorts, length(m_12)); 
subplot(4,2,3); scatter(E12,V12, 5, 'b'); title({'Efficient Combination Scatter (1 2) ';' '}); xlabel('V'); ylabel('E'); grid;
p12 = get_frontcon(m_12, C_12);
subplot(4,2,4); plotFrontier(p12, NumPorts); title('Efficient Frontier (1 2)');

% pair wise 1, 3
m_13 = m([1;3]);
C_13 = C([1,3; 7,9]);
[E13,V13] = get_eff_comb(m_13, C_13, NumPorts, length(m_13)); 
subplot(4,2,5); scatter(E13,V13, 5, 'b'); title({'Efficient Combination Scatter (1 3) ';' '}); xlabel('V'); ylabel('E'); grid;
p13 = get_frontcon(m_13, C_13);
subplot(4,2,6); plotFrontier(p13, NumPorts); title('Efficient Frontier (1 3)');

% pair wise  2,3
m_23 = m(2:3);
C_23 = C(2:3, 2:3);
[E23,V23] = get_eff_comb(m_23, C_23, NumPorts, length(m_23)); 
subplot(4,2,7); scatter(E23,V23, 5, 'b'); title({'Efficient Combination Scatter (2 3) ';' '}); xlabel('V'); ylabel('E'); grid;
p23 = get_frontcon(m_23, C_23);
subplot(4,2,8); plotFrontier(p23, NumPorts); title('Efficient Frontier (2 3)');

% Q1(c) Why linprog is used in fig 6.11 [2]



% Q1(d) CVX optmization convex 
% CVX
run C:\Users\NwoyeCID\Documents\MATLAB\ML\cvx-w64\cvx\cvx_setup.m 
%     run H:\MATLAB\cvx-w64\cvx\cvx_setup.m 
   % run H:\MATLAB\cvx-w64\cvx\cvx_startup.m

[V1, M1, PWts1] = CVX_NaiveMV(m, C, 25);
[V2, M2, PWts2] = NaiveMV(m, C, 25);

figure, clf,
subplot(1,2,1); plot(V1,M1, 'm', 'linewidth', 2); grid;
title('Mean Variance Portfolio NaiveMV ', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);

subplot(1,2,2); plot(V2, M2, 'r', 'linewidth', 2); grid;
title('Mean Variance Portfolio CVX', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);


% Comparison
figure; title('Efficient Frontier'); hold on;
plotFrontier(p, NumPorts);  
plotFrontier(p12, NumPorts);
plotFrontier(p13, NumPorts);
plotFrontier(p23, NumPorts);
legend('whole assets','[1 2] assets','[1 3] assets','[2 3] assets'); hold off;

figure;  title({'Efficient Combination Scatter';' '}); hold on;
scatter(V,E, 5, 'k');
scatter(V12,E12, 5, 'b');
scatter(V13,E13, 5, 'r');
scatter(V23,E23, 5, 'g');
legend('whole assets','[1 2] assets','[1 3] assets','[2 3] assets'); grid on; hold off;


%% FTSE100 data for the past 3 years



