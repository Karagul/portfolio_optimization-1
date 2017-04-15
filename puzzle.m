m1 = [0.15 0.2 0.08 0.1]';
C1 = [  0.20 0.05 -0.010 0.0
        0.05 0.30  0.015 0.0
       -0.01 0.015 0.100 0.0
        0.00 0.000 0.000 0.0];
    
m2 = [0.15 0.2 0.08]';
C2 = [ 0.20 0.050 -0.01
       0.05 0.300 0.015
      -0.01 0.015 0.10];
  
[V1, M1, PWts1] = NaiveMV(m1, C1, 25);
[V2, M2, PWts2] = NaiveMV(m2, C2, 25);

figure, clf,
plot(V1,M1, 'b', V2, M2, 'r', 'linewidth', 2); grid;
title('Mean Variance Portfolio LinProg/QuadProg', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);



[V11, M11, PWts11] = CVX_NaiveMV(m1, C1, 25);
[V21, M21, PWts21] = CVX_NaiveMV(m2, C2, 25);

figure, clf,
plot(V11,M11, 'b', V21, M21, 'r', 'linewidth', 2); grid;
title('Mean Variance Portfolio CVX', 'FontSize', 13)
xlabel('Portfolio Risk', 'FontSize',11)
ylabel('Portfolio Return', 'FontSize', 11);



