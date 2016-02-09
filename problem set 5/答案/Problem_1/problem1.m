clear;
close all;

sigma2 = 4;
sigma = sqrt(sigma2);

% Gauss-Quadrature I:
disp('------------------------------------------------------------');
disp('Gauss-Hermite-Integration:');
for n = [2, 3, 4, 5, 7],
    [x, w] = qnwnorm(n, 0, sigma2);
    y_1 = x.^4;
    y_2 = x.^6;
    y_3 = 1./(1+x.^2);
    F_1 = sum(w.*y_1);
    F_2 = sum(w.*y_2);
    F_3 = sum(w.*y_3);
    disp(['n = ', num2str(n), ' and y=E(x^4) is ', num2str(F_1)]); 
    disp(['n = ', num2str(n), ' and y=E(x^6) is ', num2str(F_2)]); 
    disp(['n = ', num2str(n), ' and y=E(1/(1+x^2)) is ', num2str(F_3)]); 
end;
disp('------------------------------------------------------------');
disp(' ');    

% Monte-Carlo-Integration:
disp('------------------------------------------------------------');
disp('Monte-Carlo-Integration:');
for n = [100, 1000, 10000, 50000, 100000, 1000000, 10000000],     % with 10000000 it looks partially accurate!!!
    % reset state
    randn('state', 0);
    
    % standard normal
    z = randn(n, 1);
    
    % transformation:
    x = sigma * z;
    
    y_1 = x.^4;
    y_2 = x.^6;
    y_3 = 1./(1+x.^2);
    F_1 = 1/n * sum(y_1);
    F_2 = 1/n * sum(y_2);
    F_3 = 1/n * sum(y_3);
    disp(['n = ', num2str(n), ' and y=E(x^4) is ', num2str(F_1)]); 
    disp(['n = ', num2str(n), ' and y=E(x^6) is ', num2str(F_2)]); 
    disp(['n = ', num2str(n), ' and y=E(1/(1+x^2)) is ', num2str(F_3)]); 
end;
disp('------------------------------------------------------------');
disp(' ');    