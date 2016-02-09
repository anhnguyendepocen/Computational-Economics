function [Tp] = poly1D_cheb(n,X)

for j = 1:(n+1)
    if j == 1
        Tp(:,j) = ones(size(X));                % Chebychev poly 1
    elseif j == 2
        Tp(:,j) = X;                            % Chebychev poly 2
    else
        Tp(:,j) = 2*X.*Tp(:,j-1) - Tp(:,j-2);   % Chebychev poly n>2
    end
end
