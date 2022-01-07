function [C,S,D,bestfrac] = lrSVDApprox(C,S,D,tol)
%Compute SVD approximation under criterion: eta(k) < tol < eta(k+1)
%where
% eta(k) = sum(sigma(1:k))/sum(sigma(1:k+1))

[U,SS,V] = svd(S);
sigma = diag(SS);
r1 = numel(sigma);
summ = cumsum(sigma);
frac = summ(1:end-1)./summ(2:end);

r2 = sum(frac < tol) + 2;

if r2 <= r1
    bestfrac = 0;
    C = C*U(:,1:r2);
    D = D*V(:,1:r2);
    S = SS(1:r2,1:r2);
else
    bestfrac = frac(end);
end
end

