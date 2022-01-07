function [c,ceq] = objNonlconC(Cvec,D,tol,maxabs)
%Nonlocal constraint - make all singular values greater than tol
%Get rank
r = size(D,2);
%Convert Cvec to C
C = reshape(Cvec,r,[])';
%Convert other way -- Cvec = reshape(C',[],1);

ceq = [];
%sigma = svd(C);
%c = tol - sigma;
c = [];

[Q,~] = qr(C,0);
l1 = sum(abs(Q));
c = [c;l1(:)-maxabs];
end

