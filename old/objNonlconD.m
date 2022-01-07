function [c,ceq] = objNonlconD(Dvec,C,tol,maxabs)
%Nonlocal constraint - make all singular values greater than tol
%Get rank
r = size(C,2);
%Convert Cvec to C
D = reshape(Dvec,[],r);
%Convert other way -- Cvec = reshape(C',[],1);

ceq = [];
%sigma = svd(D);
%c = tol - sigma;
c = [];


[Q,~] = qr(D,0);
l1 = sum(abs(Q));
c = [c;l1(:)-maxabs];
end

