function [fx,gradx] = objMinD(Dvec,C,U)
%Calculating \|U-CD'\|_F^2 and its gradient w.r.t D

%Get rank
r = size(C,2);
%Convert Cvec to C
D = reshape(Dvec,[],r);
%Convert other way -- Cvec = reshape(C',[],1);

%Calculate error
err = U-C*D';
%and its norm
%fx = norm(err,'fro')^2;
fx = sum(err.*err,'all');

%Now gradient
gradxmat = -2*C'*err;
gradx = reshape(gradxmat',[],1);

end

