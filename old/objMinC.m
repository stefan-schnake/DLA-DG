function [fx,gradx] = objMinC(Cvec,D,U)
%Calculating \|U-CD'\|_F^2 and its gradient w.r.t C

%Get rank
r = size(D,2);
%Convert Cvec to C
C = reshape(Cvec,r,[])';
%Convert other way -- Cvec = reshape(C',[],1);

%Calculate error
err = U-C*D';
%and its norm
%fx = norm(err,'fro')^2;
fx = sum(err.*err,'all');

%Now gradient
gradxmat = -2*err*D;
gradx = reshape(gradxmat',[],1);

end

