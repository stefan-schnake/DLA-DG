function [C,S,D] = AdaptiveDLAWithInitv5(x,v,k,C,S,D,updateDLA,updatePDE,tol,i,C0,S0,D0,FMWT)
%Adaptive algorithm for DLA update

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
[C,S,D] = updateDLA(C,S,D);
%Compute SVD of S
[S_x,Sig,S_v] = svd(S);
Sig = diag(Sig);
fprintf('-- Adaptive: Smallest singular value is %e\n',Sig(end));
if Sig(end) < tol %We don't need to bump up r
    rnew = sum(Sig > tol)+1; %See if we need to decrease r
    Sig = diag(Sig);
    if rnew ~= r
        r = rnew; 
        fprintf('-- Adaptive: New r = %d\n',r);
        S = Sig(1:r,1:r);
        C = C*S_x(:,1:r);
        D = D*S_v(:,1:r);
    end
elseif r < size(C,1) %Add and element and try again
    if i <= 1
        C = C0(:,1:r+1);
        S = S0(1:r+1,1:r+1);
        D = D0(:,1:r+1);
        [C,S,D] = AdaptiveDLAWithInitv5(x,v,k,C,S,D,updateDLA,updatePDE,tol,i,C0,S0,D0,FMWT);
    else 
        %Throw out solution and advance by PDE - truncating
        u = convertMattoVec(x,v,k,(FMWT'*C_old)*S_old*(FMWT'*D_old)');
        u = updatePDE(u);
        [C,S,D] = svd(convertVectoMat(x,v,k,u));
        r = sum(diag(S) > tol) + 1;
        fprintf('-- Adaptive: New r = %d\n',r);
        S = S(1:r,1:r);
        C = FMWT*C(:,1:r);
        D = FMWT*D(:,1:r);
    end
end

end

