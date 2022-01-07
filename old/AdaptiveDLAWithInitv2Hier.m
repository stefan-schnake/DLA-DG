function [C,S,D] = AdaptiveDLAWithInitv2Hier(C,S,D,updateDLA,plusOneDLA,tol,i,C0,S0,D0,switchval)
%Adaptive algorithm for DLA update

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
if switchval 
    for j=1:15
        [C,S,D] = plusOneDLA(C,S,D);
    end
else
    [C,S,D] = updateDLA(C,S,D);
end
%[C,S,D] = DLA(x,v,k,C_old,S_old,D_old,dt,A,A2);
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
    if i > 1
        [C,S,D] = adaptPlus1(C_old,S_old,D_old);
        [C,S,D] = AdaptiveDLAWithInitv2Hier(C,S,D,updateDLA,plusOneDLA,tol,i,C0,S0,D0,true);
    else %Add next row/col from initial SVD and try again
        %Create soln at current time step
        C = C0(:,1:r+1);
        S = S0(1:r+1,1:r+1);
        D = D0(:,1:r+1);
        [C,S,D] = AdaptiveDLAWithInitv2Hier(C,S,D,updateDLA,plusOneDLA,tol,i,C0,S0,D0,false);
    end
end

end

