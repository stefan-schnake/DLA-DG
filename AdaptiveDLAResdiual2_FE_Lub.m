function [C,S,D] = AdaptiveDLAResdiual2_FE_Lub(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update
tic

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
%C0 = C; S0 = S; D0 = D;
%Update
[C,S,D] = DLA4_HB_FE_Lub(x,v,k,C,S,D,dt,Awave,BC);
[S_C,Sig,S_D] = svd(S);
dSig = diag(Sig);
tolflag = true;
for i=1:(2*r)
    normF = norm(dSig(i+1:end),2);
    if normF < tol
        tolflag = false;
        q = i;
        break
    end
end
if tolflag 
    q = 2*r;
end
r = q;
%r = sum(diag(Sig) > tol) + 1;
fprintf('-- Adaptive: New rank: %d, tol = %e\n',r,Sig(end,end));
C = C*S_C(:,1:r);
S = Sig(1:r,1:r);
D = D*S_D(:,1:r);

toc
end
