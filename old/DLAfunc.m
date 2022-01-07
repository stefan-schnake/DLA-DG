function [C,S,D] = DLAfunc(x,v,k,C,S,D,dt,f1,f2)
%% DLA
%%%Update K = CS
K = C*S;
Kvec = reshape(K,[],1);

%%%SSP-RK3
F = @(g) g - dt*applyAonK(x,v,k,g,D,f1);
F1 = F(Kvec);
F2 = (3/4)*Kvec + (1/4)*F(F1);
Kvec = (1/3)*Kvec + (2/3)*F(F2);

K = reshape(Kvec,size(K));
%Run QR to get new K=C*S;
[C,S] = qr(K,0);
fprintf('- Norms:\n');
u = convertMattoVec(x,v,k,C*S*D');
fprintf('-- After X update: ||u|| = %e\n',norm(u));



%%%Update S
Svec = reshape(S,[],1);

%%%SSP-RK3
F = @(g) g + dt*applyAonS(x,v,k,g,C,D,f2);
F1 = F(Svec);
F2 = (3/4)*Svec + (1/4)*F(F1);
Svec = (1/3)*Svec + (2/3)*F(F2);

S = reshape(Svec,size(S));
u = convertMattoVec(x,v,k,C*S*D');
fprintf('-- After S update: ||u|| = %e\n',norm(u));



%%%Update L = D*S'
L = D*S';
Lvec = reshape(L,[],1);

%%%SSP-RK3
F = @(g) g - dt*applyAonL(x,v,k,g,C,f1);
F1 = F(Lvec);
F2 = (3/4)*Lvec + (1/4)*F(F1);
Lvec = (1/3)*Lvec + (2/3)*F(F2);

L = reshape(Lvec,size(L));
%Run QR to get new L=D*S';
[D,S] = qr(L,0);
S = S';
u = convertMattoVec(x,v,k,C*S*D');
fprintf('-- After V update: ||u|| = %e\n',norm(u));
end


function Kvec = applyAonK(x,v,k,Kvec,D,f)
%%%Updates X via the DLA-LDG update

%Get size of K (same as size of D)
[Km,Kn] = size(D);

%Y = PSI*D is constant in t in this step
%Convert sliced K to a matrix
K = reshape(Kvec,Km,Kn);

%Determine u = K*D' and make a vector
u = convertMattoVec(x,v,k,K*D');

%Apply A to u
u = f(u);

%Convert back to matrix and multiply to D
K = convertVectoMat(x,v,k,u)*D;

%Slice into a vector
Kvec = reshape(K,[],1);

end

function Svec = applyAonS(x,v,k,Svec,C,D,f)
%%%Updates X via the DLA-LDG update

%Get r (number of columns of D)
[~,r] = size(D);

%Y = PHI*D is constant in t in this step
%Convert sliced K to a matrix
S = reshape(Svec,r,r);

%Determine u = K*D and make a vector
u = convertMattoVec(x,v,k,C*S*D');

%Apply A to u
u = f(u);

%Convert back to matrix and multiply to D
S = C'*convertVectoMat(x,v,k,u)*D;

%Slice into a vector
Svec = reshape(S,[],1);

end

function Lvec = applyAonL(x,v,k,Lvec,C,f)
%%%Updates L=D*S' via the DLA-LDG update

%Get size of L (same as size of C)
[Km,Kn] = size(C);

%X = PHI*D is constant in t in this step
%Convert sliced K to a matrix
L = reshape(Lvec,Km,Kn);

%Determine u = K*D and make a vector
u = convertMattoVec(x,v,k,C*L');

%Apply A to u
u = f(u);

%Convert back to matrix and multiply to C
L = convertVectoMat(x,v,k,u)'*C;

%Slice into a vector
Lvec = reshape(L,[],1);

end

