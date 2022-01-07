function [C,S,D] = DLA_CN(x,v,k,C,S,D,dt,A,A2,BC)
%% DLA

%%%Update K = CS
K = C*S;
Kvec = reshape(K,[],1);
%%%CN
[Kvec,flag,relres,iter] = gmres(@(g) g + 0.5*dt*applyAonK(x,v,k,g,D,A),Kvec - 0.5*dt*applyAonK(x,v,k,Kvec,D,A),20,1e-13);
if flag
    fprintf('-- Warning: gmres update for K did not succeed\n');
end

K = reshape(Kvec,size(K));
%Run QR to get new K=C*S;
[C,S] = qr(K,0);
%fprintf('- Norms:\n');
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After X update: ||u|| = %e\n',norm(u));

%%%Update S
Svec = reshape(S,[],1);

%%%CN
[Svec,flag,relres,iter] = gmres(@(g)g - 0.5*dt*applyAonS(x,v,k,g,C,D,A2),Svec + 0.5*dt*applyAonS(x,v,k,Svec,C,D,A2),20,1e-13);
if flag
    fprintf('-- Warning: gmres update for S did not succeed\n');
end

S = reshape(Svec,size(S));
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After S update: ||u|| = %e\n',norm(u));

%%%Update L = D*S'
L = D*S';
Lvec = reshape(L,[],1);

%%%CN
[Lvec,flag,relres,iter] = gmres(@(g)g + 0.5*dt*applyAonL(x,v,k,g,C,A),Lvec-0.5*dt*applyAonL(x,v,k,Lvec,C,A),20,1e-13);
if flag
    fprintf('-- Warning: gmres update for L did not succeed\n');
end

L = reshape(Lvec,size(L));
%Run QR to get new L=D*S';
[D,S] = qr(L,0);
S = S';
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After V update: ||u|| = %e\n',norm(u));
end


function Kvec = applyAonK(x,v,k,Kvec,D,A)
%%%Updates X via the DLA-LDG update

%Get size of K (same as size of D)
[Km,Kn] = size(D);

%Y = PSI*D is constant in t in this step
%Convert sliced K to a matrix
K = reshape(Kvec,Km,Kn);

%Determine u = K*D' and make a vector
u = convertMattoVec(x,v,k,K*D');

%Apply A to u
u = A*u;

%Convert back to matrix and multiply to D
K = convertVectoMat(x,v,k,u)*D;

%Slice into a vector
Kvec = reshape(K,[],1);

end

function Svec = applyAonS(x,v,k,Svec,C,D,A)
%%%Updates X via the DLA-LDG update

%Get r (number of columns of D)
[~,r] = size(D);

%Y = PHI*D is constant in t in this step
%Convert sliced K to a matrix
S = reshape(Svec,r,r);

%Determine u = K*D and make a vector
u = convertMattoVec(x,v,k,C*S*D');

%Apply A to u
u = A*u;

%Convert back to matrix and multiply to D
S = C'*convertVectoMat(x,v,k,u)*D;

%Slice into a vector
Svec = reshape(S,[],1);

end

function Lvec = applyAonL(x,v,k,Lvec,C,A)
%%%Updates L=D*S' via the DLA-LDG update

%Get size of L (same as size of C)
[Km,Kn] = size(C);

%X = PHI*D is constant in t in this step
%Convert sliced K to a matrix
L = reshape(Lvec,Km,Kn);

%Determine u = K*D and make a vector
u = convertMattoVec(x,v,k,C*L');

%Apply A to u
u = A*u;

%Convert back to matrix and multiply to C
L = convertVectoMat(x,v,k,u)'*C;

%Slice into a vector
Lvec = reshape(L,[],1);

end

