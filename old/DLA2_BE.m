function [C,S,D] = DLA2_BE(x,v,k,C,S,D,dt,A,BC)
%% DLA

%Convert BCs
BCc = convertVectoMat(x,v,k,BC(:,1));

options = optimoptions(@fsolve,'Display','iter-detailed','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12);

%%%Update K = CS
K = C*S;
Kvec = reshape(K,[],1);
%%%BE
%[Kvec,flag,relres,iter] = gmres(@(g) g + dt*applyAonK(x,v,k,g,D,A),Kvec - dt*reshape(BCc*D,[],1),20,1e-13);
%[Kvec,flag] = pcg(@(g) g + dt*applyAonK(x,v,k,g,D,A),Kvec - dt*reshape(BCc*D,[],1),1e-13,5*numel(Kvec));
Kvec = fsolve(@(g) g + dt*applyAonK(x,v,k,g,D,A) - Kvec + dt*reshape(BCc*D,[],1),Kvec,options);

K = reshape(Kvec,size(K));
%Run QR to get new K=C*S;
[C1,~] = qr(K,0);
M = C1'*C;


%%%Update L = D*S'
L = D*S';
Lvec = reshape(L,[],1);

%%%BE
%[Lvec,flag,relres,iter] = gmres(@(g)g + dt*applyAonL(x,v,k,g,C,A),Lvec - dt*reshape(BCc'*C,[],1),20,1e-13);
%[Lvec,flag] = pcg(@(g)g + dt*applyAonL(x,v,k,g,C,A),Lvec - dt*reshape(BCc'*C,[],1),1e-13,5*numel(Lvec));
Lvec = fsolve(@(g) g + dt*applyAonL(x,v,k,g,C,A) - Lvec + dt*reshape(BCc'*C,[],1),Lvec,options);

L = reshape(Lvec,size(L));
%Run QR to get new L=D*S';
[D1,~] = qr(L,0);
N = D1'*D;



%%%Update S
Svec = reshape(M*S*N',[],1);

%%%BE
%[Svec,flag,relres,iter] = gmres(@(g)g + dt*applyAonS(x,v,k,g,C1,D1,A),Svec - dt*reshape(C1'*BCc*D1,[],1),numel(Svec),1e-13);
%[Svec,flag,relres,iter] = pcg(@(g)g + dt*applyAonS(x,v,k,g,C1,D1,A),Svec - dt*reshape(C1'*BCc*D1,[],1),1e-13);
Svec = fsolve(@(g)g + dt*applyAonS(x,v,k,g,C1,D1,A) - Svec + dt*reshape(C1'*BCc*D1,[],1),Svec,options);

C = C1;
S = reshape(Svec,size(S));
D = D1;
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After S update: ||u|| = %e\n',norm(u));


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

