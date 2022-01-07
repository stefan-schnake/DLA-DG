function [C,S,D] = DLA2_HB_SSP(x,v,k,C,S,D,dt,A,BC,FMWT)
%% DLA in hierarchical basis

%Convert BCs
%r = size(S,1);
BCc{1} = convertVectoMat(x,v,k,BC(:,1));
%[UB,SB,VB] = svd(BCc{1});
%BCc{1} = UB(:,1:r)*SB(1:r,1:r)*VB(:,1:r)';
BCc{2} = convertVectoMat(x,v,k,BC(:,2));
%[UB,SB,VB] = svd(BCc{2});
%BCc{2} = UB(:,1:r)*SB(1:r,1:r)*VB(:,1:r)';
BCc{3} = convertVectoMat(x,v,k,BC(:,3));
%[UB,SB,VB] = svd(BCc{3});
%BCc{3} = UB(:,1:r)*SB(1:r,1:r)*VB(:,1:r)';



%%%Update K = CS
K = C*S;
Kvec = reshape(K,[],1);

%%%SSP-RK3
F = @(g) g - dt*applyAonK(x,v,k,g,D,A,FMWT);
F1 = F(Kvec) - dt*reshape(FMWT*(BCc{1}*(FMWT'*D)),[],1);
F2 = (3/4)*Kvec + (1/4)*(F(F1)- dt*reshape(FMWT*(BCc{3}*(FMWT'*D)),[],1));
Kvec = (1/3)*Kvec + (2/3)*(F(F2)- dt*reshape(FMWT*(BCc{2}*(FMWT'*D)),[],1));

K = reshape(Kvec,size(K));
%Run QR to get new K=C*S;
[C1,~] = qr(K,0);
M = C1'*C;


%%%Update L = D*S'
L = D*S';
Lvec = reshape(L,[],1);

%%%SSP-RK3
F = @(g) g - dt*applyAonL(x,v,k,g,C,A,FMWT);
F1 = F(Lvec) - dt*reshape(FMWT*(BCc{1}'*(FMWT'*C)),[],1);
F2 = (3/4)*Lvec + (1/4)*(F(F1) - dt*reshape(FMWT*(BCc{3}'*(FMWT'*C)),[],1));
Lvec = (1/3)*Lvec + (2/3)*(F(F2) - dt*reshape(FMWT*(BCc{2}'*(FMWT'*C)),[],1));

L = reshape(Lvec,size(L));
%Run QR to get new L=D*S';
[D1,~] = qr(L,0);
N = D1'*D;



%%%Update S
Svec = reshape(M*S*N',[],1);

%%%SSP-RK3
F = @(g) g - dt*applyAonS(x,v,k,g,C1,D1,A,FMWT);
F1 = F(Svec) - dt*reshape((FMWT'*C1)'*BCc{1}*(FMWT'*D1),[],1);
F2 = (3/4)*Svec + (1/4)*(F(F1) - dt*reshape((FMWT'*C1)'*BCc{3}*(FMWT'*D1),[],1));
Svec = (1/3)*Svec + (2/3)*(F(F2) - dt*reshape((FMWT'*C1)'*BCc{2}*(FMWT'*D1),[],1));

C = C1;
S = reshape(Svec,size(S));
D = D1;
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After S update: ||u|| = %e\n',norm(u));


%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After V update: ||u|| = %e\n',norm(u));
end


function Kvec = applyAonK(x,v,k,Kvec,D,A,FMWT)
%%%Updates X via the DLA-LDG update

%Get size of K (same as size of D)
[Km,Kn] = size(D);

%Y = PSI*D is constant in t in this step
%Convert sliced K to a matrix
K = reshape(Kvec,Km,Kn);

%Determine u = K*D' and make a vector (also going from wavelet to
%realspace)
u = convertMattoVec(x,v,k,(FMWT'*K)*(FMWT'*D)');

%Apply A to u
u = A*u;

%Convert back to matrix and multiply to D
%Go back to wavelet space
K = FMWT*convertVectoMat(x,v,k,u)*(FMWT'*D);

%Slice into a vector
Kvec = reshape(K,[],1);

end

function Svec = applyAonS(x,v,k,Svec,C,D,A,FMWT)
%%%Updates X via the DLA-LDG update

%Get r (number of columns of D)
[~,r] = size(D);

%Y = PHI*D is constant in t in this step
%Convert sliced K to a matrix
S = reshape(Svec,r,r);

%Determine u = K*D and make a vector
u = convertMattoVec(x,v,k,(FMWT'*C)*S*(FMWT'*D)');

%Apply A to u
u = A*u;

%Convert back to matrix and multiply to D
S = (FMWT'*C)'*convertVectoMat(x,v,k,u)*(FMWT'*D);

%Slice into a vector
Svec = reshape(S,[],1);

end

function Lvec = applyAonL(x,v,k,Lvec,C,A,FMWT)
%%%Updates L=D*S' via the DLA-LDG update

%Get size of L (same as size of C)
[Km,Kn] = size(C);

%X = PHI*D is constant in t in this step
%Convert sliced K to a matrix
L = reshape(Lvec,Km,Kn);

%Determine u = K*D and make a vector
u = convertMattoVec(x,v,k,(FMWT'*C)*(FMWT'*L)');

%Apply A to u
u = A*u;

%Convert back to matrix and multiply to C
L = FMWT*(convertVectoMat(x,v,k,u)'*(FMWT'*C));

%Slice into a vector
Lvec = reshape(L,[],1);

end

