function [C,S,D] = DLA(x,v,k,C,S,D,dt,A,A2,BC)
%% DLA

%Convert BCs
BCc = cell(3,1);
BCc{1} = convertVectoMat(x,v,k,BC(:,1));
BCc{2} = convertVectoMat(x,v,k,BC(:,2));
BCc{3} = convertVectoMat(x,v,k,BC(:,3));





%%%Update K = CS
K = C*S;
Kvec = reshape(K,[],1);
%Project BC onto X portion
%BC_K = reshape(convertVectoMat(x,v,k,BC)*D,[],1);
%A_K = eye(numel(Kvec));
%for j=1:numel(Kvec)
%    A_K(:,j) = LDGonXFunc(x,v,k,A_K(:,j),D,Adv);
%end
%Kvec = expm(-dt*A_K)*Kvec;

%%%BE
%[Kvec,flag,relres,iter] = gmres(@(g) g + dt*applyAonK(x,v,k,g,D,A),Kvec,20,1e-13);
%if flag
%    fprintf('-- Warning: gmres update for K did not succeed\n');
%end

%%%SSP-RK3
F = @(g) g - dt*applyAonK(x,v,k,g,D,A);
F1 = F(Kvec) - dt*reshape(BCc{1}*D,[],1);
F2 = (3/4)*Kvec + (1/4)*(F(F1) - dt*reshape(BCc{3}*D,[],1));
Kvec = (1/3)*Kvec + (2/3)*(F(F2) - dt*reshape(BCc{2}*D,[],1));


K = reshape(Kvec,size(K));
%Run QR to get new K=C*S;
[C,S] = qr(K,0);
%fprintf('- Norms:\n');
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After X update: ||u|| = %e\n',norm(u));

%%%Update S
Svec = reshape(S,[],1);
%BC_S = reshape(C'*convertVectoMat(x,v,k,BC)*D,[],1);
%BC_old_S = reshape(C'*convertVectoMat(x,v,k,BC_old)*D,[],1);
%BC_hlf_S = reshape(C'*convertVectoMat(x,v,k,BC_hlf)*D,[],1);

%%%BE
%[Svec,flag,relres,iter] = gmres(@(g)g - dt*applyAonS(x,v,k,g,C,D,A2),Svec,20,1e-13);
%if flag
%    fprintf('-- Warning: gmres update for S did not succeed\n');
%end

%%%SSP-RK3
F = @(g) g + dt*applyAonS(x,v,k,g,C,D,A2);
F1 = F(Svec) + dt*reshape(C'*BCc{1}*D,[],1);
F2 = (3/4)*Svec + (1/4)*(F(F1) + dt*reshape(C'*BCc{3}*D,[],1));
Svec = (1/3)*Svec + (2/3)*(F(F2) + dt*reshape(C'*BCc{2}*D,[],1));

%Svec = bicgstabl(@(g) g - dt*applyAonS(x,v,k,g,C,D,A),Svec,1e-13,numel(Svec));
%A_S = eye(numel(Svec));
%for j=1:numel(Svec)
%    A_S(:,j) = applyAonS(x,v,k,A_S(:,j),C,D,A2);
%end
%Simpson's approximation to non-homogeneous term
%  int_{t_0}^t exp(A*(t-s))BC(s) ds
%non_h = dt/6*( expm(A_S*dt)*BC_old_S + ...
%               expm(A_S*dt/2)*BC_hlf_S + ...
%               BC_S );
%Svec = expm(dt*A_S)*Svec;
S = reshape(Svec,size(S));
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After S update: ||u|| = %e\n',norm(u));

%%%Update L = D*S'
L = D*S';
Lvec = reshape(L,[],1);
%BC_L = reshape(C'*convertVectoMat(x,v,k,BC),[],1);
%A_L = eye(numel(Lvec));
%for j=1:numel(Lvec)
%    A_L(:,j) = LDGonYFunc(x,v,k,A_L(:,j),C,Adv);
%end
%Lvec = expm(-dt*A_L)*Lvec;

%%%BE
%[Lvec,flag,relres,iter] = gmres(@(g)g + dt*applyAonL(x,v,k,g,C,A),Lvec,20,1e-13);
%if flag
%    fprintf('-- Warning: gmres update for L did not succeed\n');
%end

%%%SSP-RK3
F = @(g) g - dt*applyAonL(x,v,k,g,C,A);
F1 = F(Lvec) - dt*reshape(BCc{1}'*C,[],1);
F2 = (3/4)*Lvec + (1/4)*(F(F1) - dt*reshape(BCc{3}'*C,[],1));
Lvec = (1/3)*Lvec + (2/3)*(F(F2) - dt*reshape(BCc{2}'*C,[],1));


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

