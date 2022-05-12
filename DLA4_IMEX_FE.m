function [C,S,D] = DLA4_IMEX_FE(x,v,k,C,S,D,dt,Acell,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic
r = size(S,1);

%%%Update K = CS
K = C*S;

if BC.use
    NN = size(BC.cell1,1);
end

K_BE = BE_K(K,D,dt,Acell);
[C_BE,~] = qr(K_BE,0);

%%%SSP-RK3
F = @(g) g - dt*applyAonK(g,D,Acell);
if BC.use
    K = F(K);
    for l=1:NN
        K = K - dt*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
    end
else
    K = F(K);
end

%Run QR to get new K=C*S;
[C1,R1] = qr(K,0);
M = C1'*C;


%%%Update L = D*S'
L = D*S';

%%%SSP-RK3
F = @(g) g - dt*applyAonL(g,C,Acell);
if BC.use
    L = F(L);
    for l=1:NN
        L = L - dt*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C);
    end
else
    L = F(L);
end


%Run QR to get new L=D*S';
[D1,R2] = qr(L,0);
N = D1'*D;

%%%Update S
S_ = M*S*N';

%Build stiffness matrix by kronecker products
A = zeros(r^2);
for i=1:size(Acell,1)
    A = A + kron(D1'*Acell{i,2}*D1,C1'*Acell{i,1}*C1);
end

%RHS
F = zeros(r);
if BC.use
for i=1:NN
    F = F - dt*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
end
end

%Matrix solve BE
S1 = reshape((speye(r^2)+dt*A)\reshape(S_+F,[],1),r,r);
% CN
%S_ = reshape(S_,[],1); S1 = (speye(r^2)+0.5*dt*A)\(S_ - 0.5*dt*(A*S_)); S1 = reshape(S1,r,r);
% Two-stage 3rd order DIRK
% assert(BC.use == 0);
% S_ = reshape(S_,[],1);
% rk_A = [0.5+sqrt(3)/6, 0; -sqrt(3)/3, 0.5+sqrt(3)/6];
% rk_b = [0.5 0.5];
% k1 = (speye(r^2)+dt*rk_A(1,1)*A)\(-A*S_);
% k2 = (speye(r^2)+dt*rk_A(2,2)*A)\(-A*S_ - rk_A(2,1)*dt*A*k1);
% S1 = reshape(S_ + rk_b(1)*dt*k1 + rk_b(2)*dt*k2,r,r);

C = C1;
S = S1;
D = D1;

end


function LuD = applyAonK(K,D,Acell)
%%%Updates X via the DLA update

LuD = zeros(size(K));
for l=1:size(Acell,1)
    LuD = LuD + Acell{l,1}*K*(D'*Acell{l,2}'*D);
end

end

function CLuD = applyAonS(S,C,D,Acell)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Acell,1)
    CLuD = CLuD + (C'*Acell{l,1}*C)*S*(D'*Acell{l,2}'*D);
end

end

function CLu = applyAonL(L,C,Acell)
%%%Updates V via the DLA update

CLu = zeros(size(L));
for l=1:size(Acell,1)
    CLu = CLu + Acell{l,2}*L*(C'*Acell{l,1}'*C);
end

end

function K = BE_K(K,D,dt,Acell)
    k = pcg(@(x) matFreeAK(x,D,dt,Acell),K(:),1e-7,numel(K));
    K = reshape(k,[],size(K,2));
end

function z = matFreeAK(k,D,dt,Acell)
    K = reshape(k,[],size(D,2));
    z = reshape(K + dt*applyAonK(K,D,Acell),[],1);   
end

