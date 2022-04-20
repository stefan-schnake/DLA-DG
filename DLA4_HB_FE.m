function [U] = DLA4_HB_FE(U,dt,Acell,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic

%%%Update K = CS
K = U.C*U.S;

if BC.use
    NN = size(BC.cell1,1);
end

%%%SSP-RK3
F = @(g) g - dt*applyAonK(g,U.D,Acell);
if BC.use
    K = F(K);
    for l=1:NN
        K = K - dt*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*U.D);
    end
else
    K = F(K);
end

%Run QR to get new K=C*S;
[C1,R1] = qr(K,0);
%[C1,R1] = qr([C K],0);
%K = K - C*(C'*K);
%[C_,sigma,~] = svd(K,'econ');
%rr = sum(diag(sigma) > 1e-14);
%C1 = [C C_(:,1:rr)];
M = C1'*U.C;


%%%Update L = D*S'
L = U.D*U.S';

%%%SSP-RK3
F = @(g) g - dt*applyAonL(g,U.C,Acell);
if BC.use
    L = F(L);
    for l=1:NN
        L = L - dt*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*U.C);
    end
else
    L = F(L);
end


%Run QR to get new L=D*S';
[D1,R2] = qr(L,0);
%[D1,R2] = qr([D L],0);
%L = L - D*(D'*L);
%[D_,sigma,~] = svd(L,'econ');
%rr = sum(diag(sigma) > 1e-14);
%D1 = [D D_(:,1:rr)];
N = D1'*U.D;

%%%Update S
S_ = M*U.S*N';

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C1,D1,Acell);
if BC.use
    S1 = F(S_);
    for l=1:NN
        S1 = S1 - dt*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
    end
else
    S1 = F(S_);
end


U.C = C1;
U.S = S1;
U.D = D1;

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

