function [C,S,D] = DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic

%%%Update K = CS
K = C*S;

if BC.use
    NN = size(BC.cell1,1);
end

%%%SSP-RK3
F = @(g) g - dt*applyAonK(g,D,Awave);
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
%[C1,R1] = qr([C K],0);
%K = K - C*(C'*K);
%[C_,sigma,~] = svd(K,'econ');
%rr = sum(diag(sigma) > 1e-14);
%C1 = [C C_(:,1:rr)];
M = C1'*C;


%%%Update L = D*S'
L = D*S';

%%%SSP-RK3
F = @(g) g - dt*applyAonL(g,C,Awave);
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
%[D1,R2] = qr([D L],0);
%L = L - D*(D'*L);
%[D_,sigma,~] = svd(L,'econ');
%rr = sum(diag(sigma) > 1e-14);
%D1 = [D D_(:,1:rr)];
N = D1'*D;

%%%Update S
S_ = M*S*N';

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C1,D1,Awave);
if BC.use
    S1 = F(S_);
    for l=1:NN
        S1 = S1 - dt*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
    end
else
    S1 = F(S_);
end


C = C1;
S = S1;
D = D1;

end


function LuD = applyAonK(K,D,Awave)
%%%Updates X via the DLA update

LuD = zeros(size(K));
for l=1:size(Awave,1)
    LuD = LuD + Awave{l,1}*K*(D'*Awave{l,2}'*D);
end

end

function CLuD = applyAonS(S,C,D,Awave)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + (C'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D);
end

end

function CLu = applyAonL(L,C,Awave)
%%%Updates V via the DLA update

CLu = zeros(size(L));
for l=1:size(Awave,1)
    CLu = CLu + Awave{l,2}*L*(C'*Awave{l,1}'*C);
end

end

