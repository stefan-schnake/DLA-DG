function [C,S,D] = DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic

M_lim = 1;

%%%Update K = CS
K = C*S;

if BC.use
    NN = size(BC.cell1,1);
end

%%%SSP-RK3
F = @(g) g - dt*applyAonK(g,D,Awave);
if BC.use
    F1 = F(K);
    for l=1:NN
        F1 = F1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
        %Apply limiter
    end
    F1 = limiterWrapper(x,k,F1,M_lim);
    K = (1/2)*K + (1/2)*F(F1);
    for l=1:NN
        K = K - 0.5*dt*BC.cell3{l,1}*BC.cell3{l,2}*(BC.cell3{l,3}'*D);
        
    end
    K = limiterWrapper(x,k,K,M_lim);
else
    F1 = F(K);
    F1 = limiterWrapper(x,k,F1,M_lim);
    K = (1/2)*K + (1/2)*F(F1);
    K = limiterWrapper(x,k,K,M_lim);
end

%Run QR to get new K=C*S;
[C1,R1] = qr(K,0);
M = C1'*C;


%%%Update L = D*S'
L = D*S';

%%%SSP-RK3
F = @(g) g - dt*applyAonL(g,C,Awave);
if BC.use
    F1 = F(L);
    for l=1:NN
        F1 = F1 - dt*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C);
    end
    F1 = limiterWrapper(v,k,F1,M_lim);
    L = (1/2)*L + (1/2)*F(F1);
    for l=1:NN
        L = L - 0.5*dt*BC.cell3{l,3}*BC.cell3{l,2}'*(BC.cell3{l,1}'*C);
    end
    L = limiterWrapper(v,k,L,M_lim);
else
    F1 = F(L);
    F1 = limiterWrapper(v,k,F1,M_lim);
    L = (1/2)*L + (1/2)*F(F1);
    L = limiterWrapper(v,k,L,M_lim);
end


%Run QR to get new L=D*S';
[D1,R2] = qr(L,0);
N = D1'*D;

%%%Update S
S_ = M*S*N';

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C1,D1,Awave);
if BC.use
    F1 = F(S_);
    for l=1:NN
        F1 = F1 - dt*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
    end
    S1 = (1/2)*S_ + (1/2)*F(F1);
    for l=1:NN
        S1 = S1 - 0.5*dt*(C1'*BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D1);
    end
else
    F1 = F(S_);
    S1 = (1/2)*S_ + (1/2)*F(F1);
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

function X = limiterWrapper(x,k,X,M)
    %Wrapper to limit each vector
    for i=1:size(X,2)
        X(:,i) = slopeLimiter(x,k,X(:,i),M,1);
    end
end

