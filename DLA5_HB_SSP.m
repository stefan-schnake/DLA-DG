function [C,S,D] = DLA5_HB_SSP(x,v,k,C,S,D,dt,Awave,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

% This algorithm performs SSP-RK3 on the system
% U_t = P_H P_{T_U^0 M_r}F(U)
% where H = C_1ZD_1' and C_1,D_1 are created from the DLA


tic

if BC.use
    NN = size(BC.cell1,1);
end

%%%Update K = CS
K = C*S;

%%%SSP-RK3
F = @(g) g - dt*applyAonK(g,D,Awave);
F1 = F(K);
if BC.use
    for l=1:NN
        F1 = F1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
    end
end
F2 = (3/4)*K + (1/4)*F(F1);
if BC.use
    for l=1:NN
        F2 = F2 - 1/4*dt*BC.cell3{l,1}*BC.cell3{l,2}*(BC.cell3{l,3}'*D);
    end
end
K_ = (1/3)*K + (2/3)*F(F2);
if BC.use
    for l=1:NN
        K_ = K_ - 2/3*dt*BC.cell2{l,1}*BC.cell2{l,2}*(BC.cell2{l,3}'*D);
    end
end

%Run QR to get new K=C*S;
[C1,~] = qr(K_,0);
M = C1'*C;


%%%Update L = D*S'
L = D*S';

%%%SSP-RK3
F = @(g) g - dt*applyAonL(g,C,Awave);
F1 = F(L);
if BC.use
    for l=1:NN
        F1 = F1 - dt*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C);
    end
end
F2 = (3/4)*L + (1/4)*F(F1);
if BC.use
    for l=1:NN
        F2 = F2 - 1/4*dt*BC.cell3{l,3}*BC.cell3{l,2}'*(BC.cell3{l,1}'*C);
    end
end
L_ = (1/3)*L + (2/3)*F(F2);
if BC.use
    for l=1:NN
        L_ = L_ - 2/3*dt*BC.cell2{l,3}*BC.cell2{l,2}'*(BC.cell2{l,1}'*C);
    end
end

%Run QR to get new L=D*S';
[D1,~] = qr(L_,0);
N = D1'*D;



%%%Update S
S_ = M*S*N';

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C,D,C,D,Awave,C1,D1);
F1 = F(S_);
if BC.use
    for l=1:NN
        F1 = F1 - dt*M*(C' *BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
        F1 = F1 + dt*M*(C' *BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D )*N';
        F1 = F1 - dt*  (C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D )*N';
    end
end
F = @(g) g - dt*applyAonS(g,C1,D1,C,D,Awave,C1,D1); %Need to change because U^{1/3} lives in H
F2 = (3/4)*S_ + (1/4)*F(F1);
if BC.use
    for l=1:NN
        F2 = F2 - (1/4)*dt*M*(C' *BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D1);
        F2 = F2 + (1/4)*dt*M*(C' *BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D )*N';
        F2 = F2 - (1/4)*dt*  (C1'*BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D )*N';
    end
end
S1 = (1/3)*S_ + (2/3)*F(F2);
if BC.use
    for l=1:NN
        S1 = S1 - (2/3)*dt*M*(C' *BC.cell2{l,1})*BC.cell2{l,2}*(BC.cell2{l,3}'*D1);
        S1 = S1 + (2/3)*dt*M*(C' *BC.cell2{l,1})*BC.cell2{l,2}*(BC.cell2{l,3}'*D )*N';
        S1 = S1 - (2/3)*dt*  (C1'*BC.cell2{l,1})*BC.cell2{l,2}*(BC.cell2{l,3}'*D )*N';
    end
end

C = C1;
S = S1;
D = D1;

toc;
end


function LuD = applyAonK(K,D,Awave)
%%%Updates X via the DLA update

LuD = zeros(size(K));
for l=1:size(Awave,1)
    LuD = LuD + Awave{l,1}*K*(D'*Awave{l,2}'*D);
end

end

function CLuD = applyAonS(S,C,D,C0,D0,Awave,C1,D1)
%%%Updates S_t =  P_H P_{T_U^0} F(U) where U=C*S*D
M = C1'*C0;
N = D1'*D0;

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + M*(C0'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D1);
    CLuD = CLuD - M*(C0'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D0)*N';
    CLuD = CLuD +   (C1'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D0)*N';
end

end

function CLu = applyAonL(L,C,Awave)
%%%Updates V via the DLA update

CLu = zeros(size(L));
for l=1:size(Awave,1)
    CLu = CLu + Awave{l,2}*L*(C'*Awave{l,1}'*C);
end

end


