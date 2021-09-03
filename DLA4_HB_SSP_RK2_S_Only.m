function [S] = DLA4_HB_SSP_RK2_S_Only(x,v,k,C1,S,D1,dt,Awave,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic
%Convert BCs
if BC.use
    NN = size(BC.cell1,1);
end

S_ = S;

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

S = S1;
end

function CLuD = applyAonS(S,C,D,Awave)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + (C'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D);
end

end
