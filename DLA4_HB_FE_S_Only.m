function [S] = DLA4_HB_FE_S_Only(x,v,k,C1,S,D1,dt,Awave,BC)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up
NN = size(BC,1);
%tic

%%%Update S
S_ = S;

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

S = S1;

end

function CLuD = applyAonS(S,C,D,Awave)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + (C'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D);
end

end

