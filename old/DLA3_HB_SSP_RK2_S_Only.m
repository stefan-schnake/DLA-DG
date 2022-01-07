function [S] = DLA3_HB_SSP_RK2_S_Only(x,v,k,C1,S,D1,dt,Awave,BC,FMWT)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic
%Convert BCs
if BC.use
BCc{1} = convertVectoMat(x,v,k,BC.vec(:,1));
BCc{3} = convertVectoMat(x,v,k,BC.vec(:,3));
end

S_ = S;

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C1,D1,Awave);
if BC.use
    F1 = F(S_) - dt*(FMWT'*C1)'*BCc{1}*(FMWT'*D1);
    S1 = (1/2)*S_ + (1/2)*(F(F1) - dt*(FMWT'*C1)'*BCc{3}*(FMWT'*D1));
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
