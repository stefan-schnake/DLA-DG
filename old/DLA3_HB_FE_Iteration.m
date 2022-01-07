function [C,S,D] = DLA3_HB_FE_Iteration(x,v,k,C,S,D,dt,Awave,BC,FMWT)
%% DLA in hierarchical basis.
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%Convert BCs
BCc{1} = convertVectoMat(x,v,k,BC(:,1));
%BCc{2} = convertVectoMat(x,v,k,BC(:,2));
BCc{3} = convertVectoMat(x,v,k,BC(:,3));


C1 = C;
D1 = D;

for i=1:2000

%%%SSP-RK3
%F = @(g) g - dt*applyAonK(g,D1,Awave);
K1 = C*S*(D'*D1) - dt*applyAonK2(C*S,D,D1,Awave) - dt*FMWT*(BCc{1}*(FMWT'*D1));
%K_ = 0.5*K + 0.5*(F(K1)- dt*FMWT*(BCc{3}*(FMWT'*D1)));
K_ = K1;

%Run QR to get new K=C*S;
[C_,~] = qr(K_,0);


%%%SSP-RK3
%F = @(g) g - dt*applyAonL(g,C1,Awave);
%F = @(g) D*S'*(C'*C1) - dt*applyAonL2(g,C,C1,Awave);
L1 = D*S'*(C'*C1) - dt*applyAonL2(D*S',C,C1,Awave) - dt*FMWT*(BCc{1}'*(FMWT'*C1));
%L1 = F(L) - dt*FMWT*(BCc{1}'*(FMWT'*C1));
%L_ = 0.5*L + 0.5*(F(L1) - dt*FMWT*(BCc{3}'*(FMWT'*C1)));
L_ = L1;

%Run QR to get new L=D*S';
[D_,~] = qr(L_,0);

S = C_'*C1*S*D1'*D_;

C1 = C_;
D1 = D_;

end

C = C1;
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

function LuD = applyAonK2(K,D,D1,Awave)
%%%Updates X via the DLA update

LuD = zeros(size(K));
for l=1:size(Awave,1)
    LuD = LuD + Awave{l,1}*K*(D'*Awave{l,2}'*D1);
end

end

function CLu = applyAonL2(L,C,C1,Awave)
%%%Updates V via the DLA update

CLu = zeros(size(L));
for l=1:size(Awave,1)
    CLu = CLu + Awave{l,2}*L*(C'*Awave{l,1}'*C1);
end

end


