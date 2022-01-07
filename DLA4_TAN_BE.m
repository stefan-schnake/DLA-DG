function [C,S,D] = DLA4_TAN_BE(x,v,k,C,S,D,dt,Awave,BC)
%% DLA in hierarchical basis
%
% Update is U^{n+1} = U^n + T_{U^n}F(U^n)
%
% where T is the projection onto the tangent space of U^n

tic

obj = @(u) zfunc(u,C,S,D,dt,Awave,BC);
u0 = reshape(C*S*D',[],1);

u1 = fsolve(obj,0*u0);
U1 = reshape(u1,size(C,1),[]);
SS = svd(U1);

FU = kron(Awave{1,2},Awave{1,1});
for l=2:4
    FU = FU + kron(Awave{l,2},Awave{l,1});
end

u2 = (speye(size(FU))+dt*FU)\u0;
U2 = reshape(u2,size(C,1),[]);



toc
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

function z = zfunc(u,C,S,D,dt,Awave,BC)
    m = size(C,1); n = size(D,1);
    U = reshape(u,m,n);
    Z = U - C*S*D';
    [U,Sig,V] = svd(U);
    Sig = diag(Sig);
    rr = sum(Sig > 1e-12);
    Sig = diag(Sig(1:rr));
    U = U(:,1:rr);
    V = V(:,1:rr);
    Z = Z + dt*applyAonK(U*Sig,V,Awave)*V';
    Z = Z + dt*U*applyAonL(V*Sig',U,Awave)';
    Z = Z - dt*U*applyAonS(Sig,U,V,Awave)*V';
    z = reshape(Z,[],1);
end



