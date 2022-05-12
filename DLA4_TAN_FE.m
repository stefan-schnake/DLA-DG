function [C,S,D] = DLA4_TAN_FE(x,v,k,C,S,D,dt,Awave,BC)
%% DLA in hierarchical basis
%
% Update is U^{n+1} = U^n + T_{U^n}F(U^n)
%
% where T is the projection onto the tangent space of U^n

%Here Z is a proxy for F(U^n) but it is never calculated on its own
tic
r = size(S,1);

if BC.use
    NN = size(BC.cell1,1);
end

ZD = -applyAonK(C*S,D,Awave);
if BC.use
for l=1:NN
    ZD = ZD - BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
end
end
[C1,RC] = qr([C ZD],0);
%C1 = [C ZD];

ZtC = -applyAonL(D*S',C,Awave);
if BC.use
for l=1:NN
    ZtC = ZtC - BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C);
end
end
[D1,RD] = qr([D ZtC],0);
%D1 = [D ZtC];

S_ = -applyAonS(S,C,D,Awave);
if BC.use
for l=1:NN
    S_ = S_ - (C'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
end
end
%S1 = RC*[S-dt*S_,dt*eye(r);dt*eye(r),0*S]*RD';
%[SU,S1,SV] = svd(S1);

%C0 = C; S0 = S; D0 = D;
%C = C1; S = S1; D = D1;
%r = size(S1,1);
%r = size(S1);

%C = C1*SU(:,1:r);
%S = S1(1:r,1:r);
%D = D1*SV(:,1:r);

C = C1;
D = D1;
%S = [S-dt*S_,dt*eye(r);dt*eye(r),zeros(r)];
S = RC*[S-dt*S_,dt*eye(r);dt*eye(r),zeros(r)]*RD';

tt = toc;
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

