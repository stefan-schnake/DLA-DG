function [C,S,D,ZD,FS1,ZtC] = DLA4_TAN_AB2(x,v,k,FC0,FS0,FD0,C1,S1,D1,dt,Awave,BC)
%% Adams-Bashforth 2 in DLA
%
% Update is U^{n+2} = U^{n+1} + dt*(1.5*T_{U^{n+1}F(U^n+1) - 0.5*T_{U^n}F(U^n))
%
% where T is the projection onto the tangent space of U^n
%
% [C1,S1,D1] = U^{n+1}
%
% [FC0,FS0,FD0] = T_{U^n}F(U^n)
%

persistent step
if isempty(step)
    step = 0;
end

if step == 0
    %RUN FE
    [C,S,D,ZD,FS1,ZtC] = FE(x,v,k,C1,S1,D1,dt,Awave,BC);
    step = 1;
else
    [C,S,D,ZD,FS1,ZtC] = AB2(x,v,k,FC0,FS0,FD0,C1,S1,D1,dt,Awave,BC);
end

end


function [C,S,D,FC0,FS0,FD0] = FE(x,v,k,C,S,D,dt,Awave,BC)
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

ZtC = -applyAonL(D*S',C,Awave);
if BC.use
for l=1:NN
    ZtC = ZtC - BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C);
end
end
[D1,RD] = qr([D ZtC],0);

FS0 = applyAonS(S,C,D,Awave);
if BC.use
for l=1:NN
    FS0 = FS0 + (C'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
end
end
S1 = RC*[S + dt*FS0,dt*eye(r);dt*eye(r),0*S]*RD';
[SU,S1,SV] = svd(S1);

%C0 = C; S0 = S; D0 = D;
%C = C1; S = S1; D = D1;

FC0 = [C ZD];
FS0 = [FS0,eye(r);eye(r),zeros(r)];
FD0 = [D ZtC];

C = C1*SU(:,1:r);
S = S1(1:r,1:r);
D = D1*SV(:,1:r);

tt = toc;
end

function [C,S,D,FC1,FS1,FD1] = AB2(x,v,k,FC0,FS0,FD0,C1,S1,D1,dt,Awave,BC)
r = size(S1,1);

if BC.use
    NN = size(BC.cell1,1);
end

ZD = -applyAonK(C1*S1,D1,Awave);
if BC.use
for l=1:NN
    ZD = ZD - BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
end
end

ZtC = -applyAonL(D1*S1',C1,Awave);
if BC.use
for l=1:NN
    ZtC = ZtC - BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C1);
end
end

FS1 = applyAonS(S1,C1,D1,Awave);
if BC.use
for l=1:NN
    FS1 = FS1 + (C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
end
end
RHS = [ [S1+1.5*dt*FS1,1.5*dt*eye(r);1.5*dt*eye(r),zeros(r)], zeros(2*r); zeros(2*r), -0.5*dt*FS0 ];
[C,RC] = qr([C1 ZD  FC0],0);
[D,RD] = qr([D1 ZtC FD0],0);
[SU,S1,SV] = svd(RC*RHS*RD');

%C0 = C; S0 = S; D0 = D;
%C = C1; S = S1; D = D1;
C = C*SU(:,1:r);
S = S1(1:r,1:r);
D = D*SV(:,1:r);

FC1 = [C1 ZD];
FS1 = [FS1,eye(r);eye(r),zeros(r)];
FD1 = [D1 ZtC];

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

