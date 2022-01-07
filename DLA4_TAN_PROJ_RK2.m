function [C,S,D] = DLA4_TAN_PROJ_RK2(x,v,k,C,S,D,dt,Awave,BC)
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

%eta_1 = [C,S,D]


%% Compute FE timestep
ZD = -applyAonK(C*S,D,Awave);
if BC.use
for l=1:NN
    ZD = ZD - BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
end
end
[C_,RC_] = qr([C ZD],0);

ZtC = -applyAonL(D*S',C,Awave);
if BC.use
for l=1:NN
    ZtC = ZtC - BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*C);
end
end
[D_,RD_] = qr([D ZtC],0);

S_ = applyAonS(S,C,D,Awave);
if BC.use
for l=1:NN
    S_ = S_ + (C'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D);
end
end
S1 = RC_*[S+dt*S_,dt*eye(r);dt*eye(r),0*S]*RD_';
[SU,S1,SV] = svd(S1);

rr = 0;
dS = diag(S1);
for i=1:numel(dS)-1
    if norm(dS(i+1:end),2) < dt^2
        rr = i;
        break
    end
end
if rr == 0 
    rr = numel(dS);
end
%rr = numel(dS);
rr = r;

%FE update
C_12 = C_*SU(:,1:rr);
S_12 = S1(1:rr,1:rr);
D_12 = D_*SV(:,1:rr); %eta_2 = [C_12,S_12,D_12]

C = C_12;
D = D_12;
S = S_12;
return

Spk1_S = RC_*[S+0.5*dt*S_,0.5*dt*eye(r);0.5*dt*eye(r),0*S]*RD_'; %kappa_1 = [C_,k1_S,D_]
 %above is storing S+0.5*dt*kappa_1


%% Compute next stage
ZD = -applyAonK(C_12*S_12,D_12,Awave);
if BC.use
for l=1:NN
    ZD = ZD - BC.cell3{l,1}*BC.cell3{l,2}*(BC.cell3{l,3}'*D_12);
end
end
%[C__,RC__] = qr([C_12,ZD],0);

ZtC = -applyAonL(D_12*S_12',C_12,Awave);
if BC.use
for l=1:NN
    ZtC = ZtC - BC.cell3{l,3}*BC.cell3{l,2}'*(BC.cell3{l,1}'*C_12);
end
end
%[D__,RD__] = qr([D_12,ZtC],0);

S__ = applyAonS(S_12,C_12,D_12,Awave);
if BC.use
for l=1:NN
    S__ = S__ + (C_12'*BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D_12);
end
end

%k2_S = RC__*[S__,eye(r);eye(r),zeros(r)]*RD__'; %kappa_2 = [C__,k2_S,D__]

% [C1,RC] = qr([C_,C__],0);
% [D1,RD] = qr([D_,D__],0);
% S1 = RC*[Spk1_S        ,zeros(2*r,2*r);...
%          zeros(2*r,2*r),dt*0.5*k2_S    ]*RD';

[C1,RC] = qr([C_,ZD ],0);
[D1,RD] = qr([D_,ZtC],0);

S1 = RC*[Spk1_S+0.5*dt*SU(:,1:rr)*S__*SV(:,1:rr)',0.5*dt*SU(:,1:rr);0.5*dt*SV(:,1:rr)',zeros(rr) ]*RD';
     
[UU,S1,VV] = svd(S1);

rr = 0;
dS = diag(S1);
for i=1:numel(dS)-1
    if norm(dS(i+1:end),2) < dt^2
        rr = i;
        break
    end
end
if rr == 0 
    rr = numel(dS);
end
%rr = numel(dS);
rr = r;

S = S1(1:rr,1:rr);
C = C1*UU(:,1:rr);
D = D1*VV(:,1:rr);

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

