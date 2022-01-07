function [C,S,D] = DLA_UC_EXP(x,v,k,C,S,D,dt,Awave,BC)
%% DLA Unconventional Integrator with matrix exponential

%tic
tol = 1e-15;

%% Update K = CS
K = C*S;

%%%Matrix Exponential
k = -applyAonK(K,D,Awave);
K = K + dt*k;
for n=2:20
    k = -applyAonK(k,D,Awave);
    K = K + 1/factorial(n)*dt^n*k;
    if 1/factorial(n)*dt^n*norm(k,'fro') < tol
        break
    end
end
%Run QR to get new K=C*S;
[C1,R1] = qr(K,0);
M = C1'*C;


%% Update L = D*S'
L = D*S';

%%%Matrix Exponential
l = -applyAonL(L,C,Awave);
L = L + dt*l;
for n=2:20
    l = -applyAonL(l,C,Awave);
    L = L + 1/factorial(n)*dt^n*l;
    if 1/factorial(n)*dt^n*norm(l,'fro') < tol
        break
    end
end
%Run QR to get new L=D*S';
[D1,R2] = qr(L,0);
N = D1'*D;

%% Update S
S_ = M*S*N';

s = -applyAonS(S_,C1,D1,Awave);
S1 = S_ + dt*s;
for n=2:20
    s = -applyAonS(s,C1,D1,Awave);
    S1 = S1 + 1/factorial(n)*dt^n*s;
    if 1/factorial(n)*dt^n*norm(s,'fro') < tol
        break
    end
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

