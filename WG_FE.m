function [C,S,D] = WG_FE(x,v,k,C,S,D,dt,Awave,BC,tol)
%% Adaptive Low-Rank Method by Wei-Guo
%% Ignore BCs for now!!!

NN = size(Awave,1);
r = size(S,1);

C1 = zeros(size(C,1),r*(NN+1));
D1 = zeros(size(D,1),r*(NN+1));
S1 = zeros(r*(NN+1),r*(NN+1));

C1(:,1:r) = C;
D1(:,1:r) = D;
S1(1:r,1:r) = S;

for i=1:NN
    idx = r*i+1:r*(i+1);
    C1(:,idx) = Awave{i,1}*C;
    D1(:,idx) = Awave{i,2}*D;
    S1(idx,idx) = -dt*S;
end
if BC.use
    for i=1:size(BC.cell1)
        C1 = [C1 BC.cell1{i,1}];
        m = size(S1,1);
        n = size(BC.cell1{i,2},1);
        S1 = [S1 zeros(m,n);zeros(n,m) -dt*BC.cell1{i,2}];
        D1 = [D1 BC.cell1{i,3}];
    end
end

%Orthogonalize
[C1,RC] = qr(C1,0);
[D1,RD] = qr(D1,0);
[U_S,S1,V_S] = svd(RC*S1*RD');
sig = diag(S1);

%Want maximum tail <= tol
for rr=numel(sig)+1:-1:1
    if norm(sig(rr:numel(sig))) > tol
        break
    end
end
rr = min([rr,size(S1,1)]);

C = C1*U_S(:,1:rr);
D = D1*V_S(:,1:rr);
S = S1(1:rr,1:rr);





end




