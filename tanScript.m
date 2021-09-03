
N = 64;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 2;

r = 40;

r_cut = 1;

test = 1;

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%CFL = 1/N^2;
%dt = 0.8*CFL;
%CFL = (5/4)*(1/73536);
%dt = 0.9*CFL;
dt = pi/(8*ceil(1/CFL));
%dt = 0.05;
%dt = 0.05;

%T = .5;
T = pi;
%T = 3.5;

%%%-------------------------------------------

%L = buildLDGMatrixAlt(x,v,k);
%LDG = L'*L;
%Jmp = buildJumpMatrix(x,v,k,1);
%Adv = buildAdvectionMatrix(x,v,k,true);
%Adv_back = buildBackAdvectionMatrix(x,v,k,true);
FMWT = OperatorTwoScale_wavelet2(k+1,log2(N));

%Acell = buildAdvectionMatrixWithBlocks(x,v,k);
Acell = buildAdvectionMatrixWithBlocks2(x,v,k);
Awave = cell(size(Acell));
for i=1:numel(Acell)
    Awave{i} = FMWT*Acell{i}*FMWT';
end

plotMat = @(U) plotVec(x,v,k,convertMattoVec(x,v,k,U.mat));

init = @(x,y)  (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
%init = @(x,y) x.^2.*y + x;
u0 = buildNonSeparableSource(x,v,k,init);
U0.mat = convertVectoMat(x,v,k,u0);
U0 = LRSVD(U0);

C = U0.C;
D = U0.D;
% r = size(C,2);
% A = Acell{1,1};
% kry = cell(2,2);
% 
% V = C;
% 
% W = A*V;
% H = V'*W;
% W = W - V*H;
% [Q,R] = qr(W,0);
% VV = [C Q];
% kry{1,1} = [H;R];
% 
% [QQ,~] = qr(A*C - C*kry{1,1}(1:r,:),0);
% 
% kry{2,1} = [];
% W = A'*V;
% H = V'*W;
% kry{2,1} = [kry{2,1};H];
% W = W - V*H;
% %W = W - Q*(Q'*W);
% H = QQ'*W;
% kry{2,1} = [kry{2,1};H];
% W = W - QQ*H;
% [svdU,svdS,svdV] = svd(W,0);
% rr = sum(diag(svdS) > 1e-13);
% Q = svdU(:,1:rr);
% R = svdS(1:rr,1:rr)*svdV(:,1:rr)';
% VV = [VV Q(:,1)];
% kry{2,1} = [kry{2,1};R(1,:)];
% kry{1,1} = [kry{1,1};0*R(1,:)];

%[VV,RR] = qr(VV,0);
%kry{1,1} = RR*kry{1,1};
%kry{2,1} = RR*kry{2,1};

[CC,nodesC,kryC] = krylovData(Acell(:,1),C);

% for i=1:size(Acell,1)
%     disp(norm(Acell{i,1}*C-CC*kryC{i,1},'fro'))
% end
% disp(norm(CC'*CC - eye(size(CC,2))))


[DD,nodesD,kryD] = krylovData(Acell(:,2),D);
% for i=1:size(Acell,1)
%     disp(norm(Acell{i,2}*D-DD*kryD{i,1},'fro'))
% end
% 


FU = applyMatA(U0.mat,Acell);
FU2 = zeros(size(kryC{1},1),size(kryD{1},1));
for l=1:size(Acell,1)
    FU2 = FU2 + kryC{l}*U0.S*kryD{l}';
end

disp(norm(FU,'fro')-norm(FU2,'fro'));



% FU0.mat = applyMatA(U0.mat,Acell2);
% FU0 = LRSVD(FU0);
% TFU0.mat = calcTanProj(FU0.mat,U0.C,U0.D);
% TFU0 = LRSVD(TFU0);
% 
% [C1,S1,D1] = KLFP(U0.C,U0.S,U0.D,0.01,Acell);
% 
% U1.mat = U0.mat -0.01*TFU0.mat;
% U1 = LRSVD(U1);

%U1.mat = U0.mat -0.01*FU0.mat;
%U1 = LRSVD(U1);

%FU1.mat = applyMatA(U1.mat,Acell2);
%FU1 = LRSVD(FU1);

%TFU1.mat = calcTanProj(FU1.mat,U1.C,U1.D);
%TFU1 = LRSVD(TFU1);

function [VV,nodes,kry] = krylovData(Acell,C)

NN = size(Acell,1);
r = size(C,2);
kry = cell(NN,1);
VV = C;
nodes = zeros(NN,2);
    
for i=1:NN
    A = Acell{i};
    W = A*C;
    H = C'*W;
    kry{i,1} = H;
    W = W - C*H;
    for j=1:i-1 %Project off previous vector
%         V = A*C;
%         V = V - C*kry{j,1}(1:r,:);
%         for k=1:j-1
%             curC = Acell{k}*C - C*kry{k,1}(1:r,:);
%             for l=1:k-1
%                 curC = curC - 
%             end
%         end
%         
%         [V,~] = qr(V,0);
        V = VV(:,nodes(j,1):nodes(j,2));
        H = V'*W;
        kry{i,1} = [kry{i,1};H];
        W = W - V*H;
    end
    [svdU,svdS,svdV] = svd(W,0);
    rr = sum(diag(svdS) > 1e-10);
    Q = svdU(:,1:rr);
    R = svdS(1:rr,1:rr)*svdV(:,1:rr)';
    kry{i,1} = [kry{i,1};R];
    VV = [VV Q];
    if i == 1
        nodes(i,1) = r+1;
        nodes(i,2) = nodes(i,1)+rr-1;
    else
        nodes(i,1) = nodes(i-1,2)+1;
        nodes(i,2) = nodes(i,1)+rr-1;
    end
    for j=1:i-1
        kry{j,1} = [kry{j,1};0*R];
    end
end

end

function LU = applyMatA(U,Awave)
    LU = 0*U;
    for l=1:size(Awave,1)
        LU = LU + Awave{l,1}*U*Awave{l,2}';
    end
end

function PU = calcTanProj(U,C,D)
    PU = (U*D)*D' - C*(C'*U*D)*D' + C*(C'*U);
end

function C = LRSVD(C)
    [U,S,V] = svd(C.mat);
    S = diag(S);
    r = sum(S > 1e-14);
    C.C = U(:,1:r);
    C.S = diag(S(1:r));
    C.D = V(:,1:r);
end

function [C1,S1,D1] = KLFP(C,S,D,dt,Acell)

    NN = size(Acell,1);
    
    C1 = C;
    D1 = D;

    for i=1:100
        K = C*S;
        L = S*D';
        for l=1:NN
            K = K - dt*Acell{l,1}*C*S*(D1'*Acell{l,2}'*D1);
            L = L - dt*(C1'*Acell{l,1}*C1)*S*D'*Acell{l,2}';
        end
        [C1,RK] = qr(K,0);
        [D1,RL] = qr(L',0);
        [RK,RS,RL] = svd(RK*RL');
        C1 = C1*RK;
        D1 = D1*RL;
    end
    M = C1'*C;
    N = D1'*D;
    
    %Evolve S
    S1 = M*S*N';
    for l=1:NN
        S1 = S1 - dt*(C1'*Acell{l,1}*C1)*(M*S*N')*(D1'*Acell{l,2}'*D1);
    end
end