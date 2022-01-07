function [C_,D_] = createKrylovBasis(C,S,D,dt,Acell)
%Outputs Krylov basis for testing

FU = zeros(size(C,1),size(D,1));
for i=1:size(Acell,1)
    FU = FU - Acell{i,1}*C*S*D'*Acell{i,2}';
end

tol = 1e-12;
count = 1;
flag = 1;
C_1 = C;
D_1 = D;
C_ = C;
D_ = D;
while flag
    %K = (C*S )*(D'*D_1) - dt*applyAonK(C*S ,D,D_1,Acell);
    %L = (D*S')*(C'*C_1) - dt*applyAonL(D*S',C,C_1,Acell);
    K = (C*S )*(D'*D_1) - dt*applyAonK(C*S *(D'*D_1),D_1,D_1,Acell);
    L = (D*S')*(C'*C_1) - dt*applyAonL(D*S'*(C'*D_1),C_1,C_1,Acell);
    %Create matrix exponential
%     k = -applyAonK(C*S *(D'*D_1),D_1,D_1,Acell);
%     l = -applyAonL(D*S'*(C'*D_1),C_1,C_1,Acell);
%     K = C*S *(D'*D_1) + dt*k;
%     L = D*S'*(C'*C_1) + dt*l;
%     for j = 2:20
%         k = -applyAonK(k,D_1,D_1,Acell);
%         l = -applyAonL(l,C_1,C_1,Acell);
%         if 1/factorial(j)*dt^j*(norm(k,'fro')+norm(k,'fro')) < 1e-15
%             break
%         end
%         K = K + 1/factorial(j)*dt^j*k;
%         L = L + 1/factorial(j)*dt^j*l;
%     end
    K = K - C_*(C_'*K);
    L = L - D_*(D_'*L);
    [C_1,sigma,~] = svd(K,'econ');
    rr = sum(diag(sigma) > tol);
    C_1 = C_1(:,1:rr);
    [D_1,sigma,~] = svd(L,'econ');
    rr = sum(diag(sigma) > tol);
    D_1 = D_1(:,1:rr);
    
    C_ = [C_ C_1];
    D_ = [D_ D_1];
    
    R = dt*(FU - C_*(C_'*FU*D_)*D_');
    fprintf('iter = %02d: ||R||_F = %e. Added [%d,%d] vectors\n',count,norm(R,'fro'),size(C_1,2),size(D_1,2));
    if norm(R,'fro') < 1e-11 || count >= 10 || (size(C_1,2)+size(D_1,2)) == 0
        flag = 0;
    end
    count = count + 1;
end

end

function LuD = applyAonK(K,D,D_1,Awave)
%%%Updates X via the DLA update

LuD = zeros(size(D_1));
for l=1:size(Awave,1)
    LuD = LuD + Awave{l,1}*K*(D'*Awave{l,2}'*D_1);
end

end

function CLuD = applyAonS(S,C,D,C_1,D_1,Awave)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + (C_1'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D_1);
end

end

function CLu = applyAonL(L,C,C_1,Awave)
%%%Updates V via the DLA update

CLu = zeros(size(C_1));
for l=1:size(Awave,1)
    CLu = CLu + Awave{l,2}*L*(C'*Awave{l,1}'*C_1);
end

end

