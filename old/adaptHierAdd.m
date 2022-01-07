function [C,S,D] = adaptHierAdd(C,S,D,m,G)
%Adds attempts to find r+m rank approximation to U = G by running BE DLA on
%U_t + U = G.

%Here G is a sum of low-rank factors by a full-rank update.

%Initialize C,S,D by adding m random vectors on it
for i=1:m
    [C,S,D] = adaptPlus1(C,S,D);
end
%e5 = zeros(size(C,1),1); e5(5) = 1;
%C = [C e5];
%D = [D e5];
%S = [S zeros(size(S,1),1)];
%S = [S;zeros(size(S,2),1)'];
%S(end,end) = 1;

NN = size(G,1);

dt = 1;
for j=1:50
    %disp(sum(abs(C)))
    %disp(sum(abs(D)))
    %fprintf('------\n');
    %Update K=CS
    %K_ = 1/(1+dt)*(C*S + dt*G1*D);
    K_ = (1-dt/2)*C*S;
    for l=1:NN
        K_ = K_ + dt*(G{l,1}*G{l,2})*(G{l,3}'*D);
    end
    K_ = 1/(1+dt/2)*K_;
    [C_,~] = qr(K_,0);
    %C_ = modGS(K_);
    M = C_'*C;
    
    %Update L=DS'
    %L_ = 1/(1+dt)*(D*S' + dt*G1'*C);
    %L_ = 1/(1+dt/2)*( (1-dt/2)*D*S' + dt*G1'*C);
    L_ = (1-dt/2)*D*S';
    for l=1:NN
        L_ = L_ + dt*(G{l,3}*G{l,2}')*(G{l,1}'*C);
    end
    L_ = 1/(1+dt/2)*L_;
    [D_,~] = qr(L_,0);
    %D_ = modGS(L_);
    N = D_'*D;
    
    %Update S
    S_ = M*S*N';
    %S_ = 1/(1+dt)*(S_ + dt*C_'*G1*D_);
    %S_ = 1/(1+dt/2)*( (1-dt/2)*S_ + dt*C_'*G1*D_);
    S_ = (1-dt/2)*S_;
    for l=1:NN
        S_ = S_ + dt*(C_'*G{l,1})*G{l,2}*(G{l,3}'*D_);
    end
    S_ = 1/(1+dt/2)*S_;
    fprintf('%5.4e\n',norm(C*S*D'-C_*S_*D_','fro'));
    err = norm(C*S*D'-C_*S_*D_','fro');
    if err < 1e-12
        break
    end
    if err > 1e2
        error('Bad')
    end
    
    C = C_;
    D = D_;
    S = S_;
    
end

end

function Z = modGS(A)
%Modified gram schmidt

Z = A;
[~,n] = size(A);
Z(:,1) = Z(:,1)/norm(Z(:,1),2);
for i=2:n
    for j=1:n-1
        Z(:,i) = Z(:,i) - (Z(:,j)'*Z(:,i))*Z(:,j);
    end
    Z(:,i) = Z(:,i)/norm(Z(:,i),2);
end

end


