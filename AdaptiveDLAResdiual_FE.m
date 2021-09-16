function [C,S,D] = AdaptiveDLAResdiual_FE(x,v,k,C,S,D,dt,tol,Acell,BC)
%Adaptive algorithm for DLA update
tic

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C0 = C; S0 = S; D0 = D;
%Update using DLA
tic
[C,S,D] = DLA4_HB_FE(x,v,k,C,S,D,dt,Acell,BC);
toc
%Next we need to determine if residual is within tolarance
tic
[R] = compFEResidual(C0,S0,D0,C,S,D,dt,Acell,BC);
toc
normR = norm(R,'fro');
fprintf('-- Adaptive: Norm of residual is %e with tolerance of %e\n',normR,tol);
if normR <= tol
    fprintf('        ... looking to decrease rank\n');
    [S_U,S,S_V] = svd(S);
else
    fprintf('        ... looking to increase rank\n');
end


end

function Y = firstOrderResidual(C,S,D,Acell,dt,BC,X,trans,C1,D1)
%Compute action of first order residual
%R(U) = (I-P_H)U - dt*(L(U)-P_H L(P_H U) + (I-P_H)F) 
%where F = 0.5*(F(t_0)+F(t_1))
NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

S_ = (C1'*C)*S*(D'*D1);

if strcmp(trans,'notransp')
    %Y = zeros(size(X));
    
    Y = C*S*(D'*X);
    Y = Y - C1*S_*(D1'*X);
    
    % Y = Y - (LU)*X
    for l=1:NN
        Y = Y - dt*(Acell{l,1}*C)*S*(D'*Acell{l,2}'*X);
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - dt*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*X);
        end
    end
    %Projection
    %Y = Y + CC'*(LU)*DD'*X
    for l=1:NN
        Y = Y + dt*C1*(C1'*Acell{l,1}*C1)*S_*(D1'*Acell{l,2}'*D1)*(D1'*X);
    end
    if BC.use
        %Y = Y + C1*(C1'*F*D1)*(D1'*X);
        for l=1:MM
            Y = Y + dt*C1*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1)*(D1'*X);
        end
    end
else
    %Y = zeros(size(X));
    Y = D*S'*(C'*X);
    Y = Y - D1*S_'*(C1'*X);
    % Y = -(LU)'*X
    for l=1:NN
        Y = Y - dt*(Acell{l,2}*D)*S'*(C'*Acell{l,1}'*X);
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - dt*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*X);
        end
    end
    %Projection
    %Y = Y + DD'*(LU)*CC'*X
    for l=1:NN
        Y = Y + dt*D1*(D1'*Acell{l,2}*D1)*S_'*(C1'*Acell{l,1}'*C1)*(C1'*X);
    end
    if BC.use
        %Y = Y + D1*(D1'*F'*C1)*(C1'*X);
        for l=1:MM
            Y = Y + dt*D1*(D1'*BC.cell1{l,3})*BC.cell1{l,2}'*(BC.cell1{l,1}'*C1)*(C1'*X);
        end
    end
    
end
    

end

function FU = computeFullGrid(C0,S0,D0,Acell,dt,BC)
NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

FU = C0*S0*D0';
for l=1:NN
    FU = FU - dt*Acell{l,1}*C0*S0*D0'*Acell{l,2}';
end

if BC.use
    for l=1:MM
        FU = FU - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
    end
end
    
end

function CC = recoverKrylovBasis(C0,kry,nodes,Acell)
    n = nodes(end,2);
    NN = size(Acell,1);
    
    CC = zeros(size(C0,1),size(C0,2)+n);
    CC(:,1:size(C0,2)) = C0;
    
    for i=1:NN
        st = nodes(i,1);
        ed = nodes(i,2);
        CC(:,st:ed) = Acell{i}*C0 - C0*kry{i}(1:r,:);
        for j=1:i-1
            iter = nodes(j,1):nodes(j,2);
            CC(:,st:ed) = CC(:,st:ed) - CC(iter)*kry{i}(iter,:);
        end    
        CC(:,st:ed) = kry{i}(st:ed,:)\CC(:,st:ed);
    end
end

