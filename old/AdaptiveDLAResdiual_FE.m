function [C,S,D] = AdaptiveDLAResdiual_FE(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
[C,S,D] = DLA3_HB_FE(x,v,k,C,S,D,dt,Awave,BC,FMWT);
%Compute SVD of S
[S_x,Sig,S_v] = svd(S);
Sig = diag(Sig);
fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);
if Sig(end) < tol && 0 %We don't need to bump up r
    rnew = sum(Sig > tol)+1; %See if we need to decrease r
    Sig = diag(Sig);
    if rnew ~= r
        r = rnew; 
        fprintf('-- Adaptive: Reducing r to %d\n',r);
        S = Sig(1:r,1:r);
        C = C*S_x(:,1:r);
        D = D*S_v(:,1:r);
    end
elseif r < size(C,1) %Add most relavant vectors from residual
    %%%First order residual
    fprintf('-- Adaptive: Adding first order residual basis\n');
    Chalf = C;
    Dhalf = D;
    C0 = C_old; D0 = D_old; S0 = S_old;
    %[Chalf,~] = sRRQR([C0 Chalf],1.01,'tol',1e-12);
    %[Dhalf,~] = sRRQR([D0 Dhalf],1.01,'tol',1e-12);
    %%%Compute first order residual basis vectors
    FF = FMWT*convertVectoMat(x,v,k,BC(:,1))*FMWT';
    R0 = firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',Chalf,Dhalf);
    %[C1,S1,D1] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,FF,x,trans,Chalf,Dhalf),[1,1]*size(C,1),15);
    %[C1,S1,D1] = svd(firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',Chalf,Dhalf));
    %Prune irrelavent vectors
    %rr = sum(diag(S1) > 1e-16);
    %fprintf('-- Adaptive: Found %d non-trivial vectors to add\n',rr);
    %C1 = C1(:,1:rr); S1 = S1(1:rr,1:rr); D1 = D1(:,1:rr);
    
    [C,S,D] = svd(C*S*D'+R0);
    Sig = diag(S);

    fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);
    
    if (Sig(end) < tol) && 0 %We don't need to bump up r
        rnew = sum(Sig > tol)+1; %See if we need to decrease r
        Sig = diag(Sig);
        if rnew ~= r
            r = rnew; 
            fprintf('-- Adaptive: Reducing r to %d\n',r);
            S = Sig(1:r,1:r);
            C = C(:,1:r);
            D = D(:,1:r);
        end
    end
    
end
    
end


function Y = firstOrderResidual(C,S,D,Acell,F,dt,X,trans,C1,D1)
%Compute action of first order residual
%R(U) = (I-P_H)U - dt*(L(U)-P_H L(P_H U)) - dt*(I-P_H)F 
NN = size(Acell,1);

S_ = (C1'*C)*S*(D'*D1);
%S_ = (C1'*C)*(0*S)*(D'*D1);


if strcmp(trans,'notransp')
    % Y = U*X;
    Y = C*S*(D'*X);
    % Y = Y - P_HU*X
    Y = Y - C1*S_*(D1'*X);
    % Y = - (LU)*X
    for l=1:NN
        Y = Y - dt*(Acell{l,1}*C)*S*(Acell{l,2}*D)'*X;
    end
    Y = Y - dt*F*X;
    %Projection
    %Y = Y + CC'*(LU)*DD'*X
    for l=1:NN
        Y = Y + dt*C1*(C1'*Acell{l,1}*C1)*S_*(D1'*Acell{l,2}'*D1)*(D1'*X);
    end
    Y = Y + dt*C1*(C1'*F*D1)*(D1'*X);
else
    % Y = U'*X;
    Y = D*S'*(C'*X);
    % Y = Y - (P_HU)'*X
    Y = Y - D1*S_'*(C1'*X);
    % Y = Y - (LU)'*X
    for l=1:NN
        Y = Y - dt*(Acell{l,2}*D)*S'*(Acell{l,1}*C)'*X;
    end
    Y = Y - dt*F'*X;
    %Projection
    %Y = Y + DD'*(LU)*CC'*X
    for l=1:NN
        Y = Y + dt*D1*(D1'*Acell{l,2}*D1)*S_'*(C1'*Acell{l,1}'*C1)*(C1'*X);
    end
    Y = Y + dt*D1*(D1'*F'*C1)*(C1'*X);
end
    
    
end
