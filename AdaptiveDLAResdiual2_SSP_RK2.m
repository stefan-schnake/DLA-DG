function [C,S,D] = AdaptiveDLAResdiual2_SSP_RK2(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update
tic

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
[C,S,D] = DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC);
%Compute SVD of S
[S_x,Sig,S_v] = svd(S);
Sig = diag(Sig);
fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);
if Sig(end) < tol %We don't need to bump up r
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
    %[C1,S1,D1] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,BC,x,trans,Chalf,Dhalf),[1,1]*size(C,1),15);
    [C1,S1,D1] = svds(@(x,trans) firstOrderResidual2(C0,S0,D0,Awave,dt,BC,x,trans,Chalf,Dhalf),[1,1]*size(C,1),15);
    %[C1,S1,D1] = svd(firstOrderResidual(C0,S0,D0,Awave,BC,eye(size(Chalf,1)),'notransp',Chalf,Dhalf));
    rr = sum(diag(S1) > 1e-12);
    fprintf('-- Adaptive: Found %d non-trivial vectors to add\n',rr);
    C1 = C1(:,1:rr); S1 = S1(1:rr,1:rr); D1 = D1(:,1:rr);
    [C2,~,~] = sRRQR([Chalf C1],1.01,'tol',1e-14);
    [D2,~,~] = sRRQR([Dhalf D1],1.01,'tol',1e-14);
    %%% LUBICH ALGORITHM
    %[C2,~,~] = sRRQR([Chalf C0],1.01,'tol',1e-14);
    %[D2,~,~] = sRRQR([Dhalf D0],1.01,'tol',1e-14);
    %%%
    %[C2,~] = qr([Chalf C1],0);
    %[D2,~] = qr([Dhalf D1],0);
    %C2 = [Chalf,C1];
    %D2 = [Dhalf,D1];
    S2 = (C2'*C0)*S0*(D0'*D2);
    fprintf('-- Adaptive: New total of %d basis vectors\n',size(C2,2));
    %[C,S,D] = AdaptiveDLAWithInitv8(x,v,k,C,S,D,t,updateDLA,tol,Awave,BC,FMWT);
    %[C,S,D] = updateDLA(C,S,D,t);
    fprintf('-- Adaptive: Updating coefficients with new basis\n');
    S2 = DLA4_HB_SSP_RK2_S_Only(x,v,k,C2,S2,D2,dt,Awave,BC);
    C = C2; D = D2; S = S2;
    [S_x,Sig,S_v] = svd(S2);
    Sig = diag(Sig);
    r = size(S,1);
    fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);
    %%%Second order residual
    if Sig(end) > tol && 0
        if BC.use
            FF = FMWT*convertVectoMat(x,v,k,BC.vec(:,1))*FMWT';
        else
            FF = zeros(size(C,1));
        end
        fprintf('-- Adaptive: Adding second order residual basis\n');
        RR2 = secondOrderResidual(C0,S0,D0,Awave,FF,eye(size(C,1)),'notransp',Chalf,Dhalf,C1,S1,D1);
        [C3,S3,D3] = svd(RR2);
        %[C3,S3,D3] = svds(@(x,trans) secondOrderResidual(C0,S0,D0,Awave,FF,x,trans,Chalf,Dhalf,C1,S1,D1),[1,1]*size(C,1),50);
        %[C3,S3,D3] = svds(@(x,trans) secondOrderResidual(C0,S0,D0,Awave,FF,x,trans,C2,D2,C1,0*S1,D1),[1,1]*size(C,1),30);
        rr = sum(diag(S3) > 1e-12);
        fprintf('-- Adaptive: Found %d non-trivial vectors to add\n',rr);
        C3 = C3(:,1:rr); S3 = S3(1:rr,1:rr); D3 = D3(:,1:rr);
        [C4,~,~] = sRRQR([C2 C3],1.01,'tol',1e-14);
        [D4,~,~] = sRRQR([D2 D3],1.01,'tol',1e-14);
%         [C4,R1.mat] = qr([C2 C3],0);
%         [D4,R2.mat] = qr([D2 D3],0);
%         %Remove duplicates
%         [R1.U,R1.S,R1.V] = svd(R1.mat);
%         R1.rank = sum(diag(R1.S) > 1e-12);
%         [R2.U,R2.S,R2.V] = svd(R2.mat);
%         R2.rank = sum(diag(R2.S) > 1e-12);
%         rr = max([R1.rank R2.rank]);
%         C4 = C4*R1.U(:,1:rr);
%         D4 = D4*R2.U(:,1:rr);
        fprintf('-- Adaptive: New total of %d basis vectors\n',size(C4,2));
        S4 = (C4'*C0)*S0*(D0'*D4);
        S4 = DLA3_HB_SSP_RK2_S_Only(x,v,k,C4,S4,D4,dt,Awave,BC,FMWT);
        [S_x,Sig,S_v] = svd(S4);
        Sig = diag(Sig);
        fprintf('-- Adaptive: Smallest singular value is %e\n',Sig(end));
        C = C4; D = D4; S = S4;
    end
    
    if (Sig(end) < tol) %We don't need to bump up r
        rnew = sum(Sig > tol)+1; %See if we need to decrease r
        Sig = diag(Sig);
        if rnew ~= r
            r = rnew; 
            fprintf('-- Adaptive: Reducing r to %d\n',r);
            S = Sig(1:r,1:r);
            C = C*S_x(:,1:r);
            D = D*S_v(:,1:r);
        end
    end
    
end
toc
end

function Y = firstOrderResidual(C,S,D,Acell,BC,X,trans,C1,D1)
%Compute action of first order residual
%R(U) = L(U)-P_H L(P_H U) - (I-P_H)F 
%where F = 0.5*(F(t_0)+F(t_1))
NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

S_ = (C1'*C)*S*(D'*D1);
%S_ = (C1'*C)*(0*S)*(D'*D1);


if strcmp(trans,'notransp')
    Y = zeros(size(X));    
    % Y = (LU)*X
    for l=1:NN
        Y = Y - (Acell{l,1}*C)*S*(Acell{l,2}*D)'*X;
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - 0.5*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*X);
            Y = Y - 0.5*BC.cell3{l,1}*BC.cell3{l,2}*(BC.cell3{l,3}'*X);
        end
    end
    %Projection
    %Y = Y - CC'*(LU)*DD'*X
    for l=1:NN
        Y = Y + C1*(C1'*Acell{l,1}*C1)*S_*(D1'*Acell{l,2}'*D1)*(D1'*X);
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y + 0.5*C1*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1)*(D1'*X);
            Y = Y + 0.5*C1*(C1'*BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D1)*(D1'*X);
        end
    end
    %Y = Y + C1*(C1'*F*D1)*(D1'*X);
else
    Y = zeros(size(X));
    % Y = (LU)'*X
    for l=1:NN
        Y = Y - (Acell{l,2}*D)*S'*(Acell{l,1}*C)'*X;
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - 0.5*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*X);
            Y = Y - 0.5*BC.cell3{l,3}*BC.cell3{l,2}'*(BC.cell3{l,1}'*X);
        end
    end
    %Projection
    %Y = Y - DD'*(LU)*CC'*X
    for l=1:NN
        Y = Y + D1*(D1'*Acell{l,2}*D1)*S_'*(C1'*Acell{l,1}'*C1)*(C1'*X);
    end
    if BC.use
        %Y = Y + D1*(D1'*F'*C1)*(C1'*X);
        for l=1:MM
            Y = Y + 0.5*D1*(D1'*BC.cell1{l,3})*BC.cell1{l,2}'*(BC.cell1{l,1}'*C1)*(C1'*X);
            Y = Y + 0.5*D1*(D1'*BC.cell3{l,3})*BC.cell3{l,2}'*(BC.cell3{l,1}'*C1)*(C1'*X);
        end
    end
    
end
    

end

function Y = firstOrderResidual2(C,S,D,Acell,dt,BC,X,trans,C1,D1)
%Compute action of first order residual
%R(U) = L(U)-P_H L(P_H U) - (I-P_H)F 
%where F = 0.5*(F(t_0)+F(t_1))
NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

S_ = (C1'*C)*S*(D'*D1);
%S_ = (C1'*C)*(0*S)*(D'*D1);


if strcmp(trans,'notransp')
    %Y = zeros(size(X));
    
    Y = C*S*(D'*X);
    Y = Y - C1*S_*(D1'*X);
    
    % Y = (LU)*X
    for l=1:NN
        Y = Y - dt*(Acell{l,1}*C)*S*(Acell{l,2}*D)'*X;
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - dt*0.5*BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*X);
            Y = Y - dt*0.5*BC.cell3{l,1}*BC.cell3{l,2}*(BC.cell3{l,3}'*X);
        end
    end
    %Projection
    %Y = Y - CC'*(LU)*DD'*X
    for l=1:NN
        Y = Y + C1*(C1'*Acell{l,1}*C1)*S_*(D1'*Acell{l,2}'*D1)*(D1'*X);
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y + dt*0.5*C1*(C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1)*(D1'*X);
            Y = Y + dt*0.5*C1*(C1'*BC.cell3{l,1})*BC.cell3{l,2}*(BC.cell3{l,3}'*D1)*(D1'*X);
        end
    end
    %Y = Y + C1*(C1'*F*D1)*(D1'*X);
else
    %Y = zeros(size(X));
    Y = D*S'*(C'*X);
    Y = Y - D1*S_'*(C1'*X);
    % Y = (LU)'*X
    for l=1:NN
        Y = Y - dt*(Acell{l,2}*D)*S'*(Acell{l,1}*C)'*X;
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - dt*0.5*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*X);
            Y = Y - dt*0.5*BC.cell3{l,3}*BC.cell3{l,2}'*(BC.cell3{l,1}'*X);
        end
    end
    %Projection
    %Y = Y - DD'*(LU)*CC'*X
    for l=1:NN
        Y = Y + D1*(D1'*Acell{l,2}*D1)*S_'*(C1'*Acell{l,1}'*C1)*(C1'*X);
    end
    if BC.use
        %Y = Y + D1*(D1'*F'*C1)*(C1'*X);
        for l=1:MM
            Y = Y + dt*0.5*D1*(D1'*BC.cell1{l,3})*BC.cell1{l,2}'*(BC.cell1{l,1}'*C1)*(C1'*X);
            Y = Y + dt*0.5*D1*(D1'*BC.cell3{l,3})*BC.cell3{l,2}'*(BC.cell3{l,1}'*C1)*(C1'*X);
        end
    end
    
end
    

end


function Y = secondOrderResidual(C,S,D,Acell,F,X,trans,Chalf,Dhalf,C1,S1,D1)
%Compute action of first order residual
%R(U) = L( L(U)-P_H L(P_H U) ) + (I-P_H)L( P_H L(P_H U) ) - ( L(F) - P_H L(P_H F_0) )
%where F = F(t_0)

%Here L(U)-P_H L(P_H U)=R_0 has already been computed in low-rank fashion by
%firstOrderResidual as R_0 = C1*S*D1'

NN = size(Acell,1);

%Project U onto H
S_ = (Chalf'*C)*S*(D'*Dhalf);

%Project L(P_H U) onto H: Chalf*Z1*Dhalf'
Z_1 = zeros(size(S_));
for l=1:NN
    Z_1 = Z_1 + (Chalf'*Acell{l,1}*Chalf)*S_*(Dhalf'*Acell{l,2}*Dhalf);
end
    
%Project F onto H, so P_H F = Chalf*S_F*Dhalf'
S_F = Chalf'*F*Dhalf;

%Actual evaluation
if strcmp(trans,'notransp')
    Y = zeros(size(X));
    % Y = (LR_0)*X
    for l=1:NN
        Y = Y + (Acell{l,1}*C1)*S1*(Acell{l,2}*D1)'*X;
    end
    
    %Y = Y + L(Chalf*Z1*Dhalf')*X;
    for l=1:NN
        Y = Y + (Acell{l,1}*Chalf)*Z_1*(Acell{l,2}*Dhalf)'*X;
    end
    
    %Y = Y - (P_H L(Chalf*Z1*Dhalf'))*X;
    for l=1:NN
        Y = Y - Chalf*(Chalf'*Acell{l,1}*Chalf)*Z_1*(Dhalf'*Acell{l,2}*Dhalf)'*(Dhalf'*X);
    end
    
    %Y = Y - L(F)*X
    for l=1:NN
        Y = Y - Acell{l,1}*F*Acell{l,2}'*X;
    end
    
    %Y = Y + P_H L(P_H F)*X
    for l=1:NN
        Y = Y + Chalf*(Chalf'*Acell{l,1}*Chalf)*S_F*(Dhalf'*Acell{l,2}*Dhalf)'*(Dhalf'*X);
    end
else
    Y = zeros(size(X));
    % Y = (LR_0)'*X
    for l=1:NN
        Y = Y + (Acell{l,2}*D1)*S1'*(Acell{l,1}*C1)'*X;
    end
    
    %Y = Y + L(Chalf*Z1*Dhalf')'*X;
    for l=1:NN
        Y = Y + (Acell{l,2}*Dhalf)*Z_1'*(Acell{l,1}*Chalf)'*X;
    end
    
    %Y = Y - (P_H L(Chalf*Z1*Dhalf'))'*X;
    for l=1:NN
        Y = Y - Dhalf*(Dhalf'*Acell{l,2}*Dhalf)*Z_1'*(Chalf'*Acell{l,1}*Chalf)'*(Chalf'*X);
    end
    
    %Y = Y - L(F)'*X
    for l=1:NN
        Y = Y - Acell{l,2}*F'*Acell{l,1}'*X;
    end
    
    %Y = Y + P_H L(P_H F)'*X
    for l=1:NN
        Y = Y + Dhalf*(Dhalf'*Acell{l,2}*Dhalf)*S_F'*(Chalf'*Acell{l,1}*Chalf)'*(Chalf'*X);
    end
end



end

