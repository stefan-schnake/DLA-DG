function [C,S,D] = AdaptiveDLAResdiual3_FE(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update
tic

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C0 = C; S0 = S; D0 = D;
%Update
[C,S,D] = DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%Compute SVD of S
%[S_C,Sig,S_D] = svd(S);
%fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);
%rflag =  svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D),[1,1]*size(C,1),1);
%fprintf('-- Adaptive: Largest  singular value of residual is %e\n',rflag);
%Compute estimation to residual
PHS0 = (C'*C0)*S0*(D0'*D);

%Add on remaining residual portion
%G = dt*( PHLU(C0,S0,D0,Awave,BC,C,D) - PHLU(C,PHS0,D,Awave,BC,C,D) );
S = PHS0 + dt*PHLU(C0,S0,D0,Awave,BC,C,D);
[S_C,Sig,S_D] = svd(S);
fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);

Afun = @(X) eigAfun(X,C0,S0,D0,Awave,dt,BC,C,D,1);
tic
[CC,sigma] = eigs(Afun,size(C,1),64);
%Reorthogonalize
CC = CC - C*(C'*CC);
sigma = sqrt(diag(sigma));
DD = eigAfun(CC,C0,S0,D0,Awave,dt,BC,C,D,0);
for j=1:size(DD,2)
    DD(:,j) = DD(:,j)/norm(DD(:,j),2);
end
toc


GC = firstOrderResidual(C0,S0,D0,Awave,dt,BC,D,'notransp',C,S,D);
[C2,sig.GC1] = svd([C GC],'econ');%qr([C GC],0);
GD = firstOrderResidual(C0,S0,D0,Awave,dt,BC,C,'transp',C,S,D);
[D2,sig.GC2,~] = svd([D GD],'econ');%qr([D GD],0);
S2 = (C2'*C0)*S0*(D0'*D2) + dt*PHLU(C0,S0,D0,Awave,BC,C2,D2);

[C3,~] = qr([C0 C],0);
[D3,~] = qr([D0 D],0);
S3 = (C3'*C0)*S0*(D0'*D3) + dt*PHLU(C0,S0,D0,Awave,BC,C3,D3);


U_FE = computeFullGrid(C0,S0,D0,Awave,dt,BC);
%PH_U_FE = C*(C'*U_FE*D)*D';
norm_vec.a = norm(C*S*D'-U_FE,'fro');
%nb = norm(C*S1*D'-U_FE,'fro');
norm_vec.b = norm(C2*S2*D2'-U_FE,'fro');
na = norm_vec.a;
nb = norm_vec.b;
norm_vec.c = norm(C3*S3*D3'-U_FE,'fro');
fprintf("|U_FE - U_UC|_F = %e, |U_FE - C2*S2*D2'|_F = %e, %d\n",na,nb,nb <= na);

%S1 = (C'*C0)*S0*(D0'*D) - dt*PHLU(C0,S0,D0,Awave,BC,C,D);


if Sig(end,end) < tol
    rr = r;
    r = sum(diag(Sig) > tol) + 1;
    if rr ~= r
        fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,Sig(end,end));
        C = C*S_C(:,1:r);
        S = Sig(1:r,1:r);
        D = D*S_D(:,1:r);
    end
else
    
    rr = 0;
    tolflag = true;
    tic
    while tolflag == true
        rr = rr + 64;
        if rr >= size(C,1)
            rr = size(C,1);
            tolflag = false;
        end
        [C1,S1,D1] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,S,D),[1,1]*size(C,1),rr);
        if S1(end,end) < tol
            tolflag = false;
        end
        tolflag = false;
    end
    % % RR = computeFullGrid(C0,S0,D0,Awave,dt,BC) - C*S*D';
    % % [C1,S1,D1] = svd(RR);
    % % rr = sum(diag(S1) > tol) + 1; C1 = C1(:,1:rr); S1 = S1(1:rr,1:rr); D1 = D1(:,1:rr);
    %fprintf('-- Adaptive: Residual threshold requires %d new basis vectors\n',rr);
    [C,R_C] = qr([C C1],0);
    [D,R_D] = qr([D D1],0);
    S_fill = zeros(r,rr);
    R = R_C*[S S_fill;S_fill' S1]*R_D';
    toc;
    R1 = (C'*C0)*S0*(D0'*D) + dt*PHLU(C0,S0,D0,Awave,BC,C,D);
    [S_C,Sig,S_D] = svd(R);
    r = sum(diag(Sig) > tol) + 1;
    norm_vec.d = norm(C*R*D'-U_FE,'fro');
    norm_vec.e = norm(C*R1*D'-U_FE,'fro');
    fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,Sig(end,end));
    C = C*S_C(:,1:r);
    S = Sig(1:r,1:r);
    D = D*S_D(:,1:r);
end

toc
end

function Z = PHLU(C,S,D,Acell,BC,C1,D1)
    %Applies P_HL(U) where U=CSD'
    Z = zeros(size(C1,2));
    NN = size(Acell,1);
    
    for l=1:NN
        Z = Z - (C1'*Acell{l,1}*C)*S*(D'*Acell{l,2}'*D1);
    end
    if BC.use
        for l=1:size(BC.cell1)
            Z = Z - (C1'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*D1);
        end
    end
end

function Y = firstOrderResidual(C,S,D,Acell,dt,BC,X,trans,C1,S1,D1)
%Compute action of first order residual
%R(U) = (I-P_H)( U - dt (L(U)-F) ) 
%where F = (F(t_0))
NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

%S_ = (C1'*C)*S*(D'*D1);

if strcmp(trans,'notransp')
    %Y = zeros(size(X));
    
    %Y = U
    Y = C*S*(D'*X);
    %Y = Y - C1*S_*(D1'*X);
    
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
    %Y = Y - P_H(U-dt*(L(U)-F))
    Y = Y - (C1*S1)*(D1'*X);
else
    %Y = zeros(size(X));
    Y = (D*S')*(C'*X);
    %Y = Y - D1*S_'*(C1'*X);
    % Y = Y - (LU)'*X
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
    Y = Y - (D1*S1')*(C1'*X);
    
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

function Y = eigAfun(X,C0,S0,D0,Acell,dt,BC,C,D,swit)
% Computes 
%     (I-CC^T)Z(I-DD^T)Z^TX if swit=1
% or 
%     (I-DD^T)Z^TX if swit=0
% where Z=U_FE if 
    NN = size(Acell,1);
    if BC.use
        MM = size(BC.cell1,1);
    end
    
    %Get %Z^TX = U^TX - dt*L(U)^TX -dt*BC'*X
    Y = (D0*S0')*(C0'*X);
    for l=1:NN
        Y = Y - dt*(Acell{l,2}*D0)*S0'*(C0'*Acell{l,1}'*X);
    end
    if BC.use
        for l=1:MM
            Y = Y - dt*BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*X);
        end
    end
    
    %Now get (I-DD^T)Z^TX = Z^TX - DD^TZ^TX
    Y = Y - D*(D'*Y);
    if swit == 1
        YY = Y;

        %Now apply Z to it
        Y = (C0*S0)*(D0'*YY);
        for l=1:NN
            Y = Y - dt*(Acell{l,1}*C0)*S0*(D0'*Acell{l,2}'*YY);
        end
        if BC.use
            for l=1:MM
                Y = Y - dt*BC.cell1{l,1}*BC.cell1{l,2}'*(BC.cell1{l,3}'*YY);
            end
        end

        %Now apply (I-CC^T)
        Y = Y - C*(C'*Y);
    end
end

