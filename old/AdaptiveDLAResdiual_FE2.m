function [C,S,D] = AdaptiveDLAResdiual_FE2(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update


r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
%[Cit,Sit,Dit] = DLA3_HB_FE_Iteration(x,v,k,C_old,S_old,D_old,dt,Awave,BC,FMWT);
[C,S,D] = DLA4_HB_FE(x,v,k,C_old,S_old,D_old,dt,Awave,BC);
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
elseif r < size(C,1) || 1 %Add most relavant vectors from residual
    %%%First order residual
    fprintf('-- Adaptive: Adding first order residual basis\n');
    Chalf = C;
    Dhalf = D;
    C0 = C_old; D0 = D_old; S0 = S_old;
    %[Chalf,~] = sRRQR([C0 Chalf],1.01,'tol',1e-12);
    %[Dhalf,~] = sRRQR([D0 Dhalf],1.01,'tol',1e-12);
    %%%Compute first order residual basis vectors
    if BC.use
        %FF = FMWT*convertVectoMat(x,v,k,BC(:,1))*FMWT';
        FF = BC.cell1{1,1}*BC.cell1{1,2}*BC.cell1{1,3}';
    else
        FF = zeros(size(C,1));
    end
    
    [C1,~] = qr([Chalf C0],0);
    [D1,~] = qr([Dhalf D0],0);
    
    R0.mat = firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',Chalf,Dhalf,0);
    R0 = compSVD(R0);
    R1.mat = firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',Chalf,Dhalf,1);
    R1 = compSVD(R1);
    R2.mat = firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',R0.C(:,1:r),R0.D(:,1:r),1);
    R2 = compSVD(R2);
    %R3.mat = firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',Cit,Dit,1);
    %R3 = compSVD(R3);
    
    %rr = 5;
    
    %S1 = DLA_HB_FE_S_Only(x,v,k,Chalf      ,Chalf'*C0*S0*D0'*Dhalf            ,Dhalf      ,dt,Awave,BC,FMWT);
    %S2 = DLA_HB_FE_S_Only(x,v,k,R0.C(:,1:rr),R0.C(:,1:rr)'*C0*S0*D0'*R0.D(:,1:rr),R0.D(:,1:rr),dt,Awave,BC,FMWT);
    
    %[C2,~] = qr([Chalf R1.C(:,1:rr)],0);
    %[D2,~] = qr([Dhalf R1.D(:,1:rr)],0);
    %S3 = DLA_HB_FE_S_Only(x,v,k,C2,C2'*C0*S0*D0'*D2,D2,dt,Awave,BC,FMWT);
    
    %UU.mat = C*S*D' + R1.C(:,1:rr)*R1.S(1:rr,1:rr)*R1.D(:,1:rr)';
    %UU = compSVD(UU);
    %UU2.mat = C2*S3*D2';
    %U2 = compSVD(UU2);
    
    
    U_FE.mat = C0*S0*D0' - dt*FF;
    for l=1:size(Awave,1)
        U_FE.mat = U_FE.mat - dt*(Awave{l,1}*C0)*S0*(Awave{l,2}*D0)';
    end
    U_FE = compSVD(U_FE);
    
    fprob.objective = @(x) myobj(x,Chalf,Dhalf,C0*S0*D0',dt,Awave,FF);
    rand1 = rand(size(Chalf,1),1);
    Crand = rand1-Chalf*Chalf'*rand1; Drand = rand1-Dhalf*Dhalf'*rand1;
    fprob.x0 = reshape([Crand/norm(Crand) Drand/norm(Drand)],[],1);
    fprob.Aineq = [];
    fprob.bineq = [];
    fprob.Aeq = [Chalf' Dhalf'];
    fprob.beq = zeros(size(Chalf,2),1);
    fprob.lb = [];
    fprob.ub = [];
    fprob.nonlcon = @(x) nonlcon(x,C0*S0*D0');
    fprob.solver = 'fmincon';
    fprob.options = optimoptions(@fmincon,'display','iter-detailed');
    fprob.options = optimoptions(fprob.options,'MaxFunctionEvaluations',1e5);
    fprob.options = optimoptions(fprob.options,'OptimalityTolerance',1e-10);
    fprob.options = optimoptions(fprob.options,'ConstraintTolerance',1e-10);
    fprob.options = optimoptions(fprob.options,'FiniteDifferenceType','central');
    
    xx = fmincon(fprob);
    R3 = extractVec(xx,size(FF,1));
    [R3.C,~] = qr([Chalf R3.C],0);
    [R3.D,~] = qr([Dhalf R3.D],0);
    R4.mat = firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',R3.C,R3.D,1);
    R4 = compSVD(R4);
%     S3 = DLA_HB_FE_S_Only(x,v,k,R3.C,R3.C'*C0*S0*D0'*R3.D,R3.D,dt,Awave,BC,FMWT);
%     
%     fprob.x0 = reshape([C0(:,2) D0(:,2)],[],1);
%     fprob.Aeq = [R3.C' R3.D'];
%     fprob.beq = 0;
%     xx = fmincon(fprob);
%     R4 = extractVec(xx,size(FF,1));
%     [R4.C,~] = qr([R3.C R4.C],0);
%     [R4.D,~] = qr([R3.D R4.D],0);
%     S4 = DLA_HB_FE_S_Only(x,v,k,R4.C,R4.C'*C0*S0*D0'*R4.D,R4.D,dt,Awave,BC,FMWT);
    
    
    %[C1,S1,D1] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,FF,x,trans,Chalf,Dhalf),[1,1]*size(C,1),15);
    %[C1,S1,D1] = svd(firstOrderResidual(C0,S0,D0,Awave,FF,dt,eye(size(Chalf,1)),'notransp',Chalf,Dhalf));
    %Prune irrelavent vectors
    %rr = sum(diag(S1) > 1e-16);
    %fprintf('-- Adaptive: Found %d non-trivial vectors to add\n',rr);
    %C1 = C1(:,1:rr); S1 = S1(1:rr,1:rr); D1 = D1(:,1:rr);
    
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


function Y = firstOrderResidual(C,S,D,Acell,F,dt,X,trans,C1,D1,projflag)
%Compute action of first order residual
%R(U) = (I-P_H)U - dt*(L(U)-P_H L(P_H U)) - dt*(I-P_H)F 
NN = size(Acell,1);

if projflag
S_ = (C1'*C)*S*(D'*D1);
end

if strcmp(trans,'notransp')
    % Y = U*X;
    Y = C*S*(D'*X);
    if projflag
        % Y = Y - P_HU*X
        Y = Y - C1*S_*(D1'*X);
    end
    % Y = - (LU)*X
    for l=1:NN
        Y = Y - dt*(Acell{l,1}*C)*S*(Acell{l,2}*D)'*X;
    end
    Y = Y - dt*F*X;
    if projflag
        %Projection
        %Y = Y + CC'*(LU)*DD'*X
        for l=1:NN
            Y = Y + dt*C1*(C1'*Acell{l,1}*C1)*S_*(D1'*Acell{l,2}'*D1)*(D1'*X);
            %Y = Y + dt*C1*(C1'*Acell{l,1}*C)*S*(D'*Acell{l,2}'*D1)*(D1'*X);
        end
        Y = Y + dt*C1*(C1'*F*D1)*(D1'*X);
    end
else
    % Y = U'*X;
    Y = D*S'*(C'*X);
    if projflag
        % Y = Y - (P_HU)'*X
        Y = Y - D1*S_'*(C1'*X);
    end
    % Y = Y - (LU)'*X
    for l=1:NN
        Y = Y - dt*(Acell{l,2}*D)*S'*(Acell{l,1}*C)'*X;
    end
    Y = Y - dt*F'*X;
    if projflag
        %Projection
        %Y = Y + DD'*(LU)*CC'*X
        for l=1:NN
            Y = Y + dt*D1*(D1'*Acell{l,2}*D1)*S_'*(C1'*Acell{l,1}'*C1)*(C1'*X);
        end
        Y = Y + dt*D1*(D1'*F'*C1)*(C1'*X);
    end
end
    
    
end

function U = compSVD(U)
    [U.C,U.S,U.D] = svd(U.mat);
end

function LU = applymL(U,Acell)
    LU = 0*U;
    for i=1:size(Acell,1)
        LU = LU - Acell{i,1}*U*Acell{i,2}';
    end
end

function z = myobj(x,Chalf,Dhalf,U,dt,Acell,FF)
    n = size(U,1);
    X = reshape(x,n,[]);
    r = size(X,2)/2;
    C = [Chalf X(:,1:r)];
    D = [Dhalf X(:,r+1:2*r)];
    PhU = C*(C'*U*D)*D';
    Z = (U-PhU) - dt*(applymL(U,Acell) - C*(C'*applymL(PhU,Acell)*D)*D') - dt*(FF-C*(C'*FF*D)*D');
    z = -norm(Z,'fro')^2;
end

function [c,ceq] = nonlcon(x,U)
    n = size(U,1);
    X = reshape(x,n,[]);
    r = size(X,2)/2;
    C = X(:,1:r);
    D = X(:,r+1:2*r);
    I_C = C'*C;
    I_D = D'*D;
    c = [];
    ceq = norm(I_C-eye(r),'fro')^2+norm(I_D-eye(r),'fro')^2;
end

function R = extractVec(x,n,R)
    if nargin < 3
        R.C = [];
        R.D = [];
    end
    X = reshape(x,n,[]);
    r = size(X,2)/2;
    R.C = X(:,1:r);
    R.D = X(:,r+1:2*r);
end

function Z = tempRE(U,C,D,Acell)
    %Computes P_H L(U-P_H U)
    PhU = C*(C'*U*D)*D';
    UmPhU = U-PhU;
    Z = zeros(size(PhU));
    for l=1:size(Acell,1)
        %Z = Z + C*(C'*(Acell{l,1}*UmPhU*Acell{l,2}')*D)*D';
        Z = Z + Acell{l,1}*UmPhU*Acell{l,2}';
    end
    
end
