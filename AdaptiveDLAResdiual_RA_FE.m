function [U] = AdaptiveDLAResdiual_RA_FE(U,dt,tol,Awave,BC)
%Adaptive algorithm for DLA update

%Parameters
p = 3;
pp = p+7;
delta = 0.7;

quiet = false;

r = size(S,1);
if quiet; fprintf('-- Adaptive: r = %d\n',r); end
U0 = SOLN(U.x,U.v,U.k);
R  = SOLN(U.x,U.v,U.k);
R_tmp = SOLN(U.x,U.v,U.k);
U0.C = U.C; U0.S = U.S; U0.D = U.D;
%Update
tic
U = DLA4_HB_FE(U,dt,Awave,BC);
%[C,S,D] = DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC);
%[C,S,D] = DLA5_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%Compute SVD of S
%[S_C,Sig,S_D] = svd(S);
time1 = toc;
%fprintf(' DLA update time: %f\n',toc);
Rnt = @(x) firstOrderResidual(U0,Awave,dt,BC,x,'notransp',U);
Rt = @(x) firstOrderResidual(U0,Awave,dt,BC,x,'transp',U);
tic
Q = randn(size(U.C,1),pp);
Q = Rnt(Q);
[Q,~] = qr(Q,0);
B = Rt(Q)';
[Utmp,R.S,R.D] = svd(B,'econ');
R.C = Q*Utmp;
time1 = toc;
normF = norm(diag(R.S),2) + sqrt(max([6*r-p,1]))*R.S(end,end);
%norm2 = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D),[1,1]*size(C,1),p);
%fprintf('----!!! error is %e \n',norm(diag(S1(1:p,1:p)-norm2(1:p))));
if quiet; fprintf('-- Adaptive: Residual computation time: %f\n',time1); end
if quiet; fprintf('-- Adaptive: Residual estimate of %e with tol %e\n',normF,tol); end
if normF > tol %Add rank
    if quiet; fprintf('-- Adaptive: Increasing rank\n'); end
    tolflag = true;
    dS1 = diag(R.S);
    for i=2:p %Seeing if the vectors I created are sufficient
        normF = norm(dS1(i:pp),2) + sqrt(max([6*r-i,1]))*R.S(pp,pp);
        if normF < tol
            tolflag = false;
            break
        end
    end
    if tolflag
        q = i;
    else
        q = i-1;
    end
    rr = q;
    R.S = R.S(1:q,1:q);
    R.C = R.C(:,1:q);
    R.D = R.D(:,1:q);
    while tolflag == true
        if rr+p >= size(R.C,1)
            tolflag = false;
        end       
        Q = randn(size(R.C,1),pp);
        Q = Rnt(Q) - LReval(R,Q,'notransp');
        [Q,~] = qr(Q,0);
        B = Rt(Q)' - LReval(R,Q,'transp')';
        [Utmp,R_tmp.S,R_tmp.D] = svd(B,'econ');
        R_tmp.C = Q*Utmp;
        for i=1:p %Seeing if the vectors I created are sufficient
            dSt = diag(R_tmp.S);
            normF = norm(dSt(i:pp),2) + sqrt(max([6*r-rr-i,1]))*R_tmp.S(pp,pp);
            if normF < tol
                tolflag = false;
                break
            end
        end
        if tolflag
            q = i;
        else
            q = i-1;
        end
        rr = rr + q;
        R.S = diag([diag(R.S);diag(R_tmp.S(1:q,1:q))]);
        R.C = [R.C R_tmp.C(:,1:q)];
        R.D = [R.D R_tmp.D(:,1:q)];
    end
    if quiet; fprintf('-- Adaptive: Residual threshold requires %d new basis vectors\n',rr); end
    [U.C,R_C] = qr([U.C R.C],0);
    [U.D,R_D] = qr([U.D R.D],0);
    S_fill = zeros(r,rr);
    U.S = R_C*[S S_fill;S_fill' R.S]*R_D';
    [S_C,Sig,S_D] = svd(U.S);
    r = sum(diag(Sig) > 1e-14);
    U.C = U.C*S_C(:,1:r);
    U.S = Sig(1:r,1:r);
    U.D = U.D*S_D(:,1:r);
    if quiet; fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,S(end,end)); end
elseif normF < delta*tol %Need to decrease rank
    if quiet; fprintf('-- Adaptive: Decreasing rank\n'); end
    [S_C,Sig,S_D] = svd(S);
    Sig = diag(Sig);
    rr = r;
    for i=r:-1:1
        if sum(Sig(i:r)) + normF >= delta*tol
            rr = i;
            break
        end
    end
    if rr ~= r
        r = rr;
        %fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,Sig(end,end));
        C = C*S_C(:,1:r);
        Sig = diag(Sig);
        S = Sig(1:r,1:r);
        D = D*S_D(:,1:r);
    end
else
    if quiet; fprintf('-- Adaptive: No change in rank\n'); end
end
end

function Y = firstOrderResidual(U0,Acell,dt,BC,X,trans,U)
%Compute action of first order residual
%R(U) = (I-P_H)U - dt*(L(U)-P_H L(P_H U) + (I-P_H)F) 
%where F = 0.5*(F(t_0)+F(t_1))
NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

S_ = (U.C'*U0.C)*U0.S*(U0.D'*U.D);

if strcmp(trans,'notransp')
    %Y = zeros(size(X));
    
    Y = U0.C*U0.S*(U0.D'*X);
    Y = Y - U.C*S_*(U.D'*X);
    
    % Y = Y - (LU)*X
    for l=1:NN
        Y = Y - dt*(Acell{l,1}*U0.C)*U0.S*(U0.D'*Acell{l,2}'*X);
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
        Y = Y + dt*U.C*(U.C'*Acell{l,1}*U.C)*S_*(U.D'*Acell{l,2}'*U.D)*(U.D'*X);
    end
    if BC.use
        %Y = Y + U.C*(U.C'*F*U.D)*(U.D'*X);
        for l=1:MM
            Y = Y + dt*U.C*(U.C'*BC.cell1{l,1})*BC.cell1{l,2}*(BC.cell1{l,3}'*U.D)*(U.D'*X);
        end
    end
else
    %Y = zeros(size(X));
    Y = U0.D*U0.S'*(U0.C'*X);
    Y = Y - U.D*S_'*(U.C'*X);
    % Y = -(LU)'*X
    for l=1:NN
        Y = Y - dt*(Acell{l,2}*U0.D)*U0.S'*(U0.C'*Acell{l,1}'*X);
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
        Y = Y + dt*U.D*(U.D'*Acell{l,2}*U.D)*S_'*(U.C'*Acell{l,1}'*U.C)*(U.C'*X);
    end
    if BC.use
        %Y = Y + U.D*(U.D'*F'*U.C)*(U.C'*X);
        for l=1:MM
            Y = Y + dt*U.D*(U.D'*BC.cell1{l,3})*BC.cell1{l,2}'*(BC.cell1{l,1}'*U.C)*(U.C'*X);
        end
    end
    
end
    

end

function Y = LReval(U,X,trans)

if strcmp(trans,'notransp')
    Y = U.C*U.S*(U.D'*X);
else
    Y = U.D*U.S'*(U.C'*X);
end

end

