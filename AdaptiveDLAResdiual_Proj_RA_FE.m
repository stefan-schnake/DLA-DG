function [C,S,D] = AdaptiveDLAResdiual_Proj_RA_FE(x,v,k,C,S,D,dt,tol,Awave,BC)
%Adaptive algorithm for DLA update

%Here DLA update is U^{n+1} = P_H(U^n + dt*F(U^n))
%where H = X_1ZV_1^T:Z\in\R^{r\times r} 
%and X_1,V_1 are created vs DLA_UC update

%Parameters
p = 3;
pp = p+7;
delta = 0.2;

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C0 = C; S0 = S; D0 = D;
%Update
tic
[C,S,D] = DLA4_HB_PROJ_FE(x,v,k,C,S,D,dt,Awave,BC);
%[C,S,D] = DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC);
%[C,S,D] = DLA5_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%Compute SVD of S
%[S_C,Sig,S_D] = svd(S);
toc
%fprintf(' DLA update time: %f\n',toc);
Rnt = @(x) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,'notransp',C,D);
Rt = @(x) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,'transp',C,D);
tic
Q = rand(size(C,1),pp);
Q = Rnt(Q);
[Q,~] = qr(Q,0);
B = Rt(Q)';
[Utmp,S1,D1] = svd(B,'econ');
C1 = Q*Utmp;
time1 = toc;
normF = norm(diag(S1),2) + sqrt(max([6*r-p,1]))*S1(end,end);
%norm2 = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D),[1,1]*size(C,1),p);
%fprintf('----!!! error is %e \n',norm(diag(S1(1:p,1:p)-norm2(1:p))));
fprintf('-- Adaptive: Residual computation time: %f\n',time1);
fprintf('-- Adaptive: Residual estimate of %e with tol %e\n',normF,tol);
if normF > tol %Add rank
    fprintf('-- Adaptive: Increasing rank\n');
    tolflag = true;
    dS1 = diag(S1);
    for i=2:p %Seeing if the vectors I created are sufficient
        normF = norm(dS1(i:p),2) + sqrt(max([6*r-i,1]))*S1(p,p);
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
    S1 = S1(1:q,1:q);
    C1 = C1(:,1:q);
    D1 = D1(:,1:q);
    while tolflag == true
        if rr+p >= size(C,1)
            tolflag = false;
        end       
        Q = rand(size(C,1),pp);
        Q = Rnt(Q) - LReval(C1,S1,D1,Q,'notransp');
        [Q,~] = qr(Q,0);
        B = Rt(Q)' - LReval(C1,S1,D1,Q,'transp')';
        [Utmp,St,Dt] = svd(B,'econ');
        Ct = Q*Utmp;
        for i=1:p %Seeing if the vectors I created are sufficient
            dSt = diag(St);
            normF = norm(dSt(i:p),2) + sqrt(max([6*r-rr-i,1]))*St(p,p);
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
        S1 = diag([diag(S1);diag(St(1:q,1:q))]);
        C1 = [C1 Ct(:,1:q)];
        D1 = [D1 Dt(:,1:q)];
    end
    fprintf('-- Adaptive: Residual threshold requires %d new basis vectors\n',rr);
    [C,R_C] = qr([C C1],0);
    [D,R_D] = qr([D D1],0);
    S_fill = zeros(r,rr);
    R = R_C*[S S_fill;S_fill' S1]*R_D';
    [S_C,Sig,S_D] = svd(R);
    r = sum(diag(Sig) > 1e-14);
    C = C*S_C(:,1:r);
    S = Sig(1:r,1:r);
    D = D*S_D(:,1:r);
    fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,S(end,end));
elseif normF < delta*tol %Need to decrease rank
    fprintf('-- Adaptive: Decreasing rank\n');
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
    fprintf('-- Adaptive: No change in rank\n');
end
end

function Y = firstOrderResidual(C,S,D,Acell,dt,BC,X,trans,C1,D1)
%Compute action of first order residual R*Q or R^T*Q where
%R(U) = (I-P_H)(U + dt*F(U))
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

function Y = AtranspA(C,S,D,Acell,dt,BC,X,C1,D1)
    Y = firstOrderResidual(C,S,D,Acell,dt,BC,X,'notransp',C1,D1);
    Y = -firstOrderResidual(C,S,D,Acell,dt,BC,Y,  'transp',C1,D1);
end

function Y = LReval(C,S,D,X,trans)

if strcmp(trans,'notransp')
    Y = C*S*(D'*X);
else
    Y = D*S'*(C'*X);
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

function [c,d] = ALTSVD(C0,S0,D0,C,S,D,Acell,dt,BC)

    %Seed c and d
    c = rand(size(C,1),1); c = c/norm(c,2);
    d = rand(size(D,1),1); d = d/norm(d,2);
    
    for i=1:1000
        %Update c
        
        R = (C0*S0)*(D0'*[D d]);
        for l=1:size(Acell,1)
            R = R - dt*(Acell{l,1}*C0)*S0*(D0'*Acell{l,2}'*[D d]);
        end
        R = R - C*S*(D'*[D d]);
        %R = R - C*(C'*R);
        [U,~,~] = svd(R,0);
        c1 = U(:,1);
        %c1 = c1 - C*(C'*c1);
        
        
        %Update d
        R = (D0*S0')*(C0'*[C c]);
        for l=1:size(Acell,1)
            R = R - dt*(Acell{l,2}*D0)*S0'*(C0'*Acell{l,1}'*[C c]);
        end
        R = R - D*S'*(C'*[C c]);
        %R = R - D*(D'*R);
        [U,~,~] = svd(R,0);
        d1 = U(:,1);
        %d1 = d1 - D*(D'*d1);
        
            
        c = c1; c = c/norm(c,2);
        d = d1; d = d/norm(d,2);
    end

end

