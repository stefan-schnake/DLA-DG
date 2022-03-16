function [C,S,D] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C,S,D,dt,tol,Awave,BC)
%Adaptive algorithm for DLA update

%Parameters
p = 3;
pp = p+7;
delta = 0.5;

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C0 = C; S0 = S; D0 = D;
%Update
[C,S,D] = DLA4_TAN_FE(x,v,k,C,S,D,dt,Awave,BC);
%fprintf(' DLA update time: %f\n',toc);
Rnt = @(x) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,'notransp');
 Rt = @(x) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,'transp'  );
%% Create randomized SVD
Q = rand(size(C,1),pp);
Q = Rnt(Q);
[Q,~] = qr(Q,0);
B = Rt(Q)';
[Utmp,S1,D1] = svd(B,'econ');
C1 = Q*Utmp;
time1 = toc;

%% Create norm estimator
normF = norm(diag(S1),2) + sqrt(max([6*r-p,1]))*S1(end,end);

fprintf('-- Adaptive: Residual computation time: %f\n',time1);
fprintf('-- Adaptive: Residual estimate of %e with tol %e\n',normF,tol);
if normF > tol %Add rank
    fprintf('-- Adaptive: Increasing rank\n');
    tolflag = true;
    dS1 = diag(S1);
    for i=2:p %Seeing if the vectors I created are sufficient
        normF = norm(dS1(i:pp),2) + sqrt(max([6*r-i,1]))*S1(pp,pp);
        if normF < 0.9*tol
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
            normF = norm(dSt(i:pp),2) + sqrt(max([6*r-rr-i,1]))*St(pp,pp);
            if normF < 0.9*tol
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
    fprintf('-- Adaptive:  for a total of %d basis vectors\n',2*r+rr);
    [C,R_C] = qr([C C1],0);
    [D,R_D] = qr([D D1],0);
    S_fill = zeros(2*r,rr);
    R = R_C*[S S_fill;S_fill' S1]*R_D';
    [S_C,Sig,S_D] = svd(R);
    r = size(Sig,1);
    dSig = diag(Sig);
    for i=r:-1:1
        if norm(dSig(i:r)) > 0.1*tol
            rr = i;
            break
        end
    end
    %r = rr;
    C = C*S_C(:,1:r);
    S = Sig(1:r,1:r);
    D = D*S_D(:,1:r);
    fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,tol);
elseif normF < delta*tol %Need to decrease rank
    fprintf('-- Adaptive: Decreasing rank\n');
    [S_C,Sig,S_D] = svd(S);
    dSig = diag(Sig);
    rr = 2*r;
    for i=2*r:-1:1
        if sqrt(norm(dSig(i:min([2*r,size(C,1)])),2)^2 + normF^2) >= delta*tol
            rr = i+1;
            break
        end
    end
    rr = min([rr,2*r]);
    r = rr;
    %fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,Sig(end,end));
    C = C*S_C(:,1:r);
    S = Sig(1:r,1:r);
    D = D*S_D(:,1:r);
else
    fprintf('-- Adaptive: No change in rank\n');
end
end

function Y = firstOrderResidual(C,S,D,Acell,dt,BC,X,trans)
%Compute action of first order residual
% R(U) = dt*(I-CC^T)*F(U)*(I-DD^T)*X if trans='notransp'
% or 
% R^T(U) = dt*(I-DD^T)*F(U)*(I-CC^T)*X if trans='transp'

NN = size(Acell,1);
if BC.use
    MM = size(BC.cell1);
end

if strcmp(trans,'notransp')
    
    %Project out of span(D);
    X = dt*(X - D*(D'*X));
    
    % Y = F(U)*X
    Y = zeros(size(X));
    for l=1:NN
        Y = Y - (Acell{l,1}*C)*S*(D'*Acell{l,2}'*X);
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - BC.cell1{l,1}*BC.cell1{l,2}*(BC.cell1{l,3}'*X);
        end
    end
    
    %Project out of span(C)
    Y = Y - C*(C'*Y);
else
    
    %Project out of span(C)
    X = dt*(X - C*(C'*X));
    
    % Y = F(U)'*X
    Y = zeros(size(X));
    for l=1:NN
        Y = Y - (Acell{l,2}*D)*S'*(C'*Acell{l,1}'*X);
    end
    if BC.use
        %Y = Y - F*X;
        for l=1:MM
            Y = Y - BC.cell1{l,3}*BC.cell1{l,2}'*(BC.cell1{l,1}'*X);
        end
    end
    
    %Project out span(D)
    Y = Y - D*(D'*Y);
    
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


