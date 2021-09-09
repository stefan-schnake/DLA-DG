function [C,S,D] = AdaptiveDLAResdiual2_FE(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update
tic

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C0 = C; S0 = S; D0 = D;
%Update
tic
[C,S,D] = DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%Compute SVD of S
%[S_C,Sig,S_D] = svd(S);
toc
fprintf(' DLA update time: %f\n',toc);
tic
[C1,norm2,D1] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D),[1,1]*size(C,1),1);
toc
normF = norm2*sqrt(6*r)/2;
% Afun = @(X) AtranspA(C0,S0,D0,Awave,dt,BC,X,C,D);
% X0 = rand(size(C,1),1); X0 = X0/norm(X0);
% tic
% [V2,lm2,failureFlag] = blopex_matlab(X0,Afun,[],[],[],1e-8,20,0);
% toc
% Afun = @(X) -AtranspA(C0,S0,D0,Awave,dt,BC,X,C,D);
% tic
% [V3,lm3] = eigs(Afun,size(C,1),1);
% toc
fprintf('-- Adaptive: Residual computation time: %f\n',toc);
fprintf('-- Adaptive: Residual estimate of %e with tol %e\n',sqrt(6*r)*norm2,tol);
if normF > tol %Add rank
    fprintf('-- Adaptive: Increasing rank\n');
    S1 = norm2;
    rr = 1;
    tolflag = true;
    while tolflag == true
        if rr+5 >= size(C,1)
            tolflag = false;
        end       
        [Ct,St,Dt] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D) ...
            - LReval(C1,S1,D1,x,trans),[1,1]*size(C,1),5,'largest','FailureTreatment','drop');
        rr = rr + size(St,1);
        norm2 = St(end:end);
        normF =  norm2*sqrt(6*r)/2;
        if normF < tol
            tolflag = false;
        end
        S1 = diag([diag(S1);diag(St)]);
        C1 = [C1 Ct];
        D1 = [D1 Dt];
    end
    % % RR = computeFullGrid(C0,S0,D0,Awave,dt,BC) - C*S*D';
    % % [C1,S1,D1] = svd(RR);
    % % rr = sum(diag(S1) > tol) + 1; C1 = C1(:,1:rr); S1 = S1(1:rr,1:rr); D1 = D1(:,1:rr);
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
elseif normF < 0.5*tol %Need to decrease rank
    fprintf('-- Adaptive: Decreasing rank\n');
    [S_C,Sig,S_D] = svd(S);
    Sig = diag(Sig);
    rr = r;
    for i=r:-1:1
        if sum(Sig(i:r)) + normF >= 0.5*tol
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

