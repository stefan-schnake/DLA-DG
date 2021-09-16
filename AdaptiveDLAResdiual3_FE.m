function [C,S,D] = AdaptiveDLAResdiual3_FE(x,v,k,C,S,D,dt,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update
tic
r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C0 = C; S0 = S; D0 = D;
%Update
[C,S,D] = DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%Compute SVD of S
[S_C,Sig,S_D] = svd(S);
fprintf('-- Adaptive: Smallest singular value is %e with tolerance %e\n',Sig(end),tol);
%rflag =  svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D),[1,1]*size(C,1),1);
%fprintf('-- Adaptive: Largest  singular value of residual is %e\n',rflag);
%Compute estimation to residual
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
    C1 = [];
    S1 = [];
    D1 = [];
    tolflag = true;
    while tolflag == true
        rr = rr + 5;
        if rr >= size(C,1)
            rr = size(C,1);
            tolflag = false;
        end
        if rr <= 5
            [Ct,St,Dt] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D),[1,1]*size(C,1),rr);
        else
            [Ct,St,Dt] = svds(@(x,trans) firstOrderResidual(C0,S0,D0,Awave,dt,BC,x,trans,C,D) ...
            - LReval(C1,S1,D1,x,trans),[1,1]*size(C,1),5,'largest','FailureTreatment','drop');
        end
        if St(end,end) < tol
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
    r = sum(diag(Sig) > tol) + 1;
    fprintf('-- Adaptive: Cutting down to %d vectors, tol = %e\n',r,Sig(end,end));
    C = C*S_C(:,1:r);
    S = Sig(1:r,1:r);
    D = D*S_D(:,1:r);
end

toc
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

function Y = LReval(C,S,D,X,trans)

if strcmp(trans,'notransp')
    Y = C*S*(D'*X);
else
    Y = D*S'*(C'*X);
end

end
