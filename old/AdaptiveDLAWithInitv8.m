function [C,S,D] = AdaptiveDLAWithInitv8(x,v,k,C,S,D,t,dt,updateDLA,tol,Awave,BC,FMWT)
%Adaptive algorithm for DLA update

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
[C,S,D] = updateDLA(C,S,D);
%Compute SVD of S
[S_x,Sig,S_v] = svd(S);
Sig = diag(Sig);
fprintf('-- Adaptive: Smallest singular value is %e\n',Sig(end));
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
    Chalf = C;
    Dhalf = D;
    C = C_old; D = D_old; S = S_old;
    FF = FMWT*convertVectoMat(x,v,k,BC(:,1))*FMWT';
    [C1,S1,D1] = svds(@(x,trans) computeResidual2(C,S,D,Awave,FF,x,trans,Chalf,Dhalf),[1,1]*size(C,1),10);
    %Prune vectors with lower tolerance
    rr = sum(diag(S1) > tol);
    C1 = C1(:,1:rr); S1 = S1(1:rr,1:rr); D1 = D1(:,1:rr);
    %Rfunc = @(X,trans) computeResidual(C,S,D,Awave,X,trans);
    %[C1,~,D1] = svt(Rfunc,'m',size(C,1),'n',size(C,1),'lambda',tol,'k',r+15);
    [C1,~] = qr([Chalf C1],0);
    [D1,~] = qr([Dhalf D1],0);
    S = (C1'*C)*S*(D'*D1);
    C = C1;
    D = D1;
    fprintf('-- Adaptive: Bumping up r to %d\n',size(C,2));
    %[C,S,D] = AdaptiveDLAWithInitv8(x,v,k,C,S,D,t,updateDLA,tol,Awave,BC,FMWT);
    %[C,S,D] = updateDLA(C,S,D,t);
    S = DLA_HB_FE_S_Only(x,v,k,C,S,D,dt,Awave,BC,FMWT);
    [S_x,Sig,S_v] = svd(S);
    Sig = diag(Sig);
    fprintf('-- Adaptive: Smallest singular value is %e\n',Sig(end));
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
    end
    
end

end

