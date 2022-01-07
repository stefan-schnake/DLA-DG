function [C,S,D,t] = AdaptiveDLAWithInitv4(x,v,k,C,S,D,t0,t1,updateDLA,tol,i,C0,S0,D0,Acell,RHS,FMWT)
%Adaptive algorithm for DLA update

r = size(S,1);
fprintf('-- Adaptive: r = %d\n',r);
C_old = C; S_old = S; D_old = D;
%Update
[C,S,D] = updateDLA(C,S,D,t1);
t = t1;
%[C,S,D] = DLA(x,v,k,C_old,S_old,D_old,dt,A,A2);
%Compute SVD of S
[S_x,Sig,S_v] = svd(S);
Sig = diag(Sig);
fprintf('-- Adaptive: Smallest singular value is %e\n',Sig(end));
if Sig(end) < tol %We don't need to bump up r
    rnew = sum(Sig > tol)+1; %See if we need to decrease r
    Sig = diag(Sig);
    if rnew ~= r
        r = rnew; 
        fprintf('-- Adaptive: New r = %d\n',r);
        S = Sig(1:r,1:r);
        C = C*S_x(:,1:r);
        D = D*S_v(:,1:r);
    end
elseif r < size(C,1) %Add and element and try again
    %if i <= 1
    %    C = C0(:,1:r+1);
    %    S = S0(1:r+1,1:r+1);
    %    D = D0(:,1:r+1);
    %    [C,S,D,t] = AdaptiveDLAWithInitv4(x,v,k,C,S,D,t,updateDLA,tol,i,C0,S0,D0,Acell,RHS,FMWT);
    %else 
        C = C_old; S = S_old; D = D_old;
        %Throw out solution and advance by full-rank PDE
        newdt = 1e-10;
        G = cell(size(Acell,1)+2,3);
        G{1,1} = C;
        G{1,2} = S;
        G{1,3} = D;
        for l=1:size(Acell,1)
            G{l+1,1} = Acell{l,1}*C;
            G{l+1,2} = -newdt*S;
            G{l+1,3} = Acell{l,2}*D;
        end
        %Do SVD truncation of full-rank RHS vector
        RHS3 = RHS(t0,0);
        [RHS_U,RHS_S,RHS_V] = svd(convertVectoMat(x,v,k,RHS3(:,1)));
        RHS_r = sum(diag(RHS_S) > eps);
        G{end,1} = FMWT*RHS_U(:,1:RHS_r);
        G{end,2} = -newdt*RHS_S(1:RHS_r,1:RHS_r);
        G{end,3} = FMWT*RHS_V(:,1:RHS_r);
        
        %Find r+1 DLA to U^n-dt*A*U^n-dt*F^n
        %[C,S,D] = adaptHierAdd(C,S,D,5,G);
        afun = @(x,transflag) myafun(G,transflag,x);
        [C,S,D] = svds(afun,[size(C,1) size(C,1)],r+5);
        
        %Recursive call
        [C,S,D,t] = AdaptiveDLAWithInitv4(x,v,k,C,S,D,t0+newdt,t1+newdt,updateDLA,tol,i,C0,S0,D0,Acell,RHS,FMWT);
    %end
end

end

function y = myafun(G,transflag,x)
    NN = size(G,1);
    y = zeros(size(x));
    if strcmp(transflag,'notransp')
        for l=1:NN
            y = y + G{l,1}*(G{l,2}*(G{l,3}'*x));
        end
    else
        for l=1:NN
            y = y + G{l,3}*(G{l,2}'*(G{l,1}'*x));
        end
    end
end

