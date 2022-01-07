function [C_,S_,D_] = DLA2_CONS_BE2(x,v,k,C,S,D,p,dt,A,BC)
%%DLA with SSP_RK3 timestepping; mass conservative update

% [C_,S_,D_] = EulerUpdate(x,v,k,C ,S ,D ,p,dt,A,BC(:,1));
% C1 = C_;
% S1 = S_;
% D1 = D_;
% [C_,S_,D_] = EulerUpdate(x,v,k,C1,S1,D1,p,dt,A,BC(:,3));
% C2 = (3/4)*C + (1/4)*C_;
% S2 = (3/4)*S + (1/4)*S_;
% D2 = (3/4)*D + (1/4)*D_;
% [C_,S_,D_] = EulerUpdate(x,v,k,C2,S2,D2,p,dt,A,BC(:,2));
% C = (1/3)*C + (2/3)*C_;
% S = (1/3)*S + (2/3)*S_;
% D = (1/3)*D + (2/3)*D_;

[C_,S_,D_] = EulerUpdate(x,v,k,C,S,D,p,dt,A,BC(:,2));
end

function [C_,S_,D_] = EulerUpdate(x,v,k,C,S,D,p,dt,A,BC)

%Get u and F(u)=Au+BC; u_t = -(Au-BC);
%u = convertMattoVec(x,v,k,C*S*D');
%F = -convertVectoMat(x,v,k,A*u+BC);
%CFD = C'*F*D;

options = optimoptions(@fsolve,'Display','iter-detailed','FunctionTolerance',1e-12,'OptimalityTolerance',1e-12);


%Evolve C(S)
n = numel(C)-size(C,1);
CSvec = [reshape(C(:,p+1:end),[],1);reshape(S,[],1)];
%Kvec_ = gmres(@(g) g - dt*updateK(x,v,k,g,C,S,D,p,A),Kvec,10,1e-10,numel(Kvec)/10,[],[],Kvec);
CSvec_ = fsolve(@(g) g - dt*updateK2(x,v,k,g,n,C,D,p,A)-CSvec,CSvec,options);
%Kvec_ = fsolve(@(g) g - 0.5*dt*updateK(x,v,k,g,C,S,D,p,A) - 0.5*dt*updateK(x,v,k,Kvec,C,S,D,p,A)-Kvec,zeros(numel(Kvec),1),options);
C_ = reshape(CSvec_(1:n),size(C,1),[]);
C_ = modGS([C(:,1:p),C_]);
M = C_'*C;

%Evolve L=D*(S'*S);
n = numel(D)-size(D,1);
DSvec = [reshape(D(:,p+1:end),[],1);reshape(S,[],1)];
%Lvec_ = gmres(@(g) g - dt*updateL(x,v,k,g,C,S,D,p,A),Lvec,10,1e-10);
DSvec_ = fsolve(@(g) g - dt*updateL2(x,v,k,g,n,C,D,p,A)-DSvec,DSvec,options);
%Lvec_ = fsolve(@(g) g - 0.5*dt*updateL(x,v,k,g,C,S,D,p,A) - 0.5*dt*updateL(x,v,k,Lvec,C,S,D,p,A)-Lvec,zeros(numel(Lvec),1),options);
D_ = reshape(DSvec_(1:n),size(D,1),[]);
D_ = modGS([D(:,1:p),D_]);
N = D_'*D;


%Update S
SS = M*S*N';
Svec = reshape(SS,[],1);
%Svec_ = gmres(@(g) g - dt*updateS(x,v,k,g,C_,D_,A),Svec,10,1e-10);
Svec_ = fsolve(@(g) g - dt*updateS(x,v,k,g,C_,D_,A)-Svec,Svec,options);
%Svec_ = fsolve(@(g) g - 0.5*dt*updateS(x,v,k,g,C_,D_,A) - 0.5*dt*updateS(x,v,k,Svec,C_,D_,A)-Svec,zeros(numel(Svec),1),options);
S_ = reshape(Svec_,size(S));

% F2 = @(C,S,D) -convertVectoMat(x,v,k,A*convertMattoVec(x,v,k,C*S*D'));
% updateS2 = @(S) dt*C_'*F2(C_,S,D_)*D_;
% S = M*S*N';
% S1 = S + updateS2(S);
% S2 = (3/4)*S + (1/4)*(S1 + updateS2(S1));
% S_ = (1/3)*S + (2/3)*(S2 + updateS2(S2));

end

function Kvec = updateK(x,v,k,Kvec,C,S,D,p,A)
%Convert to matrix
[row,~] = size(C);
K = reshape(Kvec,row,[]);

%Get C factors
%[C,~] = qr([C(:,1:p),K],0);
C = modGS([C(:,1:p),K]);
u = convertMattoVec(x,v,k,C*S*D');
F = -convertVectoMat(x,v,k,A*u);
K = F*D*S(p+1:end,:)' - C*(C'*F*D)*S(p+1:end,:)';

%Convert to vector
Kvec = reshape(K,[],1);

end

function Datvec = updateK2(x,v,k,Datvec,n,C_c,D,p,A)
%Convert to matrix
Cvec = Datvec(1:n);
Svec = Datvec(n+1:end);
[row,col] = size(C_c);
C = reshape(Cvec,row,[]);
C = [C_c(:,1:p) C];
S = reshape(Svec,col,[]);

u = convertMattoVec(x,v,k,C*S*D');
F = -convertVectoMat(x,v,k,A*u);

%Action of K where K=CSS';
K = F*D*S(p+1:end,:)' - C(:,1)*(C(:,1)'*F*D)*S(p+1:end,:)' + C(:,p+1:end)*S(p+1:end,:)*(C(:,p+1:end)'*F*D)';
%Update of S
S = C'*F*D;
%Get C
R = S(:,p+1:end)'*S(:,p+1:end);
C = K/R;

Datvec = [reshape(C,[],1);reshape(S,[],1)];
end

function Lvec = updateL(x,v,k,Lvec,C,S,D,p,A)
%Convert to matrix
[row,~] = size(D);
L = reshape(Lvec,row,[]);

%Get D factors
%[D,~] = qr([D(:,1:p),L],0);
D = modGS([D(:,1:p),L]);
u = convertMattoVec(x,v,k,C*S*D');
F = -convertVectoMat(x,v,k,A*u);
L = F'*C*S(:,p+1:end) - D*((C'*F*D)'*S(:,p+1:end));

%Convert to vector
Lvec = reshape(L,[],1);

end

function Datvec = updateL2(x,v,k,Datvec,n,C,D_c,p,A)
%Convert to matrix
Dvec = Datvec(1:n);
Svec = Datvec(n+1:end);
[row,col] = size(D_c);
D = reshape(Dvec,row,[]);
D = [D_c(:,1:p) D];
S = reshape(Svec,col,[]);

u = convertMattoVec(x,v,k,C*S*D');
F = -convertVectoMat(x,v,k,A*u);

%Action ofL=DS'S
L = F'*C*S(:,p+1:end) - D(:,1)*((C'*F*D(:,1))'*S(:,p+1:end)) + D(:,p+1:end)*S(:,p+1:end)'*((C'*F*D(:,p+1:end)));
%Update of S
S = C'*F*D;
%Get C
R = S(:,p+1:end)'*S(:,p+1:end);
D = L/R;

%Convert to vector
Datvec = [reshape(D,[],1);reshape(S,[],1)];
end

function Svec = updateS(x,v,k,Svec,C,D,A)
%Convert to matrix
[~,col] = size(C);
S = reshape(Svec,[],col);

u = convertMattoVec(x,v,k,C*S*D');
F = -convertVectoMat(x,v,k,A*u);
S = C'*F*D;

%Convert to vector
Svec = reshape(S,[],1);

end

function Z = modGS(A)
%Modified gram schmidt

Z = A;
[~,n] = size(A);
Z(:,1) = Z(:,1)/norm(Z(:,1),2);
for i=2:n
    for j=1:n-1
        Z(:,i) = Z(:,i) - (Z(:,j)'*Z(:,i))*Z(:,j);
    end
    Z(:,i) = Z(:,i)/norm(Z(:,1),2);
end


end

