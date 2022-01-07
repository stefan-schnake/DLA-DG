function [C,S,D] = DLA2_CONS_HYB(x,v,k,C,S,D,p,dt,A,BC)
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

[C,S,D] = EulerUpdate(x,v,k,C,S,D,p,dt,A,BC(:,2));
end

function [C,S,D] = EulerUpdate(x,v,k,C,S,D,p,dt,A,BC)

options = optimoptions(@fsolve,'Display','iter-detailed');


%Get u and F(u)=Au+BC; u_t = -(Au-BC);
%u = convertMattoVec(x,v,k,C*S*D');
%F = -convertVectoMat(x,v,k,A*u+BC);
%CFD = C'*F*D;

F = @(C,S,D) -convertVectoMat(x,v,k,A*convertMattoVec(x,v,k,C*S*D'));

%Evolve K = C*(S*S')
updateC = @(C) dt*F(C,S,D)*D*S(p+1:end,:)' - dt*C*(C'*F(C,S,D)*D)*S(p+1:end,:)';
K = C(:,p+1:end)*(S(p+1:end,:)*S(p+1:end,:)');
K1 = K + updateC(C);
[C1,~] = qr([C(:,1:p),K1],0);
K2 = (3/4)*K + (1/4)*(K1 + updateC(C1));
[C2,~] = qr([C(:,1:p),K2],0);
K_  = (1/3)*K + (2/3)*(K2 + updateC(C2));
[C_,~] = qr([C(:,1:p),K_],0);
M = C_'*C;

%Update L=D*(S'*S);

updateD = @(D) dt*F(C,S,D)'*C*S(:,p+1:end) - dt*D*((C'*F(C,S,D)*D)'*S(:,p+1:end));
%L = D(:,p+1:end)*(S(:,p+1:end)'*S(:,p+1:end));
%L = L + dt*F'*C*S(:,p+1:end) - dt*D*(CFD'*S(:,p+1:end));
%[D_,~] = qr(L,0);
%D_ = [D(:,1:p) D_];
%N = D_'*D;

L = D(:,p+1:end)*(S(:,p+1:end)'*S(:,p+1:end));
L1 = L + updateD(D);
[D1,~] = qr([D(:,1:p),L1],0);
L2 = (3/4)*L + (1/4)*(L1 + updateD(D1));
[D2,~] = qr([D(:,1:p),L2],0);
L_ = (1/3)*L + (2/3)*(L2 + updateD(D2));
[D_,~] = qr([D(:,1:p),L_],0);
N = D_'*D;


%Update S
SS = M*S*N';
Svec = reshape(SS,[],1);
%Svec_ = gmres(@(g) g - dt*updateS(x,v,k,g,C_,D_,A),Svec,10,1e-10);
Svec_ = fsolve(@(g) g - dt*updateS(x,v,k,g,C_,D_,A)-Svec,zeros(numel(Svec),1),options);
%Svec_ = fsolve(@(g) g - 0.5*dt*updateS(x,v,k,g,C_,D_,A) - 0.5*dt*updateS(x,v,k,Svec,C_,D_,A)-Svec,zeros(numel(Svec),1),options);
S_ = reshape(Svec_,size(S));

S = S_;
C = C_;
D = D_;

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

