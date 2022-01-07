
N = 64;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 2;

r = 3;

r_cut = 3;

plotbool = true;

adapt = true;
adapt_tol = 1e-9;

hier_adapt = false;
hier_tol = 0.5;

test = 2;

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%CFL = 1/N^2;
%dt = 0.8*CFL;
%CFL = (5/4)*(1/73536);
%dt = 0.9*CFL;
dt = pi/(8*ceil(1/CFL));
%dt = 0.05;

T = .5;

%%%-------------------------------------------

FMWT = OperatorTwoScale_wavelet2(k+1,log2(N));

init = @(x,y) 1 + x.*y + x.^2.*y.^2 + x.^3.*y.^3;
%u0 = buildSeparableSource(x,v,k,@(x) cos(pi*x),@(y) cos(pi*y));
%u0 = buildSeparableSource(x,v,k,@(x) exp(x),@(y) cos(y));
%u0 = buildSeparableSource(x,v,k,@(x) (x > -1/2).*(x < 1/2),@(y) (y > -1/2).*(y < 1/2));
u0 = buildNonSeparableSource(x,v,k,init);
u = u0;  

one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);

%Run a few timesteps to smooth the solution out
%for i=1:5
%u = u + dt*LDG*u;
%end

uu = u0;

matu = convertVectoMat(x,v,k,u);
[C0,Sig,D0] = svd(matu);
S0 = Sig;


C = C0(:,1:r);
S = S0(1:r,1:r);
D = D0(:,1:r);
tol = min(diag(S))/10;
%tol = 1e-6;

U = C*S*D';
C = FMWT*C;
D = FMWT*D;

[C1,S1,D1] = adaptHierAdd(x,v,k,C,S,D,1,convertMattoVec(x,v,k,matu),FMWT);
%[C1,S1,D1] = adaptHierAdd(x,v,k,C,S,D,1,convertMattoVec(x,v,k,U),FMWT);
UU = (FMWT'*C1)*S1*(FMWT'*D1)';

% options = optimoptions(@fmincon,'Display','iter-detailed');
% %options = optimoptions(options,'SpecifyObjectiveGradient',true,'CheckGradients',true);
% options = optimoptions(options,'MaxFunctionEvaluations',1e6);
% optionsGrad = optimoptions(options,'SpecifyObjectiveGradient',true);
% 
% maxC = max(sum(abs(C)));
% maxD = max(sum(abs(D)));
% 
% %[C_old,S_old,D_old] = adaptPlus1(C,S,D);
% %[Sc,S_sig,Sd] = svd(S_old);
% %S_sig(end,end) = tol;
% %C_ = C_old*Sc*sqrt(S_sig);
% %D_ = D_old*Sd*sqrt(S_sig);
% C_ = C0(:,1:r+1)*sqrt(S0(1:r+1,1:r+1));
% D_ = D0(:,1:r+1)*sqrt(S0(1:r+1,1:r+1));
% 
% for i=1
% Cvec = reshape(C_',[],1);
% objC = @(Cvec) objMinC(Cvec,D_,U);
% noncon = @(Cvec) objNonlconC(Cvec,D_,sqrt(tol),maxC);
% 
% %Cvec_ = fmincon(objC,rand(size(Cvec)),[],[],[],[],[],[],noncon,optionsGrad);
% Cvec_ = fmincon(objC,Cvec,[],[],[],[],[],[],noncon,optionsGrad);
% C_ = reshape(Cvec_,r+1,[])';
% 
% Dvec = reshape(D_,[],1);
% objD = @(Dvec) objMinD(Dvec,C_,U);
% noncon = @(Dvec) objNonlconD(Dvec,C_,sqrt(tol),maxD);
% 
% %Dvec_ = fmincon(objD,rand(size(Dvec)),[],[],[],[],[],[],noncon,options);
% Dvec_ = fmincon(objD,Dvec,[],[],[],[],[],[],noncon,optionsGrad);
% D_ = reshape(Dvec_,[],r+1);
% end
% 
% UU = C_*D_';
% [C1,S1a] = qr(C_,0);
% [D1,S1b] = qr(D_,0);
% S1 = S1a*S1b';



