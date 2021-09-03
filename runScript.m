
N = 32;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 2;

r = 15;

T = pi;

plotbool = true;

adapt = false;
adapt_tol = 1e-4;

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%CFL = 1/N^2;
%dt = 0.8*CFL;
%CFL = (5/4)*(1/73536);
%dt = 0.9*CFL;
dt = pi/(8*ceil(1/CFL));
%dt = 0.05;

%T = 10*dt;
T = 0.5;

%%%-------------------------------------------

L = buildLDGMatrixAlt(x,v,k);
LDG = L'*L;
Jmp = buildJumpMatrix(x,v,k,1);
Adv = buildAdvectionMatrix(x,v,k,true);
Adv_back = buildBackAdvectionMatrix(x,v,k,true);

%Use diffusion or advection
A = Adv;
%A2 = Adv_back;
%A = LDG+N*Jmp;
%A2 = LDG+N*Jmp;
%A = LDG;
%A2 = LDG;

init = @(x,y)  (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
%init = @(x,y) 16*(x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2).*(1/4-x.^2).*(1/4-y.^2);
%init = @(x,y) (1-x.^2).*(1-y.^2);
%init = @(x,y) sin(pi*x).*sin(pi*y);
%init = @(x,y) cos(pi*x).*cos(pi*y);
%init = @(x,y) exp(-10*(x.^2+y.^2));
%init = @(x,y) x.*y;
%init = @(x,y) (x > -.25).*(x < .25).*(y > .25).*(y < .75).*x;
%soln = @(x,y,t) init(x+t*y,y-t*x);
soln = @(x,v,t) soln_func(x,v,t,init);
%soln = @(x,v,t) 0*x;
%soln = @(x,y,t) cos(pi*x).*cos(pi*y)*exp(-2*pi^2*t);

%u0 = buildSeparableSource(x,v,k,@(x) cos(pi*x),@(y) cos(pi*y));
%u0 = buildSeparableSource(x,v,k,@(x) exp(x),@(y) cos(y));
%u0 = buildSeparableSource(x,v,k,@(x) (x > -1/2).*(x < 1/2),@(y) (y > -1/2).*(y < 1/2));
u0 = buildNonSeparableSource(x,v,k,init);
u = u0;  

one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);
cons = one'*u0;
    
%Run a few timesteps to smooth the solution out
%for i=1:5
%u = u + dt*LDG*u;
%end

uu = u0;

matu = convertVectoMat(x,v,k,u);
[U,Sig,V] = svd(matu);

Sig = diag(Sig);
if adapt
    r = sum(Sig > adapt_tol)+1;
    if r > size(U,1); r = size(U,1); end
    fprintf('-- Initial Adaptive r: r = %d\n',r);
    R = [r];
else
   R = r; 
end
Sig = diag(Sig);

%Get rank5 approximation
C = U(:,1:r);
D = V(:,1:r);
S = Sig(1:r,1:r);

rank1 = log10(Sig(1))-log10(Sig(2));


u = convertMattoVec(x,v,k,C*S*D');
i = 0;
while (dt*(i+1) <= T+1e-9) 
i = i+1;
fprintf("-------------------------------------------------------\n");
fprintf("i=%d; t = %f\n",i,dt*i);
BC = [buildDirichletBC(x,v,k,@(x,y) soln(x,y,(i-1)*dt)), ...
      buildDirichletBC(x,v,k,@(x,y) soln(x,y,(i-0.5)*dt)), ...
      buildDirichletBC(x,v,k,@(x,y) soln(x,y,i*dt))];
%BC = buildDirichletBC(x,v,k,@(x,y) soln(x,y,i*dt));
      
%BC_old = buildDirichletBC(x,v,k,@(x,y) soln(x,y,(i-1)*dt));
%BC_hlf = buildDirichletBC(x,v,k,@(x,y) soln(x,y,(i-1/2)*dt));

%updateDLA = @(C,S,D) DLA(x,v,k,C,S,D,dt,A,A2,BC);
%updateDLA = @(C,S,D) DLA_CN(x,v,k,C,S,D,dt,A,A2,BC);
%updateDLA = @(C,S,D) DLA2_BE(x,v,k,C,S,D,dt,A,BC);
updateDLA = @(C,S,D) DLA2_SSP(x,v,k,C,S,D,dt,A,BC);

if adapt
    [C,S,D] = AdaptiveDLA(C,S,D,updateDLA,adapt_tol);
    R = [R size(S,1)];
else
    %[C,S,D] = DLA(x,v,k,C,S,D,dt,A,Adv_back);
    [C,S,D] = updateDLA(C,S,D);
end

sing = svd(S);
rank1 = [rank1 log10(sing(1))-log10(sing(2))];
fprintf('Min singular value of S: %e\n',min(sing));

%Get new u
u = convertMattoVec(x,v,k,C*S*D');


%% Traditonal Backward Euler
%uu = pcg(speye(size(LDG))+dt*LDG,uu,1e-13);
%%%BE
%uu = (speye(size(LDG))+dt*A)\(uu);
%%%CN
%uu = (speye(size(LDG))+dt/2*A)\( (speye(size(LDG))-dt/2*A)*uu - 0.5*dt*(BC(:,1)+BC(:,3)));
%uu = (speye(size(LDG))+dt/2*A)\( (speye(size(LDG))-dt/2*A)*uu );
%uu = uu - dt*A*uu;
%uu = expm(-A*(i*dt))*u0;
%%%SSP-RK3
u1 = uu - dt*A*uu - dt*BC(:,1);
u2 = (3/4)*uu + (1/4)*(u1-dt*A*u1-dt*BC(:,3));
uu = (1/3)*uu + (2/3)*(u2-dt*A*u2-dt*BC(:,2));
[UU,SS,VV] = svd(convertVectoMat(x,v,k,uu));
UU_lr = UU(:,1:R(end));
SS_lr = SS(1:R(end),1:R(end));
VV_lr = VV(:,1:R(end));
uu_lr = convertMattoVec(x,v,k,UU_lr*SS_lr*VV_lr');

%% Error Calc and plotting
u_sol = buildNonSeparableSource(x,v,k,@(x,y) soln(x,y,i*dt));
%fprintf('-- Norm of BE: %e\n',norm(uu));
%fprintf('- Errors:\n');
%fprintf('-- Error of DLA   and L2 Proj: %e\n',norm(u-u_sol))
%fprintf('-- Error of BE    and L2 Proj: %e\n',norm(uu-u_sol))
%fprintf('-- Error of BE_lr and L2 Proj: %e\n',norm(uu_lr-u_sol))
if plotbool
    figure(5)
    plotVec(x,v,k,u,@(x,y) soln(x,y,i*dt));
    %plotVec(x,v,k,u);
    sgtitle('DLA Solution')
    %imagesc(flipud(reshape(u,N,[])));colorbar;
    figure(6)
    %imagesc(flipud(reshape(uu,N,[])));colorbar;
    %plotVec(x,v,k,uu);
    plotVec(x,v,k,uu,@(x,y) soln(x,y,i*dt));
    sgtitle('Full-Grid Solution')
    drawnow
end
end

if adapt
    figure(7)
    plot((1:numel(R))-1,R)
    title('Adaptive Rank plot');
    xlabel('iteration')
    ylabel('rank')
end

function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end


