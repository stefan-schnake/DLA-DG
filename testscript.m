N = 256;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 1;

LF = 0.1;

T = pi/4;

%Create permuation vectors
u0.real = buildSeparableSource(x,v,k,@(x) (abs(x) < 1/2),@(v) (abs(v) < 1/2) );

perm = convertVectoMat(x,v,k,1:numel(u0.real));
perm = perm(:);
iperm = zeros(size(perm));
iperm(perm) = 1:numel(perm);

u.real = u0.real;
u.perm = u0.real(perm);

CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%dt = 1/2*CFL;
%dt = pi/(8*ceil(1/CFL));
dt = 4.382802251101832e-04;

Acell = buildAdvectionMatrixWithBlocks2(x,v,k,LF);
A = sparseKronAdd(Acell);
A = A(iperm,iperm);

t = 0;i = 0;
while (t+dt <= T+1e-9) 
%while 0
    i = i + 1;
    t = t + dt;
    if mod(i,100) == 0
        fprintf('i = %5d: t = %.4f |u|_2 = %f\n',i,t,norm(u.real));
    end
    %u.perm = SSP_RK3(A,dt,u.perm);
    u.real = SSP_RK3_lim(x,v,k,A,dt,u.real);
    %u.real = SSP_RK3(A,dt,u.real);
end
%u.real = u.perm(iperm);

plotVec(x,v,k,u.real);
u.perm = u.real(perm);
[u.X,u.S,u.V] = svd(reshape(u.perm,N*(k+1),[]));
u.mat = u.X*u.S*u.V';
u.S = diag(u.S);


function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end

function z = SSP_RK3(A,dt,z)
    z1 = z - dt*A*z;
    z2 = (3/4)*z + (1/4)*(z1 - dt*A*z1);
    z  = (1/3)*z + (2/3)*(z2 - dt*A*z2);
end

function z = SSP_RK3_lim(x,v,k,A,dt,z)
    z1 = z - dt*A*z;
    z1 = slopeLimiter2D(x,v,k,z1,1,1);
    z2 = (3/4)*z + (1/4)*(z1 - dt*A*z1);
    z2 = slopeLimiter2D(x,v,k,z2,1,1);
    z  = (1/3)*z + (2/3)*(z2 - dt*A*z2);
    z = slopeLimiter2D(x,v,k,z,1,1);
end

