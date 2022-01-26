N = 1024;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 0;

T = pi/4;

%Create permuation vectors
u0.real = buildSeparableSource(x,v,k,@(x) (abs(x) < 1/2),@(v) (abs(v) < 1/2) );

perm = convertVectoMat(x,v,k,1:numel(u0.real));
perm = perm(:);
iperm = zeros(size(perm));
iperm(perm) = 1:numel(perm);

u.perm = u0.real(perm);

CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%dt = 1/2*CFL;
%dt = pi/(8*ceil(1/CFL));
dt = 4.382802251101832e-04;

Acell = buildAdvectionMatrixWithBlocks2(x,v,k);
A = sparseKronAdd(Acell);

t = 0;i = 0;
while (t+dt <= T+1e-9) 
%while 0
    i = i + 1;
    t = t + dt;
    if mod(i,100) == 0
        fprintf('i = %5d\n',i);
    end
    u.perm = SSP_RK3(A,dt,u.perm);
end
u.real = u.perm(iperm);

plotVec(x,v,k,u.real);
[u.X,u.S,u.V] = svd(reshape(u.perm,N*(k+1),[]));
u.mat = u.X*u.S*u.V';
u.S = diag(u.S);


%% Create projection
box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
soln = @(x,y) soln_func(x,y,pi/4,box);

proj.vec = buildNonSeparableSource(x,v,k,soln);
proj.mat = convertVectoMat(x,v,k,proj.vec);
proj.S = svd(proj.mat);
%plotVec(x,v,k,proj.vec);



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

