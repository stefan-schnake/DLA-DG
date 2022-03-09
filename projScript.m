N = 256;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 1;

%% Create projection
box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
soln = @(x,y) soln_func(x,y,pi/4,box);

proj.vec = buildNonSeparableSource(x,v,k,soln);
proj.mat = convertVectoMat(x,v,k,proj.vec);
proj.S = svd(proj.mat);
lim.vec = slopeLimiter2D(x,v,k,proj.vec,5,1);
lim.mat = convertVectoMat(x,v,k,lim.vec);
lim.S = svd(lim.mat);




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

