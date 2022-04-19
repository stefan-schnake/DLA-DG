N = 1024;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 1;

per = reshape(convertVectoMat(x,v,k,1:(numel(x)-1)*(numel(v)-1)*(k+1)^2),[],1);
iper(per) = 1:numel(per);

box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
theta = 0.719747; %%N = 256
init = @(x,y) soln_func(x,y,theta,box);

opt = @(theta) optfunc(x,v,k,per,theta);

opts = optimset('Display','iter');
[theta,fx] = fminbnd(opt,0,pi/4,opts); 

% u = buildNonSeparableSource(x,v,k,init);
% Lu = slopeLimiter2D(x,v,k,u,1);
% 
% U.mat = reshape(u(per),(numel(x)-1)*(k+1),[]);
% LU.mat = reshape(Lu(per),(numel(x)-1)*(k+1),[]);
% 
% U.S = svd(U.mat);
% LU.S = svd(LU.mat);
% semilogy(1:numel(U.S),U.S,1:numel(U.S),LU.S); legend({'L2 Proj','Limited'})



function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end

function z = optfunc(x,v,k,per,theta)
    box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
    init = @(x,y) soln_func(x,y,theta,box);
    u = buildNonSeparableSource(x,v,k,init);
    U = reshape(u(per),(numel(x)-1)*(k+1),[]);
    S = svd(U);
    
    z = -sum(S > 1e-14);
end