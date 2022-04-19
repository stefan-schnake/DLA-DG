N = 256;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 1;

per = reshape(convertVectoMat(x,v,k,1:(numel(x)-1)*(numel(v)-1)*(k+1)^2),[],1);
iper(per) = 1:numel(per);

box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
theta = 0.719747; %%N = 256
init = @(x,y) soln_func(x,y,theta,box);

%opt = @(theta) optfunc(x,v,k,per,theta);

%opts = optimset('Display','iter');
%[theta,fx] = fminbnd(opt,0,pi/4,opts); 

U.vec = buildNonSeparableSource(x,v,k,init);
LU.vec = slopeLimiter2D(x,v,k,U.vec,1);

U = vec2mat(U,x,v,k,per);
LU = vec2mat(LU,x,v,k,per);

LU = lrSVD(LU);
r = 250;
LU.S(r:end,r:end) = zeros(size(LU.S(r:end,r:end)));
LU = getMat(LU);
LU = mat2vec(LU,iper);

% 
% 
% U.S = svd(U.mat);
% LU.S = svd(LU.mat);
% semilogy(1:numel(U.S),U.S,1:numel(U.S),LU.S); legend({'L2 Proj','Limited'})

%[LU1D.C,LU1D.S,LU1D.D] = svd(U.mat);
% U = lrSVD(U);
% LU1D = U;
% for i=1:r
%     LU1D.C(:,i) = slopeLimiter(x,k,LU1D.C(:,i),1,1);
%     LU1D.D(:,i) = slopeLimiter(v,k,LU1D.D(:,i),1,1);
% end
% LU1D = getMat(LU1D);
% LU1D = mat2vec(LU1D,iper);




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

function U = lrSVD(U,tol)
    if nargin < 2
        tol = 1e-14;
    end
    [U.C,U.S,U.D] = svd(U.mat);
    r = sum(diag(U.S) > tol);
    U.C = U.C(:,1:r);
    U.S = U.S(1:r,1:r);
    U.D = U.D(:,1:r);
end

function U = getMat(U)
    U.mat = U.C*U.S*U.D';
end

function U = mat2vec(U,iper)
    U.vec = reshape(U.mat,[],1);
    U.vec = U.vec(iper);
end

function U = vec2mat(U,x,v,k,per)
    U.tmp = U.vec(per);
    U.mat = reshape(U.tmp,(numel(x)-1)*(k+1),(numel(v)-1)*(k+1));
    U.tmp = [];
end