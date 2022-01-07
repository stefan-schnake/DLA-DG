
N = 64;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 2;

a = [1,1];

one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);

A = buildAdvectionMatrix(x,v,k,1);
Acell = buildAdvectionMatrixWithBlocks(x,v,k);

%[A,Lx] = buildAdvectionMatrixConst(x,v,k,a);
%A2 = kron(Lx{1},Lx{2}) + kron(Lx{2},Lx{1});
%A = full(A);
%A2 = full(A2);

%A1 = full(buildAdvectionMatrixConst(x,v,k,[1,0]));
%A2 = full(buildAdvectionMatrixConst(x,v,k,[0,1]));
%A3 = full(buildAdvectionMatrixConst(x,v,k,[1,1]));

%u0 = buildSeparableSource(x,v,k,@(x) (x > -5/6).*(x < -4/6),@(y) (y > -5/6).*(y < -4/6));
%u = u0;
%dt = 0.01;

%u = buildSeparableSource(x,v,k,@(x) sin(pi*x),@(y) cos(pi*y));
%u = buildSeparableSource(x,v,k,@(x) sin(pi*x),@(y) sin(pi*y));
u = rand(size(A,1),1);
U = convertVectoMat(x,v,k,u);
[X,S,V] = svd(U);

Fval = A*u;
Fmat = 0*U;
for i=1:size(Acell,1)
    Fmat = Fmat + Acell{i,1}*U*Acell{i,2}';
end
Fval2 = convertMattoVec(x,v,k,Fmat);


%Fval = A*u;
%Fmat = convertVectoMat(x,v,k,Fval);
%Fval2 = convertMattoVec(x,v,k,Lx{1}*U*Lx{2} + localTranspose(x,v,k,Lx{1}*locTransU*Lx{2}));
%Fval2 = convertMattoVec(x,v,k,Lx{1}*U*Lx{2}' +  Lx{2}*U*Lx{1}');
%FF = Lx{1}*locTransU*Lx{2};

%FX = Lx{2}*X;
%FV = Lx{1}*V;
%FF = FX*S*FV';
%FX1 = FX(:,1);
%FV1 = FV(:,1);

%plotVec(x,v,k,u);
%for i=1:100
%    u = (speye(size(A))+0.5*dt*A)\(u-0.5*dt*A*u);
%    plotVec(x,v,k,u);
%    drawnow
%end

function locTransU = localTranspose(x,v,k,U)
num_x = numel(x)-1;
num_v = numel(v)-1;

locTransU = zeros(size(U));
%Compose local transpose of U
for i=1:num_x
   for j=1:num_v
       %Get (k+1)\times(k+1) S chunk we care about
       locTransU((i-1)*(k+1)+1:i*(k+1),(j-1)*(k+1)+1:j*(k+1)) = U((i-1)*(k+1)+1:i*(k+1),(j-1)*(k+1)+1:j*(k+1))';
   end
end

end