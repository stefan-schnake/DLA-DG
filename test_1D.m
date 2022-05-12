
N = 32;
x = [-1:2/N:1];
v = [-1:2/N:1];

k = 0;
a_vec = [1,1];

Acell = buildConstantAdvectionMatrixWithBlocks(x,v,k,a_vec);

init = @(x) (x > 0).*(x < 1/4);

u = buildSeparableSourceX(x,k,init);
u0 = u;

dt = (x(2)-x(1));
T = 1;

Nt = ceil(T/dt);
for it = 1:Nt
    %u = u - dt*Acell{1,1}*u;
    u = (speye(size(Acell{1,1})) + dt*Acell{1,1})\u;
end