N = 20;
k = 1;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);

Acell = buildIPDGMatrixWithBlocks(x,v,k,10);
A = sparseKronAdd(Acell);
L = buildLDGMatrix(x,v,k);
Jmp = buildJumpMatrix(x,v,k,1);
LDG = L'*L + Jmp;

f.real = buildSeparableSource(x,v,k,@(x) x.^2,@(y) 0*y+1);

perm = convertVectoMat(x,v,k,1:numel(f.real));
perm = perm(:);
iperm = zeros(size(perm));
iperm(perm) = 1:numel(perm);

f.perm = f.real(perm);


u.perm = A\f.perm;
u.real = u.perm(iperm);

u_ldg = LDG\f.real;

plotVec(x,v,k,u.real);