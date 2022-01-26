N = 64;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 1;

p.n  = @(x) 1*(x<-0.5) + 1/8*(x >= -0.5).*(x <= 0.5) + 1*(x>0.5);
p.u  = @(x) 0;
p.th = @(x) 1*(x<-0.5) + 5/4*(x >= -0.5).*(x <= 0.5) + 1*(x>0.5);

maxwell = @(x,v) p.n(x)./sqrt(2*pi*p.th(x)).*exp(-(v-p.u(x)).^2./(2*p.th(x)));

u0 = buildNonSeparableSource(x,v,k,maxwell);
U0.mat = convertVectoMat(x,v,k,u0);

[U0.C,U0.S,U0.V] = svd(U0.mat);
