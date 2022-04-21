N = 128;
k = 2;

x = -1:2/N:1;

init = @(x) x > 1/3;

u = buildSeparableSourceX(x,k,init);

Lu = slopeLimiter(x,k,u,1,1);
Wu = wenoLimiter(x,k,u,1);

fprintf("L^2 Error: Pu = %e, Lu = %e, Wu = %e\n",getL2Error(x,k,u,init),getL2Error(x,k,Lu,init),getL2Error(x,k,Wu,init));