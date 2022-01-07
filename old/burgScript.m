
N = 32;

xx = [-3,3];vv = [-3,3];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 2;

r = 5;

T = 1;

adapt = true;
adapt_tol = 1e-1;

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = 1/(2*k+1)*(1/N);
dt = 0.8*CFL;

%%%-------------------------------------------

%init = @(x,v) exp(-10*((x-3).^2+(v-3).^2));
init = @(x,v) (x+v > 0);
soln = @(x,v,t) (x+v > 0).*( (x+v) < (2*t) ).*(x+v)/(2*t) + (x+v >= (2*t)); 

u0 = buildNonSeparableSource(x,v,k,init);
uu = u0;

one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);

myx = buildSeparableSource(x,v,k,@(x) x,@(v) v);

%Get LA or u
[C0,S0,D0] = svd(convertVectoMat(x,v,k,u0)); 

Sig = diag(S0);
if adapt
    r = sum(Sig > adapt_tol)+1;
    if r > size(C0,1); r = size(C0,1); end
    fprintf('-- Initial Adaptive r: r = %d\n',r);
    R = [r];
else
    R = [r];
end

C = C0(:,1:r);
S = S0(1:r,1:r);
D = D0(:,1:r);

%evalScripts
F1 = @(g) evalBurgers2D(x,v,k,g);
F2 = @(g) evalBurgers2Dback(x,v,k,g);

updateDLA = @(C,S,D) DLAfunc(x,v,k,C,S,D,dt,F1,F2);

t = 0;
i = 0;
while t < T
fprintf("-------------------------------------------------------\n");
fprintf("i=%d; t = %f\n",i+1,dt*(i+1));
i = i + 1;
t = i*dt;

if adapt
    [C,S,D] = AdaptiveDLA(C,S,D,updateDLA,adapt_tol);
    R = [R size(S,1)];
else
    [C,S,D] = updateDLA(C,S,D);
end
u = convertMattoVec(x,v,k,C*S*D');

%% SSP-RK3
F = @(g) g - dt*evalBurgers2D(x,v,k,g);
f1 = F(uu);
f2 = 0.75*uu + 0.25*F(f1);
uu = (1/3)*uu + (2/3)*F(f2);
%uu = fsolve(@(g) g + dt*evalBurgers2D(x,v,k,g) - uu,uu);

%% Error Calc and plotting
fprintf('-- Norm of DLA : %e\n',norm(u));
figure(5)
plotVec(x,v,k,u,@(x,v) soln(x,v,t),1);
sgtitle("DLA Solution r = "+num2str(R(end)))
fprintf('-- Norm of Full: %e\n',norm(uu));
figure(6)
plotVec(x,v,k,uu,@(x,v) soln(x,v,t),1);
sgtitle('Full-grid Solution')
drawnow
end

if adapt
    figure(7)
    plot((1:numel(R))-1,R)
    title('Adaptive Rank plot');
    xlabel('iteration')
    ylabel('rank')
end


