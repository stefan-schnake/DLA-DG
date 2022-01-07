%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Script for Dynamic Low Rank Approximation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Parameters

%Number of cells in x and v direction
N = 512;

%X and V domain 
xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);

%Polynomial degree
k = 0;

%Rank
r = 20;

%Plot and save results
plotbool = true;
savebool = false;

%%Which test to run:
% test = 1 is spinning box
% test = 2 is smooth solution
test = 2;

%Set CFL
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%Making dt and integer division of pi
dt = pi/(8*ceil(1/CFL));

%Final time
T = pi;

%%%Timestepping method
% 'FE','RK2','RK3'
tstep = 'FE';

%%%-------------------------------------------

%Clear figures
if plotbool
figure(5)
clf('reset')
end

%Create wavelet to real space operator
FMWT = OperatorTwoScale_wavelet2(k+1,log2(N));

Acell = buildAdvectionMatrixWithBlocks2(x,v,k);
Awave = cell(size(Acell));
for i=1:numel(Acell)
    Awave{i} = FMWT*Acell{i}*FMWT';
end

if test == 1
init = @(x,y)  (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
BCsoln = @(x,v,t) soln_func(x,v,t,init);
source = @(x,v,t) 0*x;
elseif test == 2
u_sht = @(x,y) x+y+x.*y.^2; %u_short
L_usht = @(x,y) -y.*(1+y.^2) + x.*(1+2*x.*y); %L_u
u_lng = @(x,y) sin(pi*x).*sin(pi*y);
L_ulng = @(x,y) -pi*y.*cos(pi*x).*sin(pi*y) + pi*x.*sin(pi*x).*cos(pi*y);
%u_lng = @(x,y) 0*x+1;
%L_ulng = @(x,y) 0*x;
%tmdp = @(t) 1/(t+1);
%dt_tmdp = @(t) -1/(t+1)^2;
tmdp = @(t) exp(-5*t);
dt_tmdp = @(t) -5*exp(-5*t); 
BCsoln = @(x,y,t) (1-tmdp(t))*u_lng(x,y)+tmdp(t)*u_sht(x,y);
init = @(x,y) BCsoln(x,y,0);
soln   = {@(x,t) sin(pi*x), @(t) 1-tmdp(t), @(y,t) sin(pi*y);...
          @(x,t) x, @(t) tmdp(t), @(y,t) 0*y+1; ...
          @(x,t) 0*x+1, @(t) tmdp(t), @(y,t) y;...
          @(x,t) x, @(t) tmdp(t), @(y,t) y.^2};
source_vec = @(x,y,t) -dt_tmdp(t)*u_lng(x,y)+dt_tmdp(t)*u_sht(x,y)+(1-tmdp(t))*L_ulng(x,y)+tmdp(t)*L_usht(x,y);
%%%Coding negative of source!!!!
source = {@(x,t) sin(pi*x),    @(t) dt_tmdp(t) , @(y,t) sin(pi*y);...
          @(x,t) x,            @(t) -dt_tmdp(t), @(y,t) 0*y+1; ...
          @(x,t) 0*x+1,        @(t) -dt_tmdp(t), @(y,t) y;...
          @(x,t) x,            @(t) -dt_tmdp(t), @(y,t) y.^2;...
          @(x,t) -pi*cos(pi*x),@(t) -1+tmdp(t),  @(y,t) y.*sin(pi*y);...
          @(x,t) x.*sin(pi*x), @(t) -1+tmdp(t),  @(y,t) pi*cos(pi*y);...
          @(x,t) 0*x-1,        @(t) -tmdp(t),    @(y,t) y.*(1+y.^2);...
          @(x,t) x,            @(t) -tmdp(t),    @(y,t) 0*y+1;...
          @(x,t) 2*x.^2,       @(t) -tmdp(t),    @(y,t) y};
elseif test == 3
u_lng = @(x,y) x+y+x.*y.^2 +... x.^3.*y.^3 + x.^4.*y.^4 + ...
    sin(pi*(x-1/4)).*sin(pi*(y-1/4)) + cos(pi*(x-1/4)).*cos(pi*(y-1/4));
L_ulng = @(x,y) -y.*(1+y.^2) + x.*(1+2*x.*y) + ... -y.*3.*x.^2.*y.^3 + x.*3*x.^3.*y.^2 - y.*4.*x.^3.*y.^4 + x.*4.*x.^4.*y.^3;% + ...
    -y.*pi.*cos(pi*(x-1/4)).*sin(pi*(y-1/4)) + x.*pi.*sin(pi*(x-1/4)).*cos(pi*(y-1/4)) + ...
     y.*pi.*sin(pi*(x-1/4)).*cos(pi*(y-1/4)) - x.*pi.*cos(pi*(x-1/4)).*sin(pi*(y-1/4));
u_sht = @(x,y) sin(pi*x).*sin(pi*y);
L_usht = @(x,y) -pi*y.*cos(pi*x).*sin(pi*y) + pi*x.*sin(pi*x).*cos(pi*y);
%u_lng = @(x,y) 0*x+1;
%L_ulng = @(x,y) 0*x;
%tmdp = @(t) 1/(t+1);
%dt_tmdp = @(t) -1/(t+1)^2;
tmdp = @(t) exp(-10*t);
dt_tmdp = @(t) -10*exp(-10*t); 
soln = @(x,y,t) (1-tmdp(t))*u_lng(x,y)+tmdp(t)*u_sht(x,y);
init = @(x,y) soln(x,y,0);
source = @(x,y,t) -dt_tmdp(t)*u_lng(x,y)+dt_tmdp(t)*u_sht(x,y)+(1-tmdp(t))*L_ulng(x,y)+tmdp(t)*L_usht(x,y);
else
disp("Please select a test");
return 
end

%Build initial condition
u0 = buildNonSeparableSource(x,v,k,init);
u = u0;  

%Used to check conservation
one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);
cons = one'*u0;

%Take SVD of initial condition
matu = convertVectoMat(x,v,k,u);
[U,Sig,V] = svd(matu);
C0 = FMWT*U;
S0 = Sig;
D0 = FMWT*V;
R = [r]; 
C = C0(:,1:r);
S = S0(1:r,1:r);
D = D0(:,1:r);

U = C*S*D';
UU = U;
uu = u;
i = 0;
t = 0;
while (t+dt <= T+1e-9) 
i = i+1;
t0 = t;
t = t + dt;
lastplot = 0;
if t+dt > T+1e-9; lastplot=1; end
fprintf("-------------------------------------------------------\n");
fprintf("i=%d; t = %f\n",i,dt*i);
if test == 1
    BC.use = 0;    
else
    BC.use = 1;
    BC.cell1 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t-dt));buildSeparableSourceMat(x,v,k,t-dt,source)];
    BC.cell2 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t-0.5*dt));buildSeparableSourceMat(x,v,k,t-0.5*dt,source)];
    BC.cell3 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t));buildSeparableSourceMat(x,v,k,t,source)];
    %Convert BCs to Wavespace
    for l=1:size(BC.cell1,1)
       BC.cell1{l,1} = FMWT*BC.cell1{l,1}; 
       BC.cell1{l,3} = FMWT*BC.cell1{l,3}; 
       BC.cell2{l,1} = FMWT*BC.cell2{l,1}; 
       BC.cell2{l,3} = FMWT*BC.cell2{l,3}; 
       BC.cell3{l,1} = FMWT*BC.cell3{l,1}; 
       BC.cell3{l,3} = FMWT*BC.cell3{l,3}; 
    end
end
 

if strcmp(tstep,'FE')
    updateDLA = @(C,S,D) DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
elseif strcmp(tstep,'RK2')
    updateDLA = @(C,S,D) DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC);
elseif strcmp(tstep,'RK3')
    updateDLA = @(C,S,D) DLA4_HB_SSP_RK3(x,v,k,C,S,D,dt,Awave,BC);
end

%Update DLA
[C,S,D] = updateDLA(C,S,D);
    
R = [R size(S,1)];
r = size(S,1);
fprintf("r = %d\n",r);

%% Full-grid run
fullgrid = 0;
if fullgrid
    if BC.use
        U1 = UU - dt*applyMatA(UU,Awave);
        for l=1:size(BC.cell1,1)
            U1 = U1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
        %UU = (1/2)*UU + (1/2)*(U1-dt*applyMatA(U1,Awave)-dt*BC.vec(:,3));
        UU = U1;
    else
        U1 = UU - dt*applyMatA(UU,Awave);
        %UU = (1/2)*UU + (1/2)*(U1-dt*applyMatA(U1,Awave));
        UU = U1;
    end

    if BC.use
       LTEU1 = U - dt*applyMatA(U,Awave) - dt*BC.cell1{1,1}*BC.cell1{1,2}*BC.cell1{1,3}';
       for l=1:size(BC.cell1,1)
            LTEU1 = LTEU1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
       %LTEU = (1/2)*U + (1/2)*(LTEU1-dt*applyMatA(LTEU1,Awave)-dt*BC.vec(:,3));
       LTEU = LTEU1;
    else
       LTEU1 = U - dt*applyMatA(U,Awave);
       %LTEU = (1/2)*U + (1/2)*(LTEU1-dt*applyMatA(LTEU1,Awave));
       LTEU = LTEU1;
    end
    U = C*S*D';
    LTE = norm(LTEU-U,'fro');

end



%% Error Calc and plotting
%fprintf('-- Norm of BE: %e\n',norm(uu));
%fprintf('- Errors:\n');
%fprintf('-- Error of DLA   and L2 Proj: %e\n',norm(u-u_sol))
%fprintf('-- Error of BE    and L2 Proj: %e\n',norm(uu-u_sol))
%fprintf('-- Error of BE_lr and L2 Proj: %e\n',norm(uu_lr-u_sol))
if plotbool && ( mod(i,ceil(T/(10*dt))) == 0 || i == 1 || lastplot)
    u = convertMattoVec(x,v,k,FMWT'*C*S*D'*FMWT);
    figure(5)
    plotVec(x,v,k,u,@(x,y) BCsoln(x,y,t));
    sgtitle('DLA Solution')
 
    drawnow
    aa = 2;
end

end

%Save plots
if savebool
    fold = "figures/"+datestr(now,'mm-dd');
    
    if ~isfolder(fold)
        mkdir(fold);
    end
    cd(fold)
    figure(5)
    saveas(gcf,"N="+num2str(N)+"-k="+num2str(k)+"-test="+num2str(test)+"-tol="+num2str(adapt_tol)+"_DLA.eps",'epsc');
    %figure(6)
    %saveas(gcf,"test="+num2str(test)+"-tol="+num2str(adapt_tol)+"_FUL.eps",'epsc');
    figure(8)
    saveas(gcf,"N="+num2str(N)+"-k="+num2str(k)+"-test="+num2str(test)+"-tol="+num2str(adapt_tol)+"_MET.eps",'epsc');
    cd ..
    cd ..  
    
end

function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end


