function [C,S,D,UU] = hierFunc(N,k,dtfrac,xx,vv,r,test,fullgrid)

x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%CFL = 1/N^2;
%dt = 0.8*CFL;
%CFL = (5/4)*(1/73536);
%dt = 0.9*CFL;
%dt = pi/(8*ceil(1/CFL));
%dt = 1/32*CFL;
dt = dtfrac*CFL;
%dt = 0.05;
%dt = 0.05;

moviebool = false;
plotbool  = false;
savebool  = false;

adapt  = false;
adapt_tol = dt;

T = 1;
%T = pi;
%T = 3.5;

%%%-------------------------------------------

%Clear figures
if plotbool
figure(5)
clf('reset')
figure(6)
clf('reset')
%figure(7)
%clf('reset')
figure(8)
clf('reset')
end

if moviebool
    myVideo = VideoWriter('movie.mp4','MPEG-4');
    myVideo.FrameRate = 1;
    open(myVideo);
end

%L = buildLDGMatrixAlt(x,v,k);
%LDG = L'*L;
%Jmp = buildJumpMatrix(x,v,k,1);
%Adv = buildAdvectionMatrix(x,v,k,true);
%Adv_back = buildBackAdvectionMatrix(x,v,k,true);
FMWT = OperatorTwoScale_wavelet2(k+1,log2(N));

Acell = buildAdvectionMatrixWithBlocks2(x,v,k);
Awave = cell(size(Acell));
for i=1:numel(Acell)
    Awave{i} = FMWT*Acell{i}*FMWT';
end

%Use diffusion or advection
%A = Adv;
%A2 = Adv_back;
%A = LDG+N*Jmp;
%A2 = LDG+N*Jmp;
%A = LDG;
%A2 = LDG;

if test == 1
init = @(x,y)  (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
%box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
%init = @(x,y) soln_func(x,y,pi/4,box);
%init = @(x,y) x.^2.*y + x;
%init = @(x,y) exp(10*(-(x-0.7).^2-(y-0.4)^2));
%init = @(x,y) 16*(x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2).*(1/4-x.^2).*(1/4-y.^2);
%init = @(x,y) (1-x.^2).*(1-y.^2);
%init = @(x,y) exp(1/4*(x-y));
%init = @(x,y) sin(pi*x).*sin(pi*y);
%init = @(x,y) cos(pi*x).*cos(pi*y);
%init = @(x,y) exp(-10*(x.^2+y.^2));
%init = @(x,y) x.*y;
%init = @(x,y) (x > -.25).*(x < .25).*(y > .25).*(y < .75).*x;
%soln = @(x,y,t) init(x+t*y,y-t*x);
BCsoln = @(x,v,t) soln_func(x,v,t,init);
source = @(x,v,t) 0*x;
%soln = @(x,v,t) 0*x;
%soln = @(x,y,t) cos(pi*x).*cos(pi*y)*exp(-2*pi^2*t);
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

%u0 = buildSeparableSource(x,v,k,@(x) cos(pi*x),@(y) cos(pi*y));
%u0 = buildSeparableSource(x,v,k,@(x) exp(x),@(y) cos(y));
%u0 = buildSeparableSource(x,v,k,@(x) (x > -1/2).*(x < 1/2),@(y) (y > -1/2).*(y < 1/2));
u0 = buildNonSeparableSource(x,v,k,init);
u = u0;  
one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);



matu = convertVectoMat(x,v,k,u);
[U,Sig,V] = svd(matu);
C0 = FMWT*U;
S0 = Sig;
D0 = FMWT*V;

Sig = diag(Sig);
if adapt
    r = sum(Sig > adapt_tol)+1;
    if r > size(U,1); r = size(U,1); end
    fprintf('-- Initial Adaptive r: r = %d\n',r);
    R = [r];
else
    R = [r]; 
end
%r = 5;
C = C0(:,1:r);
S = S0(1:r,1:r);
D = D0(:,1:r);

%Used as storage for multistep

U = C*S*D';
UU = U;
%u = convertMattoVec(x,v,k,FMWT'*C*S*D'*FMWT);
uu = u;
i = 0;
t = 0;
%while i < 1
while (t+dt <= T+1e-9) 
i = i+1;
t0 = t;
t = t + dt;
lastplot = 0;
if t+dt > T+1e-9; lastplot=1; end
if test == 1
    BC.use = 0;
    BC2.use = 0;
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
    
    BC2.use = 1;
    MM = zeros(size(BC.cell1{1,1},1),size(BC.cell1,1));
    NN = zeros(size(BC.cell1{1,3},1),size(BC.cell1,1));
    SS = zeros(size(BC.cell1,1));
    for l=1:size(BC.cell1,1)
       MM(:,l)  = BC.cell1{l,1}; 
       NN(:,l)  = BC.cell1{l,3}; 
       SS(l,l) = BC.cell1{l,2}; 
    end
    BC2.cell1 = {MM,SS,NN};
    for l=1:size(BC.cell1,1)
       MM(:,l)  = BC.cell2{l,1}; 
       NN(:,l)  = BC.cell2{l,3}; 
       SS(l,l) = BC.cell2{l,2}; 
    end
    BC2.cell2 = {MM,SS,NN};
    for l=1:size(BC.cell1,1)
       MM(:,l)  = BC.cell3{l,1}; 
       NN(:,l)  = BC.cell3{l,3}; 
       SS(l,l) = BC.cell3{l,2}; 
    end
    BC2.cell3 = {MM,SS,NN};
        
end
BC = BC2;
clear BC2

%updateDLA = @(C,S,D) DLA(x,v,k,C,S,D,dt,A,A2,BC);
%updateDLA = @(C,S,D) DLA_CN(x,v,k,C,S,D,dt,A,A2,BC);
%updateDLA = @(C,S,D) DLA2_HB_SSP(x,v,k,C,S,D,dt,A,BC-SC,FMWT);
%updateDLA = @(C,S,D) DLA3_HB_SSP(x,v,k,C,S,D,dt,Awave,BC,FMWT);
%updateDLA = @(C,S,D) DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA_Mot_FE(C,S,D,dt,Awave,BC,2*dt);
%updateDLA = @(C,S,D) DLA4_TAN_PROJ_RK2(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_TAN_FE(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_TAN_AB2(x,v,k,C,S,D,FC,FS,FD,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_TAN_BE(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_HB_SSP_RK3(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA_UC_SYS_RK2(x,v,k,C,S,D,dt,Awave,BC);
updateDLA = @(C,S,D) DLA_UC_EXP(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA5_HB_SSP(x,v,k,C,S,D,dt,Awave,BC);


if adapt
    %[C,S,D] = AdaptiveDLAResdiual_FE(x,v,k,C,S,D,dt,adapt_tol,Awave,BC);
    %[C,S,D] = AdaptiveDLAResdiual2_FE(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    %[C,S,D] = AdaptiveDLAResdiual2_FE_Lub(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    [C,S,D] = AdaptiveDLAResdiual_RA_FE(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    %[C,S,D] = AdaptiveDLAResdiual3_FE(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    %[C,S,D] = updateDLA(C,S,D);
    %[C,S,D] = initHierAdapt(C,S,D,hier_tol);
    %R = [R size(S,1)];
else
    [C,S,D] = updateDLA(C,S,D);  
    %[C,S,D,FC,FS,FD] = DLA4_TAN_AB2(x,v,k,FC,FS,FD,C,S,D,dt,Awave,BC);
    %[C1,S1,D1] = DLA4_TAN_FE(x,v,k,C,S,D,dt,Awave,BC);
    %[C2,S2,D2] = DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
end
R = [R size(S,1)];
r = size(S,1);

%% Full-grid run
if fullgrid
%UU_FE = UU - dt*applyMatA(UU,Awave);
%UU_RK2 = MAT_SSP_RK2(UU,dt,Awave,BC);
%UU_RK3 = MAT_SSP_RK3(UU,dt,Awave,BC);
%FU = -applyMatA(UU,Awave);
UU = MAT_SSP_RK3(UU,dt,Awave,BC);
%UU = computeMatrixExpon(UU,dt,Awave);
end

end

end


function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end

function U = MAT_FE(U,dt,Acell,BC)
    %Update by SSP-RK3
    U = U - dt*applyMatA(U,Acell);
    if BC.use
        for l=1:size(BC.cell1,1)
        U = U - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
    end
end

function U = MAT_SSP_RK2(U,dt,Acell,BC)
    %Update by SSP-RK3
    U1 = U - dt*applyMatA(U,Acell);
    if BC.use
        for l=1:size(BC.cell1,1)
        U1 = U1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
    end
    U = (1/2)*U + (1/2)*(U1-dt*applyMatA(U1,Acell));
    if BC.use
        for l=1:size(BC.cell1,1)
        U = U - (1/2)*dt*BC.cell3{l,1}*BC.cell3{l,2}*BC.cell3{l,3}';
        end
    end
end

function U = MAT_SSP_RK3(U,dt,Acell,BC)
    %Update by SSP-RK3
    U1 = U - dt*applyMatA(U,Acell);
    if BC.use
        for l=1:size(BC.cell1,1)
        U1 = U1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
    end
    U2 = (3/4)*U + (1/4)*(U1-dt*applyMatA(U1,Acell));
    if BC.use
        for l=1:size(BC.cell1,1)
        U2 = U2 - (1/4)*dt*BC.cell3{l,1}*BC.cell3{l,2}*BC.cell3{l,3}';
        end
    end
    U  = (1/3)*U + (2/3)*(U2-dt*applyMatA(U2,Acell));
    if BC.use
        for l=1:size(BC.cell1,1)
        U  = U  - (2/3)*dt*BC.cell2{l,1}*BC.cell2{l,2}*BC.cell2{l,3}';
        end
    end
end

function LU = applyMatA(U,Awave)
    LU = 0*U;
    for l=1:size(Awave,1)
        LU = LU + Awave{l,1}*U*Awave{l,2}';
    end
end

function sparseSVD(fcell,X,trans)
    NN = size(fcell,1);
    Y = zeros(size(X));
    if strcmp(trans,'notransp')
        for l=1:NN
            Y = Y + (fcell{l,1}*fcell{l,2})*(fcell{l,3}'*X);
        end
    else
        for l=1:NN
            Y = Y + (fcell{l,3}*fcell{l,2}')*(fcell{l,1}'*X);
        end        
    end

end

function M = assembleLRform(fcell)
    NN = size(fcell,1);
    M = zeros(size(fcell{1,1},1),size(fcell{1,3},1));
    for l=1:NN
        M = M + fcell{l,1}*fcell{l,2}*fcell{l,3}';
    end
end

function PU = calcTanProj(U,C,D)
    PU = (U*D)*D' - C*(C'*U*D)*D' + C*(C'*U);
end

function U = computeMatrixExpon(U,t,Awave)
    u = U;
    for i=1:20
        u_temp = zeros(size(u));
        for l=1:size(Awave,1)
            u_temp = u_temp - Awave{l,1}*u*Awave{l,2}';
        end
        u = u_temp;
        U = U + 1/factorial(i)*t^i*u;
        if 1/factorial(i)*t^i*norm(u,'fro') < 1e-15
            break
        end
    end
end

