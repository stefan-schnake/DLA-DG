% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% PDE:   
%         u_t + c\cdot\grad u + alpha(u-Pu) = 0
% 
% where
% 
%   c = (-y,x)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N = 256;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 0;

r = 1;

test = 1;

if exist('metaflag','var') == 0
    metaflag = false;
end

if ~metaflag
    frac = 1/2;
end
CFL = 1/(2*k+1)*(1/N);
dt = frac*CFL;

%Dissipation parameter
alpha = 5;

fy = @(y) 0*y+1;

if ~metaflag
    alg = 'RARA_UC';
end


fullgrid  = false;
moviebool = false;
plotbool  = false;
savebool  = false;

if metaflag
    adapt_tol = adapt_val*dt^2;
else
    adapt_tol = 50*dt^2;
end
adapt  = true;
%    T = 1;  50dt^2 for RA; 12*dt^2 for WG
%    T = pi; 50dt^2 for RA;  5*dt^2 for WG
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%  Tolerance Table  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%  err  | RARA_UC | RARA_TAN  | RARA_PJ | RA_UC  | WG
% 3e-3  | 880dt^2 | 11255dt^2 | 857dt^2 | 25dt^2 | 19dt^2
%       | 2.9546  | 2.5973    | 3.0877  | 3.0385 | 3.237
% 
% 
%  err  | RARA_UC | RARA_TAN  | RARA_PJ | RA_UC  | WG
% 2.138 | 400dt^2 | 500dt^2   | 200dt^2 |   dt^2 | 5dt^2
%   e-3 | 2.1348  | 2.1414    | 2.1642  | 2.5915 | 2.1082

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%T = 1;
%T = 2;
T = 6;
%T = 3.5;


time.DLA = 0;
time.FUL = 0;

%%-------------------------------------------

%Clear figures
if plotbool
figure(5)
clf('reset')
end

if moviebool
    myVideo = VideoWriter('movie.mp4','MPEG-4');
    myVideo.FrameRate = 1;
    open(myVideo);
end

Acell = buildAdvectionPeriodicMatrixWithBlocks(x,v,k);
%Acell(5,:) = {sqrt(alpha)*speye(size(Acell{1,1})),...
%             sqrt(alpha)*speye(size(Acell{1,2}))};
%Acell(6,:) = buildSeperableProjectionBlocks(x,v,k,fx,fy);
%Acell{6,1} = -sqrt(alpha)*Acell{6,1}; Acell{6,2} =  sqrt(alpha)*Acell{6,2};
Acell(5,:) = buildProjectionBlocks1D(x,v,k,fy,alpha);

if test == 1
box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
init = @(x,y) box(x,y);
BCsoln = @(x,v,t) soln_func(x,v,t,init);
source = @(x,v,t) 0*x;
end

u0 = buildNonSeparableSource(x,v,k,init);
u = u0;  

sigk = []; %smallest singular value

one = buildSeparableSource(x,v,k,@(x) 0*x+1,@(v) 0*v+1);
cons = one'*u0;

matu = convertVectoMat(x,v,k,u);
[U,Sig,V] = svd(matu);
C0 = U;
S0 = Sig;
D0 = V;

Sig = diag(Sig);
R = [r];
C = C0(:,1:r);
S = S0(1:r,1:r);
D = D0(:,1:r);

U = C*S*D';
UU = U;
UU_FE = U;

if fullgrid
    myhist = [0;r;S(end,end);0;norm(U-UU,'fro');0];
else
    myhist = [0;r;S(end,end);norm(U,'fro');1e-5;0];
end

%Create source
F.mat = convertVectoMat(x,v,k,buildNonSeparableSource(x,v,k,@(x,y) exp(-(x.^2+y.^2)./(abs(x)+1))));
[F.C,F.S,F.D] = svd(F.mat);
F.r = sum(diag(F.S) > 1e-14);
F.C = F.C(:,1:F.r); F.S = F.S(1:F.r,1:F.r); F.D = F.D(:,1:F.r);

i = 0;
t = 0;
% Time iteration
%while i < 1
while (t+dt <= T+1e-9) 
i = i+1;
t0 = t;
t = t + dt;
lastplot = 0;
if t+dt > T+1e-9; lastplot=1; end
if mod(i,ceil(T/(20*dt))) == 0 || i == 1
    fprintf("-------------------------------------------------------\n");
    fprintf("i=%d; t = %f\n",i,dt*i);
    fprintf("r = %d\n",r);
end

BC.use = 1;
BC.cell1 = {F.C,F.S,F.D};
BC.cell2 = {F.C,F.S,F.D};
BC.cell3 = {F.C,F.S,F.D};


% Specify method
updateDLA = @(C,S,D) DLA4_HB_FE(x,v,k,C,S,D,dt,Acell,BC);
%updateDLA = @(C,S,D) DLA4_IMEX_FE(x,v,k,C,S,D,dt,Acell,BC);
%updateDLA = @(C,S,D) DLA4_TAN_FE(x,v,k,C,S,D,dt,Acell,BC);
%updateDLA = @(C,S,D) DLA4_HB_SSP_RK3(x,v,k,C,S,D,dt,Acell,BC);
%updateDLA = @(C,S,D) DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Acell,BC);

tic
if adapt
    switch alg
        case 'LUB'
            [C,S,D] = AdaptiveDLAResdiual2_FE_Lub(x,v,k,C,S,D,dt,adapt_tol,Acell,BC);
        case 'WG'
            [C,S,D] = WG_FE(x,v,k,C,S,D,dt,Acell,BC,adapt_tol);
        case 'RARA_UC'
            [C,S,D] = AdaptiveDLAResdiual_RA_FE(x,v,k,C,S,D,dt,adapt_tol,Acell,BC);
        case 'RARA_TAN'
            [C,S,D] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C,S,D,dt,adapt_tol,Acell,BC);
        case 'RARA_PROJ'
            [C,S,D] = AdaptiveDLAResdiual_Proj_RA_FE(x,v,k,C,S,D,dt,adapt_tol,Acell,BC);
    end
else
    [C,S,D] = updateDLA(C,S,D);  
end
time.DLA = time.DLA + toc;
R = [R size(S,1)];
r = size(S,1);


% Full-grid run
if fullgrid


UU_FE = MAT_FE(UU_FE,dt,Acell,BC);
%UU_FE = UU_FE - dt*applyMatA(UU_FE,Acell);
%UU_RK2 = MAT_SSP_RK2(UU,dt,Acell,BC);
%UU_RK3 = MAT_SSP_RK3(UU,dt,Acell,BC);
%FU = -applyMatA(UU,Acell);
%UU = UU - dt*applyMatA(UU,Acell);
tic
UU = MAT_SSP_RK3(UU,dt,Acell,BC);
time.FUL = time.FUL + toc;
%UU = MAT_SSP_RK2(UU,dt,Acell,BC);
%UU = computeMatrixExpon(UU,dt,Acell);

%LTEU = MAT_SSP_RK2(U,dt,Acell,BC);
U = C*S*D';
%LTE = norm(LTEU-U,'fro');
LTE = 0;

end

sing = svd(S);

if fullgrid
    %[UUU,SSS,VVV] = svd(UU);
    %UU_lr = UUU(:,1:r)*SSS(1:r,1:r)*VVV(:,1:r)';
    myhist = [myhist [i;r;sing(end);LTE;norm(U-UU,'fro');t]];
    %myhist = [myhist [i;r;sing(end);LTE;norm(U-UU,'fro');t;norm(UU-UU_lr,'fro')]];
else
    myhist = [myhist [i;r;sing(end);norm(U,'fro');1e-5;t]];
end

% Error Calc and plotting
if plotbool && ( mod(i,ceil(T/(20*dt))) == 0 || i == 1 || lastplot)
%if plotbool && lastplot
    


    u = convertMattoVec(x,v,k,C*S*D');
    uu = convertMattoVec(x,v,k,UU);
    figure(5)
    plotVecDualDiscrete(x,v,k,u,uu,myhist);
    sgtitle("DLRA vs Full Rank - t = " + num2str(t));
    
    if moviebool
        frame = getframe(gcf);
        writeVideo(myVideo,frame);
    end
    %figure(6)
    %plotVec(x,v,k,uu,@(x,y) BCsoln(x,y,t));
    %sgtitle('Full-Grid Solution')
    %figure(7)
    %plotVec(x,v,k,uu_lr,@(x,y) soln(x,y,t));
    %sgtitle('Full-Grid-Low-Rank Solution')
    %plotVec(x,v,k,convertMattoVec(x,v,k,(FMWT'*C(:,1:2))*S(1:2,1:2)*(FMWT'*D(:,1:2))'),@(x,y) soln(x,y,i*dt));
    %sgtitle('Truncated DLA Solution, r=2')
    
    %figure(6)
    %semilogy(sigk);
    
    
    %figure(8)
    %cla reset
    %plot(log10(abs(C)))
    %absC = sum(abs(C));
    %absD = sum(abs(D));
    %Xval = 1:numel(absC);
    %plot(Xval,absC,Xval,absD);   
    %semilogy(1:size(myhist,2),myhist(5,:),1:size(myhist,2),myhist(4,:))
%     if fullgrid
%         %semilogy(myhist(6,:),myhist(4,:),myhist(6,:),myhist(5,:))
%         semilogy(myhist(6,:),myhist(4,:),myhist(6,:),myhist(5,:),myhist(6,:),myhist(7,:))
%     end
%     yyaxis right
%     plot(myhist(6,:),myhist(2,:))
%     yyaxis left
%     %legend({'LTE','u-uu','rank'},'location','best')
%     legend({'LTE','u-uu','uu-uu_lr','rank'},'location','best')
%     
% 
%     semilogy(sing);
%  
     drawnow
     aa = 2;
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

function LU = applyMatA(U,Acell)
    LU = 0*U;
    for l=1:size(Acell,1)
        LU = LU + Acell{l,1}*U*Acell{l,2}';
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

function U = computeMatrixExpon(U,t,Acell)
    u = U;
    for i=1:20
        u_temp = zeros(size(u));
        for l=1:size(Acell,1)
            u_temp = u_temp - Acell{l,1}*u*Acell{l,2}';
        end
        u = u_temp;
        U = U + 1/factorial(i)*t^i*u;
        if 1/factorial(i)*t^i*norm(u,'fro') < 1e-15
            break
        end
    end
end


