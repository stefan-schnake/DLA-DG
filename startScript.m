%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DLA RUN SCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


N = 256;

xx = [-1,1];vv = [-1,1];
x = MESH(xx(1),xx(2),N);
v = MESH(vv(1),vv(2),N);
k = 0;

per = reshape(convertVectoMat(x.nodes,v.nodes,k,1:x.N*x.N*(k+1)^2),[],1);
iper = zeros(size(per));
iper(per) = 1:numel(per);

r = 20;

test = 1;

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%CFL = 1/N^2;
%dt = 0.8*CFL;
%CFL = (5/4)*(1/73536);
%dt = 0.9*CFL;
%dt = pi/(8*ceil(1/CFL));
if exist('frac','var') == 0
    frac = 1/2;
end
dt = frac*CFL;
%dt = 0.05;
%dt = 0.05;

if exist('frac','var') == 0
    alg = 'RARA_UC';
end

fullgrid  = true;
moviebool = false;
plotbool  = true;
savebool  = false;

adapt  = false;
adapt_tol = 880*dt^2;
    %T = 1;  50dt^2 for RA; 12*dt^2 for WG
    %T = pi; 50dt^2 for RA;  5*dt^2 for WG
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%  Tolerance Table  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  err  | RARA_UC | RARA_TAN  | RARA_PJ | RA_UC  | WG
% 3e-3  | 880dt^2 | 11255dt^2 | 857dt^2 | 25dt^2 | 19dt^2
%       | 2.9546  | 2.5973    | 3.0877  | 3.0385 | 3.237
%
%
%  err  | RARA_UC | RARA_TAN  | RARA_PJ | RA_UC  | WG
% 2.138 | 400dt^2 | 500dt^2   | 200dt^2 |   dt^2 | 5dt^2
%   e-3 | 2.1348  | 2.1414    | 2.1642  | 2.5915 | 2.1082
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if exist('T','var') == 0
%T = 1;
T = pi;
%T = pi;
%T = 3.5;
end

time = 0;

%%%-------------------------------------------

%Clear figures
if plotbool
figure(5)
clf('reset')
%figure(6)
%clf('reset')
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
%FMWT = OperatorTwoScale_wavelet2(k+1,log2(N));

Acell = buildAdvectionMatrixWithBlocks2(x.nodes,v.nodes,k);
%Acell = cell(size(Acell));
%for i=1:numel(Acell)
%    Acell{i} = FMWT*Acell{i}*FMWT';
%end

%Use diffusion or advection
%A = Adv;
%A2 = Adv_back;
%A = LDG+N*Jmp;
%A2 = LDG+N*Jmp;
%A = LDG;
%A2 = LDG;

if test == 1
%init = @(x,y)  (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
box = @(x,y) (x > -1/2).*(x < 1/2).*(y > -1/2).*(y < 1/2);
init = @(x,y) soln_func(x,y,0,box);
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

u0 = buildNonSeparableSource(x.nodes,v.nodes,k,init);
u = u0;  

sigk = []; %smallest singular value

U = SOLN(x,v,k,per,iper);
U.vec = u0;
U.vec2Mat();
U.getLRFactors(r);
U.getMatFromLR();

R = [r];


UU = SOLN(x,v,k,per,iper); UU.mat = U.mat;
%UU_FE = SOLN(x,v,k,per,iper); UU.mat = U.mat;

if fullgrid
    myhist = [0;r;U.S(end,end);0;norm(U.mat-UU.mat,'fro');0;0];
else
    myhist = [0;r;U.S(end,end);norm(U.S,'fro');1e-5;0];
end

i = 0;
t = 0;

%% Time iteration
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

%% Create BC
if test == 1
    BC.use = 0;
    BC2.use = 0;
else
    BC.use = 1;
    BC.cell1 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t-dt));buildSeparableSourceMat(x,v,k,t-dt,source)];
    BC.cell2 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t-0.5*dt));buildSeparableSourceMat(x,v,k,t-0.5*dt,source)];
    BC.cell3 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t));buildSeparableSourceMat(x,v,k,t,source)];
    
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


%% Specify method
updateDLA = @(U) DLA4_HB_FE(U,dt,Acell,BC);

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
    U = updateDLA(U);  
end
time = time + toc;
R = [R size(U.S,1)];
r = size(U.S,1);


%% Full-grid run
if fullgrid



%UU_FE = UU_FE - dt*applyMatA(UU_FE,Acell);
%UU_RK2 = MAT_SSP_RK2(UU,dt,Acell,BC);
%UU_RK3 = MAT_SSP_RK3(UU,dt,Acell,BC);
%FU = -applyMatA(UU,Acell);
%UU = UU - dt*applyMatA(UU,Acell);
UU = MAT_SSP_RK3(UU,dt,Acell,BC);
%UU = MAT_SSP_RK2(UU,dt,Acell,BC);
%UU = computeMatrixExpon(UU,dt,Acell);

end

sing = svd(U.S);

if fullgrid
    UU.getLRFactors();
    %myhist = [myhist [i;r;sing(end);LTE;norm(U-UU,'fro');t]];
    myhist = [myhist [i;r;sing(end);0;norm(U.mat-UU.mat,'fro');t;0]];
else
    myhist = [myhist [i;r;sing(end);norm(U,'fro');1e-5;t]];
end

%% Error Calc and plotting
if plotbool && ( mod(i,ceil(T/(20*dt))) == 0 || i == 1 || lastplot)
%if plotbool && lastplot
    
    figure(5)
    U.getMatFromLR();
    U.mat2Vec();
    UU.mat2Vec();
    plotVecDualDiscrete(x.nodes,v.nodes,k,U.vec,UU.vec,myhist);
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
    
    
    figure(8)
    cla reset
    %plot(log10(abs(C)))
    %absC = sum(abs(C));
    %absD = sum(abs(D));
    %Xval = 1:numel(absC);
    %plot(Xval,absC,Xval,absD);   
    %semilogy(1:size(myhist,2),myhist(5,:),1:size(myhist,2),myhist(4,:))
    if fullgrid
        %semilogy(myhist(6,:),myhist(4,:),myhist(6,:),myhist(5,:))
        semilogy(myhist(6,:),myhist(4,:),myhist(6,:),myhist(5,:),myhist(6,:),myhist(7,:))
    end
    yyaxis right
    plot(myhist(6,:),myhist(2,:))
    yyaxis left
    %legend({'LTE','u-uu','rank'},'location','best')
    legend({'LTE','u-uu','uu-uu_lr','rank'},'location','best')
    

    %semilogy(sing);
 
    drawnow
    aa = 2;
end
end


if moviebool
    close(myVideo);
end


if adapt && plotbool  && 0
    figure(9)
    plot((1:numel(R))-1,R)
    title('Adaptive Rank plot');
    xlabel('iteration')
    ylabel('rank')
    
    figure(5)
    sgtitle("DLA Solution - Test = " + num2str(test) + "; tol = " + num2str(adapt_tol));
    %figure(6)
    %sgtitle("Full-Grid Solution - Test = " + num2str(test) + "; tol = " + num2str(adapt_tol));
    figure(8)
    title("Metrics - Test = " + num2str(test) + "; tol = " + num2str(adapt_tol) + ", N = "+num2str(N)+", k = "+num2str(k));
    
  
    
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

if ~exist('sigK','var')
    sigK = sigk;
elseif size(sigk,2) == size(sigK,2)
    sigK = [sigK;sigk];
end



function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end

function U = MAT_FE(U,dt,Acell,BC)
    %Update by SSP-RK3
    U.mat = U.mat - dt*applyMatA(U.mat,Acell);
    if BC.use
        for l=1:size(BC.cell1,1)
        U.mat = U.mat - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
    end
end

function U = MAT_SSP_RK2(U,dt,Acell,BC)
    %Update by SSP-RK3
    U1 = U.mat - dt*applyMatA(U.mat,Acell);
    if BC.use
        for l=1:size(BC.cell1,1)
        U1 = U1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
    end
    U.mat = (1/2)*U.mat + (1/2)*(U1-dt*applyMatA(U1,Acell));
    if BC.use
        for l=1:size(BC.cell1,1)
        U.mat = U.mat - (1/2)*dt*BC.cell3{l,1}*BC.cell3{l,2}*BC.cell3{l,3}';
        end
    end
end

function U = MAT_SSP_RK3(U,dt,Acell,BC)
    %Update by SSP-RK3
    U1 = U.mat - dt*applyMatA(U.mat,Acell);
    if BC.use
        for l=1:size(BC.cell1,1)
        U1 = U1 - dt*BC.cell1{l,1}*BC.cell1{l,2}*BC.cell1{l,3}';
        end
    end
    U2 = (3/4)*U.mat + (1/4)*(U1-dt*applyMatA(U1,Acell));
    if BC.use
        for l=1:size(BC.cell1,1)
        U2 = U2 - (1/4)*dt*BC.cell3{l,1}*BC.cell3{l,2}*BC.cell3{l,3}';
        end
    end
    U.mat = (1/3)*U.mat + (2/3)*(U2-dt*applyMatA(U2,Acell));
    if BC.use
        for l=1:size(BC.cell1,1)
        U.mat  = U.mat  - (2/3)*dt*BC.cell2{l,1}*BC.cell2{l,2}*BC.cell2{l,3}';
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