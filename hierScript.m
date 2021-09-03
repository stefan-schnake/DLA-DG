
N = 64;

xx = [-1,1];vv = [-1,1];
x = xx(1):(xx(2)-xx(1))/N:xx(2);
v = vv(1):(vv(2)-vv(1))/N:vv(2);
k = 2;

r = 40;

r_cut = 1;

plotbool = true;
savebool = false;

adapt  = false;
adapt2 = false;
adapt_tol = 1e-3;

adapt3 = false;

adapt4 = false;

hier_adapt = false;
hier_tol = 0.55;

test = 1;

%dt = 0.5*(2/N^2);
%dt = .1;
CFL = (max(abs(xx))+max(abs(vv)))/(2*k+1)*(1/N);
%CFL = 1/N^2;
%dt = 0.8*CFL;
%CFL = (5/4)*(1/73536);
%dt = 0.9*CFL;
dt = pi/(8*ceil(1/CFL));
%dt = 0.05;
%dt = 0.05;

%T = .5;
T = pi;
%T = 3.5;

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

%L = buildLDGMatrixAlt(x,v,k);
%LDG = L'*L;
%Jmp = buildJumpMatrix(x,v,k,1);
Adv = buildAdvectionMatrix(x,v,k,true);
%Adv_back = buildBackAdvectionMatrix(x,v,k,true);
FMWT = OperatorTwoScale_wavelet2(k+1,log2(N));

Acell = buildAdvectionMatrixWithBlocks(x,v,k);
Acell2 = buildAdvectionMatrixWithBlocks2(x,v,k);
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
cons = one'*u0;

%Run a few timesteps to smooth the solution out
%for i=1:5
%u = u + dt*LDG*u;
%end


matu = convertVectoMat(x,v,k,u);
[U,Sig,V] = svd(matu);
C0 = FMWT*U;
S0 = Sig;
D0 = FMWT*V;

Sig = diag(Sig);
if (adapt || adapt2 || adapt3)
    r = sum(Sig > adapt_tol)+1;
    if r > size(U,1); r = size(U,1); end
    fprintf('-- Initial Adaptive r: r = %d\n',r);
    R = [r];
elseif hier_adapt
    [C,S,D] = initHierAdapt(C0,S0,D0,hier_tol);
    r = size(S,1);
    fprintf('-- Initial Adaptive r: r = %d\n',r);
    R = [r];
elseif adapt3
    [~,Stemp,~,~] = lrSVDApprox(C0,S0,D0,adapt3_tol);
    r = size(Stemp,1);
    R = [r];
else
    R = [r]; 
end
%r = 10;
C = C0(:,1:r);
S = S0(1:r,1:r);
D = D0(:,1:r);

%myhist = [0;r;S0(r,r);norm(u-uu)];
myhist = [];

% if hier_adapt
%     [C,S,D] = initHierAdapt(FMWT*U,Sig,FMWT*V,hier_tol);
%     r = size(S,1);
%     fprintf('-- Initial Adaptive r: r = %d\n',r);
%     R = [r];
%     %Sig = svd(S);
%     %adapt_tol = S(end);
% else
%     R = [r]; 
%     C = FMWT*U(:,1:r);
%     D = FMWT*V(:,1:r);
%     S = Sig(1:r,1:r);
% end

%Get rank5 approximation


rank1 = log10(Sig(1))-log10(Sig(2));

%Convert from realspace to waveletspace
%C = FMWT*C;
%D = FMWT*D;

U = C*S*D';
UU = U;
%u = convertMattoVec(x,v,k,FMWT'*C*S*D'*FMWT);
uu = u;
i = 0;
t = 0;
%while i < 5
while (t+dt <= T+1e-9) 
i = i+1;
t0 = t;
t = t + dt;
lastplot = 0;
if t+dt > T+1e-9; lastplot=1; end
fprintf("-------------------------------------------------------\n");
fprintf("i=%d; t = %f\n",i,dt*i);
%BC = [buildDirichletBC(x,v,k,@(x,y) soln(x,y,t-dt)), ...
%      buildDirichletBC(x,v,k,@(x,y) soln(x,y,t-0.5*dt)), ...
%      buildDirichletBC(x,v,k,@(x,y) soln(x,y,t))];
%SC  = [buildNonSeparableSource(x,v,k,@(x,y) source(x,y,t-dt)),...
%      buildNonSeparableSource(x,v,k,@(x,y) source(x,y,t-0.5*dt)),...
%      buildNonSeparableSource(x,v,k,@(x,y) source(x,y,t))];
if test == 1
    BC.use = 0;    
else
    BC.use = 1;
%     RHS = @(t,dt) ...
%     [buildDirichletBC(x,v,k,@(x,y) BCsoln(x,y,t-dt)), ...
%      buildDirichletBC(x,v,k,@(x,y) BCsoln(x,y,t-0.5*dt)), ...
%      buildDirichletBC(x,v,k,@(x,y) BCsoln(x,y,t))] - ...
%     [buildNonSeparableSource(x,v,k,@(x,y) source_vec(x,y,t-dt)),...
%      buildNonSeparableSource(x,v,k,@(x,y) source_vec(x,y,t-0.5*dt)),...
%      buildNonSeparableSource(x,v,k,@(x,y) source_vec(x,y,t))];
%     BC.vec = RHS(t,dt);
    BC.cell1 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t-dt));buildSeparableSourceMat(x,v,k,t-dt,source)];
    BC.cell2 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t-0.5*dt));buildSeparableSourceMat(x,v,k,t-0.5*dt,source)];
    BC.cell3 = [buildDirichletMatBC(x,v,k,@(x,y) BCsoln(x,y,t));buildSeparableSourceMat(x,v,k,t,source)];
    %Convert BCs to Wavespace
    for l=1:size(BC.cell1,1)
       BC.cell1{l,1} = FMWT*BC.cell1{l,1}; 
       BC.cell1{l,3} = FMWT*BC.cell1{l,3}; 
       BC.cell3{l,1} = FMWT*BC.cell3{l,1}; 
       BC.cell3{l,3} = FMWT*BC.cell3{l,3}; 
    end
end
 
%Get indicator at this timestep
%BC = RHS(t,dt);
%sourceFMat = FMWT*convertVectoMat(x,v,k,BC(:,1))*FMWT';
%tic;
%indi = svd(computeResidualMatrix(C,S,D,Awave,sourceFMat));
%%indi = svds(@(x,trans) computeResidual(C,S,D,Awave,sourceFMat,x,trans),[1,1]*size(C,1),size(C,1));
%time1 = toc;
%fprintf('Indicator is %e and took %f seconds\n',indi(1),time1);

%BC = buildDirichletBC(x,v,k,@(x,y) soln(x,y,i*dt));
      
%BC_old = buildDirichletBC(x,v,k,@(x,y) soln(x,y,(i-1)*dt));
%BC_hlf = buildDirichletBC(x,v,k,@(x,y) soln(x,y,(i-1/2)*dt));

%updateDLA = @(C,S,D) DLA(x,v,k,C,S,D,dt,A,A2,BC);
%updateDLA = @(C,S,D) DLA_CN(x,v,k,C,S,D,dt,A,A2,BC);
%updateDLA = @(C,S,D) DLA2_HB_SSP(x,v,k,C,S,D,dt,A,BC-SC,FMWT);
%updateDLA = @(C,S,D) DLA3_HB_SSP(x,v,k,C,S,D,dt,Awave,BC,FMWT);
updateDLA = @(C,S,D) DLA4_HB_SSP_RK2(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_HB_FE(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA5_HB_SSP(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA4_HB_SSP_RK2_woproj(x,v,k,C,S,D,dt,Awave,BC);
%updateDLA = @(C,S,D) DLA3_HB_FE(x,v,k,C,S,D,dt,Awave,BC,FMWT);
%updateDLA = @(C,S,D) DLA3_HB_FE_Iteration(x,v,k,C,S,D,dt,Awave,BC,FMWT);
%updateDLA = @(C,S,D) DLA2_HB_BE(x,v,k,C,S,D,dt,A,BC-SC,FMWT);
%updateDLA = @(C,S,D) DLA2_HB_CN(x,v,k,C,S,D,dt,A,BC-SC,FMWT);
%updateDLA = @(C,S,D) DLA2_SSP(x,v,k,C,S,D,dt,A,BC);
%updateDLA = @(C,S,D) DLA2_BE(x,v,k,C,S,D,dt,A,BC);


if adapt
    [C,S,D] = AdaptiveDLAWithInit(C,S,D,updateDLA,adapt_tol,i,C0,S0,D0);
    %[C,S,D] = updateDLA(C,S,D);
    %[C,S,D] = initHierAdapt(C,S,D,hier_tol);
    %R = [R size(S,1)];
elseif adapt2
    %u = convertMattoVec(x,v,k,(FMWT'*C)*S*(FMWT'*D)');
    %plusOneDLA = @(C,S,D) DLA2_HB_BE(x,v,k,C,S,D,5,speye(size(A))+dt*A,-repmat(u,1,3) + dt*(BC-SC),FMWT);
    %[C,S,D] = AdaptiveDLAWithInitv2Hier(C,S,D,updateDLA,plusOneDLA,adapt_tol,i,C0,S0,D0,false);
    %[C,S,D] = updateDLA(C,S,D);
    %[C,S,D,t] = AdaptiveDLAWithInitv4(x,v,k,C,S,D,t-dt,t,updateDLA,adapt_tol,i,C0,S0,D0,Acell,RHS,FMWT);
    %[C,S,D] = AdaptiveDLAWithInitv7(x,v,k,C,S,D,t,updateDLA,adapt_tol,Awave,RHS(t,dt),FMWT);
    [C,S,D] = AdaptiveDLAWithInitv8(x,v,k,C,S,D,t,dt,updateDLA,adapt_tol,Awave,BC,FMWT);
elseif adapt3
    %[C,S,D,t] = AdaptiveDLAWithInitv6(x,v,k,C,S,D,t-dt,t,updateDLA,adapt3_tol,i,C0,S0,D0,Acell,RHS,FMWT);
    %[C,S,D] = AdaptiveDLAResdiual2_SSP_RK2(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    %[C,S,D] = AdaptiveDLAResdiual2_FE(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    [C,S,D] = AdaptiveDLAResdiual3_FE(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    %[C,S,D] = AdaptiveDLAResdiual2_FE_Lub(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
    %[C,S,D] = AdaptiveDLAResdiual_FE2(x,v,k,C,S,D,dt,adapt_tol,Awave,BC,FMWT);
elseif adapt4
    [C,S,D] = DLA4_HB_SSP(x,v,k,C,S,D,dt,Awave,RHS(t,dt),FMWT,1e-7);
elseif hier_adapt
    [C,S,D] = updateDLA(C,S,D,t);
    [C,S,D] = initHierAdapt(C,S,D,hier_tol);
else
    [C,S,D] = updateDLA(C,S,D);
    
end
R = [R size(S,1)];
r = size(S,1);
fprintf("r = %d\n",r);

sing = svd(S);
%rank1 = [rank1 log10(sing(1))-log10(sing(2))];
%fprintf('r = %d\n',r);
%fprintf('Min singular value of S: %e\n',min(sing));

%% Full-grid run
fullgrid = 1;
if fullgrid
%uu = pcg(speye(size(LDG))+dt*LDG,uu,1e-13);
%%%BE
%uu = (speye(size(A))+dt*A)\(uu - dt*(BC(:,3)-SC(:,3)));
%%%CN
%uu = (speye(size(A))+dt/2*A)\( (speye(size(A))-dt/2*A)*uu - 0.5*dt*(BC(:,1)-SC(:,1)+BC(:,3)-SC(:,3)));
%uu = (speye(size(LDG))+dt/2*A)\( (speye(size(LDG))-dt/2*A)*uu );
%uu = uu - dt*A*uu;
%uu = expm(-A*(i*dt))*u0;
%%%SSP-RK3
%u1 = uu - dt*A*uu - dt*BC(:,1);
%u2 = (3/4)*uu + (1/4)*(u1-dt*A*u1-dt*BC(:,3));
%uu = (1/3)*uu + (2/3)*(u2-dt*A*u2-dt*BC(:,2));
%uu = u1;

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

%

end



%LTEu1 = u - dt*A*u - dt*BC(:,1);
%LTEu2 = (3/4)*u + (1/4)*(LTEu1-dt*A*LTEu1-dt*BC(:,3));
%LTEu  = (1/3)*u + (2/3)*(LTEu2-dt*A*LTEu2-dt*BC(:,2));
%LTEu = LTEu1;
%if BC.use
%    LTEu1 = u - dt*A*u - dt*BC.vec(:,1);
%    LTEu = (1/2)*u + (1/2)*(LTEu1-dt*A*LTEu1-dt*BC.vec(:,3));
%else
%    LTEu1 = u - dt*A*u;
%    LTEu = (1/2)*u + (1/2)*(LTEu1-dt*A*LTEu1);
%end
%LTEu = LTEu1;

%%%FE
%uu = uu - dt*A*uu - dt*RHS_t(:,1);

%u = convertMattoVec(x,v,k,(FMWT'*C)*S*(FMWT'*D)');
%LTE = norm(LTEu-u);


%dis_l2 = norm(u-uu);

%Full Grid Low Rank
%[UU,SS,VV] = svd(convertVectoMat(x,v,k,uu));
%r_cut = r;
%UU_lr = UU(:,1:r_cut);
%SS_lr = SS(1:r_cut,1:r_cut);
%VV_lr = VV(:,1:r_cut);
%uu_lr = convertMattoVec(x,v,k,UU_lr*SS_lr*VV_lr');
%SS = diag(SS);

%myhist [i;rank;smallest singular value;LTE;residual;time];
if fullgrid
myhist = [myhist [i;r;sing(end);LTE;norm(U-UU,'fro');t]];
else
myhist = [myhist [i;r;sing(end);norm(U,'fro');1e-5;t]];
end

%% Error Calc and plotting
%fprintf('-- Norm of BE: %e\n',norm(uu));
%fprintf('- Errors:\n');
%fprintf('-- Error of DLA   and L2 Proj: %e\n',norm(u-u_sol))
%fprintf('-- Error of BE    and L2 Proj: %e\n',norm(uu-u_sol))
%fprintf('-- Error of BE_lr and L2 Proj: %e\n',norm(uu_lr-u_sol))
if plotbool && ( mod(i,ceil(T/(10*dt))) == 0 || i == 1 || lastplot)
%if plotbool
    u = convertMattoVec(x,v,k,FMWT'*C*S*D'*FMWT);
    uu = convertMattoVec(x,v,k,FMWT'*UU*FMWT);
    %u_sol = buildNonSeparableSource(x,v,k,@(x,y) soln(x,y,i*dt));
    figure(5)
    plotVec(x,v,k,u,@(x,y) BCsoln(x,y,t));
    sgtitle('DLA Solution')
    %figure(6)
    %plotVec(x,v,k,uu,@(x,y) BCsoln(x,y,t));
    %sgtitle('Full-Grid Solution')
    %figure(7)
    %plotVec(x,v,k,uu_lr,@(x,y) soln(x,y,t));
    %sgtitle('Full-Grid-Low-Rank Solution')
    %plotVec(x,v,k,convertMattoVec(x,v,k,(FMWT'*C(:,1:2))*S(1:2,1:2)*(FMWT'*D(:,1:2))'),@(x,y) soln(x,y,i*dt));
    %sgtitle('Truncated DLA Solution, r=2')
    
    
    figure(8)
    cla reset
    %plot(log10(abs(C)))
    %absC = sum(abs(C));
    %absD = sum(abs(D));
    %Xval = 1:numel(absC);
    %plot(Xval,absC,Xval,absD);   
    %semilogy(1:size(myhist,2),myhist(5,:),1:size(myhist,2),myhist(4,:))
    semilogy(myhist(6,:),myhist(4,:),myhist(6,:),myhist(5,:))
    yyaxis right
    plot(myhist(6,:),myhist(2,:))
    yyaxis left
    legend({'LTE','u-uu','rank'},'location','best')
    

    %semilogy(sing);
 
    drawnow
    aa = 2;
end

%err = u-u_sol;
%[Cerr,Serr,Derr] = svd(convertVectoMat(x,v,k,err));
%Serr = diag(Serr);
%Cerr = FMWT*Cerr;
%Derr = FMWT*Derr;
%figure(10)
%NN = size(Cerr,1);
%plot(1:NN,sum(abs(Cerr)),1:NN,sum(abs(Derr)));
%disp(Serr(1:10)');
end


if (hier_adapt+adapt+adapt2+adapt3+adapt4) && plotbool
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

function z = soln_func(x,v,t,u0)
    %Follow characteristics backwards to initial condition        
    x0 =  cos(t)*x+sin(t)*v;
    v0 = -sin(t)*x+cos(t)*v;
    
    z = u0(x0,v0);
end

function u = update_RK3_SSP(u,dt,A,SF)
    %Update by SSP-RK3
    u1 = u - dt*A*u - dt*SF(:,1);
    u2 = (3/4)*u + (1/4)*(u1-dt*A*u1-dt*SF(:,3));
    u  = (1/3)*u + (2/3)*(u2-dt*A*u2-dt*SF(:,2));
end

function LU = applyMatA(U,start,Awave)
    LU = 0*U;
    for l=start:size(Awave,1)
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
