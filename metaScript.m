

%% Order computations
% for mi=4
%     frac = 1/2^mi;
%     hierScript
%     tocvec(mi) = time;
%     U = C*S*D';
%     load("fullrank/frac"+num2str(2^mi)+"T"+num2str(T)+".mat");
%     err(mi) = norm(U-UU,'fro');
%     myrank{mi} = [myhist(6,:);myhist(2,:)];
%     %FMWT = OperatorTwoScale_wavelet2(1,8);
%     %plotVecDualDiscrete(x,v,k,convertMattoVec(x,v,k,FMWT'*U*FMWT),convertMattoVec(x,v,k,FMWT'*UU*FMWT));
%     %save("fullrank/frac"+num2str(2^mi)+"T"+num2str(T)+".mat",'UU');
% end

%% UC error plot routine
% for r=1:35
%     frac = 1/16;
%     hierScript
%     U = C*S*D';
%     load("fullrank/frac"+num2str(16)+"T"+num2str(T)+".mat");
%     err(r) = norm(U-UU,'fro');
% end

%% Randomized Test

%mytol = [880,11255,857];
mytol = [400,500,200];
alg_list = {'RARA_UC','RARA_TAN','RARA_PROJ'};

for ii = 1:3
for mi=1:50
    alg = alg_list{ii};
    disp(mi)
    frac = 1/16;
    hierScript
    U = C*S*D';
    %load("fullrank/frac16T1.mat");
    load("fullrank/frac16T3.1416.mat");
    myrank(mi) = r;
    err(mi) = norm(U-UU,'fro'); 
end
    randval{ii} = err;
end



%tol = 50dt^2 produces first order method

% log2(err(1:end-1)./err(2:end))

% plot(myrank{1}(2,:),myrank{1}(1,:),'-',myrank{2}(2,:),myrank{2}(1,:),':',myrank{3}(2,:),myrank{3}(1,:),'--',myrank{4}(2,:),myrank{4}(1,:),'.-','LineWidth',1.5);
% legend({'dt = \omega/2','dt = \omega/4','dt = \omega/8','dt = \omega/16'})
% %title('SSP-RK2 Adaptive TAN Rank')
% title('Plot 3 -- RA Residual Adaptivity - Unconventional Integrator -- tol = 50dt^2, \omega = 1/256')
% xlim([0,pi])
% xticks(0:pi/4:pi)
% xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
% xlabel('time')
% ylabel('rank')