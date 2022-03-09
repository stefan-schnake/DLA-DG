T = pi;

for mi=1:4
    frac = 1/2^mi;
    hierScript
    tocvec(mi) = time;
    U = C*S*D';
    load("fullrank/frac"+num2str(2^mi)+"T"+num2str(T)+".mat");
    err(mi) = norm(U-UU,'fro');
    myrank{mi} = [myhist(2,:);myhist(6,:)];
    %save("fullrank/frac"+num2str(2^mi)+"T"+num2str(T)+".mat",'UU');
end


%tol = 50dt^2 produces first order method

log2(err(1:end-1)./err(2:end))

plot(myrank{1}(2,:),myrank{1}(1,:),'-',myrank{2}(2,:),myrank{2}(1,:),':',myrank{3}(2,:),myrank{3}(1,:),'--',myrank{4}(2,:),myrank{4}(1,:),'.-','LineWidth',1.5);
legend({'dt = CFL/2','dt = CFL/4','dt = CFL/8','dt = CFL/16'})
%title('SSP-RK2 Adaptive TAN Rank')
xlim([0,pi])
xticks(0:pi/4:pi)
xticklabels({'0','\pi/4','\pi/2','3\pi/4','\pi'})
xlabel('time')
ylabel('rank')