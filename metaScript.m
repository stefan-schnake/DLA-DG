
frac = 1/2;
hierScript
err = norm(U-UU,'fro');
myrank{1} = [myhist(2,:);myhist(6,:)];

frac = 1/4;
hierScript
err = [err norm(U-UU,'fro')];
myrank{2} = [myhist(2,:);myhist(6,:)];

frac = 1/8;
hierScript
err = [err norm(U-UU,'fro')];
myrank{3} = [myhist(2,:);myhist(6,:)];

frac = 1/16;
hierScript
err = [err norm(U-UU,'fro')];
myrank{4} = [myhist(2,:);myhist(6,:)];

log2(err(1:end-1)./err(2:end))

plot(myrank{1}(2,:),myrank{1}(1,:),myrank{2}(2,:),myrank{2}(1,:),myrank{3}(2,:),myrank{3}(1,:),myrank{4}(2,:),myrank{4}(1,:));
legend({'1/2','1/4','1/8','1/16'})
title('SSP-RK2 Adaptive TAN Rank')
xlabel('time')
ylabel('rank')