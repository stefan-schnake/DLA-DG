

%% Order computations
% mivec = 1:4;
% for mii = 1:numel(mivec)
%     mi = mivec(mii);
%     frac = 1/2^mi;
%     dissScript
%     tocvec(mii) = time;
%     U = C*S*D';
%     %load("fullrank/constLR-k0-N256-T2.mat",'UU');
%     load("fullrank/diss-a5-T"+num2str(T)+"-F-1.mat",'UU');
%     err_FE(mii) = norm(UU_FE-UU,'fro');
%     err(mii) = norm(U-UU,'fro');
%     myrank{mii} = [myhist(6,:);myhist(2,:)];
%     %FMWT = OperatorTwoScale_wavelet2(1,8);
%     %plotVecDualDiscrete(x,v,k,convertMattoVec(x,v,k,FMWT'*U*FMWT),convertMattoVec(x,v,k,FMWT'*UU*FMWT));
%     %save("fullrank/frac"+num2str(2^mi)+"T"+num2str(T)+".mat",'UU');
%     
% end
% 
% MM = numel(mivec);
% fprintf('Error\t');
% for mii=1:MM
%     fprintf('%e',err(mii));
%     if mii ~= MM
%         fprintf(',');
%     else
%         fprintf('\n');
%     end
% end
% fprintf('Rank\t');
% for mii=1:MM
%     fprintf('%i',max(myrank{mii}(2,:)));
%     if mii ~= MM
%         fprintf(',');
%     else
%         fprintf('\n');
%     end
% end

alglist = {"RARA_UC","RARA_TAN","RARA_PROJ","WG","LUB"};
tollist = [15,100,15,9,50];

for ii=1:numel(alglist)
    metaflag = true;
    frac = 1/16;
    alg = alglist{ii};
    adapt_val = tollist(ii);
    dissScript
    U = C*S*D';
    load("fullrank/diss-a5-T"+num2str(T)+"-F-1.mat",'UU');
    err(ii) = norm(U-UU,'fro');
    myrank{ii} = [myhist(6,:);myhist(2,:)];
end
pvals = 1:400:numel(myrank{1}(1,:));
p = plot(myrank{1}(1,pvals),myrank{1}(2,pvals),'-*',...
     myrank{2}(1,pvals),myrank{2}(2,pvals),'-b',...
     myrank{3}(1,pvals),myrank{3}(2,pvals),'-d',...
     myrank{4}(1,pvals),myrank{4}(2,pvals),'-+',...
     myrank{5}(1,pvals),myrank{5}(2,pvals),'-o','LineWidth',1);
p(1).Color = '#EDB120';
p(2).Color = '#7E2F8E';
p(3).Color = '#77AC30';
p(4).Color = '#0072BD';
p(5).Color = '#D95319';
legend({strjoin(["RARA\_UC ",num2str(err(1))]),...
        strjoin(["RARA\_Tan",num2str(err(2))]),...
        strjoin(["RARA\_Proj",num2str(err(3))]),...
        strjoin(["FRA",num2str(err(4))]),...
        strjoin(["RA\_UC",num2str(err(5))])});
xlabel('t');
ylabel('rank');