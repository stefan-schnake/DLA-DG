

%% Order computations
mivec = 1:4;
for mii = 1:numel(mivec)
    mi = mivec(mii);
    frac = 1/2^mi;
    constAdvLR
    tocvec(mii) = time;
    U = C*S*D';
    load("fullrank/constLR-k0-N256-T2.mat",'UU');
    err_FE(mii) = norm(UU_FE-UU,'fro');
    err(mii) = norm(U-UU,'fro');
    myrank{mii} = [myhist(6,:);myhist(2,:)];
    %FMWT = OperatorTwoScale_wavelet2(1,8);
    %plotVecDualDiscrete(x,v,k,convertMattoVec(x,v,k,FMWT'*U*FMWT),convertMattoVec(x,v,k,FMWT'*UU*FMWT));
    %save("fullrank/frac"+num2str(2^mi)+"T"+num2str(T)+".mat",'UU');
    
end

MM = numel(mivec);
fprintf('Error\t');
for mii=1:MM
    fprintf('%e',err(mii));
    if mii ~= MM
        fprintf(',');
    else
        fprintf('\n');
    end
end
fprintf('Rank\t');
for mii=1:MM
    fprintf('%i',max(myrank{mii}(2,:)));
    if mii ~= MM
        fprintf(',');
    else
        fprintf('\n');
    end
end
