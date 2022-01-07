%%% Script for Error and Rate Calculations using DLA

N = 256;
k = 0;
xx = [-1,1];
vv = [-1,1];
test = 1;
fullgrid = 1;

r = 20;

dtVec = 2.^(-1:-1:-5);
load UU_EXP
UU = UU_EXP;

errVec.DLA = 0*dtVec;
errVec.FUL = 0*dtVec;

for i=1:numel(dtVec)
    fprintf("i = %2d\n",i);
    if fullgrid
        [C,S,D,U_FUL] = hierFunc(N,k,dtVec(i),xx,vv,r,test);
        errVec.FUL(i) = norm(U_FUL-UU,'fro');
    else
        [C,S,D] = hierFunc(N,k,dtVec(i),xx,vv,r,test);
    end
    errVec.DLA(i) = norm(C*S*D'-UU,'fro');
end

rateVec.DLA = log2(errVec.DLA(1:end-1)./errVec.DLA(2:end));
rateVec.FUL = log2(errVec.FUL(1:end-1)./errVec.FUL(2:end));