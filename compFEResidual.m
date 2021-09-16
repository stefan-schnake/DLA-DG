function [R] = compFEResidual(C0,S0,D0,C,S,D,dt,Acell,BC)
%Compute norm of low-memory obtained resdiual
%  R = U_FE - U_UC where 
%  U_UC = CSD' and U_FE = U - dt*F(U)
%  U = C0*S0*D0'
% Order of basis is [C0 kry(C0) C]

[CC,nodesC,kryC] = krylovData(Acell(:,1),C0);
[DD,nodesD,kryD] = krylovData(Acell(:,2),D0);
r = size(S,1);

mm = size(kryC{1},1);
nn = size(kryD{1},1);

R = zeros(mm+r,nn+r);
%Add U
R(1:r,1:r) = R(1:r,1:r) + S0;
%Compute F(U)
for l=1:size(Acell,1)
    R(1:mm,1:nn) = R(1:mm,1:nn) -dt*kryC{l}*S0*kryD{l}';
end


RC = CC'*C;
PC = C - CC*RC;
[PC,RRC] = qr(PC,0);

RD = DD'*D;
PD = D - DD*RD;
[PD,RRD] = qr(PD,0);

%Subtract off DLA soln
R(1:mm,1:nn) = R(1:mm,1:nn) - RC*S*RD';
R(mm+1:end,1:nn) = R(mm+1:end,1:nn) - RRC*S*RD';
R(1:mm,nn+1:end) = R(1:mm,nn+1:end) - RC*S*RRD';
R(mm+1:end,nn+1:end) = R(mm+1:end,nn+1:end) - RRC*S*RRD';


%R = [CC PC]*R*[DD PD]';



end

function [VV,nodes,kry] = krylovData(Acell,C)

NN = size(Acell,1);
r = size(C,2);
kry = cell(NN,1);
VV = C;
nodes = zeros(NN,2);
    
for i=1:NN
    A = Acell{i};
    W = A*C;
    H = C'*W;
    kry{i,1} = H;
    W = W - C*H;
    for j=1:i-1 %Project off previous vector
%         V = A*C;
%         V = V - C*kry{j,1}(1:r,:);
%         for k=1:j-1
%             curC = Acell{k}*C - C*kry{k,1}(1:r,:);
%             for l=1:k-1
%                 curC = curC - 
%             end
%         end
%         
%         [V,~] = qr(V,0);
        V = VV(:,nodes(j,1):nodes(j,2));
        H = V'*W;
        kry{i,1} = [kry{i,1};H];
        W = W - V*H;
    end
    [svdU,svdS,svdV] = svd(W,0);
    rr = sum(diag(svdS) > 1e-10);
    Q = svdU(:,1:rr);
    R = svdS(1:rr,1:rr)*svdV(:,1:rr)';
    kry{i,1} = [kry{i,1};R];
    VV = [VV Q];
    if i == 1
        nodes(i,1) = r+1;
        nodes(i,2) = nodes(i,1)+rr-1;
    else
        nodes(i,1) = nodes(i-1,2)+1;
        nodes(i,2) = nodes(i,1)+rr-1;
    end
    for j=1:i-1
        kry{j,1} = [kry{j,1};0*R];
    end
end

end

function CC = recoverKrylovBasis(C0,kry,nodes,Acell)
    n = nodes(end,2);
    r = size(C0,2);
    NN = size(Acell,1);
    
    CC = zeros(size(C0,1),n);
    CC(:,1:size(C0,2)) = C0;
    
    for i=1:NN
        st = nodes(i,1);
        ed = nodes(i,2);
        tmp = Acell{i}*C0 - C0*kry{i}(1:r,:);
        for j=1:i-1
            iter = nodes(j,1):nodes(j,2);
            tmp = tmp - CC(:,iter)*kry{i}(iter,:);
        end    
        CC(:,st:ed) = tmp/kry{i}(st:ed,:);
    end
end

