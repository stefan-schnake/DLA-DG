function [C,S,D] = DLA3_HB_SSP(x,v,k,C,S,D,dt,Awave,BC,FMWT)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%tic
%Convert BCs
BCc{1} = convertVectoMat(x,v,k,BC(:,1));
BCc{2} = convertVectoMat(x,v,k,BC(:,2));
BCc{3} = convertVectoMat(x,v,k,BC(:,3));

%%%Update K = CS
K = C*S;

%%%SSP-RK3
F = @(g) g - dt*applyAonK(g,D,Awave);
F1 = F(K) - dt*FMWT*(BCc{1}*(FMWT'*D));
F2 = (3/4)*K + (1/4)*(F(F1)- dt*FMWT*(BCc{3}*(FMWT'*D)));
K = (1/3)*K + (2/3)*(F(F2)- dt*FMWT*(BCc{2}*(FMWT'*D)));

%Run QR to get new K=C*S;
[C1,R1] = qr(K,0);
M = C1'*C;


%%%Update L = D*S'
L = D*S';

%%%SSP-RK3
F = @(g) g - dt*applyAonL(g,C,Awave);
F1 = F(L) - dt*FMWT*(BCc{1}'*(FMWT'*C));
F2 = (3/4)*L + (1/4)*(F(F1) - dt*FMWT*(BCc{3}'*(FMWT'*C)));
L = (1/3)*L + (2/3)*(F(F2) - dt*FMWT*(BCc{2}'*(FMWT'*C)));

%Run QR to get new L=D*S';
[D1,R2] = qr(L,0);
N = D1'*D;



%%%Update S
S_ = M*S*N';

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C1,D1,Awave);
F1 = F(S_) - dt*(FMWT'*C1)'*BCc{1}*(FMWT'*D1);
F2 = (3/4)*S_ + (1/4)*(F(F1) - dt*(FMWT'*C1)'*BCc{3}*(FMWT'*D1));
S1 = (1/3)*S_ + (2/3)*(F(F2) - dt*(FMWT'*C1)'*BCc{2}*(FMWT'*D1));

C = C1;
S = S1;
D = D1;
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After S update: ||u|| = %e\n',norm(u));


%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After V update: ||u|| = %e\n',norm(u));
%toc;
end


function LuD = applyAonK(K,D,Awave)
%%%Updates X via the DLA update

LuD = zeros(size(K));
for l=1:size(Awave,1)
    LuD = LuD + Awave{l,1}*K*(D'*Awave{l,2}'*D);
end

end

function CLuD = applyAonS(S,C,D,Awave)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + (C'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D);
end

end

function CLu = applyAonL(L,C,Awave)
%%%Updates V via the DLA update

CLu = zeros(size(L));
for l=1:size(Awave,1)
    CLu = CLu + Awave{l,2}*L*(C'*Awave{l,1}'*C);
end

end

