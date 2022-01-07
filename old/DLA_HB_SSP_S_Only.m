function [S] = DLA_HB_SSP_S_Only(x,v,k,C1,S,D1,dt,Awave,BC,FMWT)
%% DLA in hierarchical basis
% Does not convert soln to full basis but evaluates the PDE on
% the x and v components and then sums up

%This iteration has an flag to not compute the basis update, but to use the
%basis given

%tic
%Convert BCs
BCc{1} = convertVectoMat(x,v,k,BC(:,1));
BCc{2} = convertVectoMat(x,v,k,BC(:,2));
BCc{3} = convertVectoMat(x,v,k,BC(:,3));

%%%Update S
S_ = S;

%%%SSP-RK3
F = @(g) g - dt*applyAonS(g,C1,D1,Awave);
F1 = F(S_) - dt*(FMWT'*C1)'*BCc{1}*(FMWT'*D1);
F2 = (3/4)*S_ + (1/4)*(F(F1) - dt*(FMWT'*C1)'*BCc{3}*(FMWT'*D1));
S1 = (1/3)*S_ + (2/3)*(F(F2) - dt*(FMWT'*C1)'*BCc{2}*(FMWT'*D1));

S = S1;
%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After S update: ||u|| = %e\n',norm(u));


%u = convertMattoVec(x,v,k,C*S*D');
%fprintf('-- After V update: ||u|| = %e\n',norm(u));
%toc;
end

function CLuD = applyAonS(S,C,D,Awave)
%%%Updates S via the DLA update

CLuD = zeros(size(S));
for l=1:size(Awave,1)
    CLuD = CLuD + (C'*Awave{l,1}*C)*S*(D'*Awave{l,2}'*D);
end

end
