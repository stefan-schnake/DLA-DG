function [C,S,D] = DLA_UC_Adapt_RK2(x,v,k,C,S,D,dt,Awave,BC)
%% DLA using UC integrator by RK2 as a system in [C,S,D]
%%% using Adaptive DLA in the FE step

%tic
assert( BC.use == 0 );

%% First step

tol = dt^2;

%[C1,S1,D1] = AdaptiveDLAResdiual_RA_FE(x,v,k,C ,S ,D ,dt,tol,Awave,BC);
[C1,S1,D1] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C ,S ,D ,dt,tol,Awave,BC);

%% Second step

%[C1,S1,D1] = AdaptiveDLAResdiual_RA_FE(x,v,k,C1,S1,D1,dt,tol,Awave,BC);
[C1,S1,D1] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C1,S1,D1,dt,tol,Awave,BC);


[C2,RC] = qr([C C1],0);
[D2,RD] = qr([D D1],0);
S_fill = zeros(size(S,1),size(S1,1));
S2.mat = RC*(0.5*[S,S_fill;S_fill',S1])*RD';
[S2.U,S2.S,S2.V] = svd(S2.mat);

r = sum(diag(S2.S) > 0.5*dt^2);

C = C2*S2.U(:,1:r);
D = D2*S2.V(:,1:r);
S = S2.S(1:r,1:r);

end