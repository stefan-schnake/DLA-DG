function [C1,S1,D1] = DLA_UC_Adapt_RK3(x,v,k,C,S,D,dt,Awave,BC)
%% DLA using UC integrator by RK2 as a system in [C,S,D]
%%% using Adaptive DLA in the FE step

%tic
assert( BC.use == 0 );

%% First step

tol = 0.5*dt^3;

[C1,S1,D1] = AdaptiveDLAResdiual_RA_FE(x,v,k,C ,S ,D ,dt,tol,Awave,BC);
%[C1,S1,D1] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C ,S ,D ,dt,tol,Awave,BC);

%% Second step

[C1,S1,D1] = AdaptiveDLAResdiual_RA_FE(x,v,k,C1,S1,D1,dt,tol,Awave,BC);
%[C1,S1,D1] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C1,S1,D1,dt,tol,Awave,BC);


[C1,RC] = qr([C C1],0);
[D1,RD] = qr([D D1],0);
S_fill = zeros(size(S,1),size(S1,1));
Sf.mat = RC*([0.75*S,S_fill;S_fill',0.25*S1])*RD';
[Sf.U,Sf.S,Sf.V] = svd(Sf.mat);

r = sum(diag(Sf.S) > tol);

C1 = C1*Sf.U(:,1:r);
D1 = D1*Sf.V(:,1:r);
S1 = Sf.S(1:r,1:r);

%% Third step

[C1,S1,D1] = AdaptiveDLAResdiual_RA_FE(x,v,k,C1,S1,D1,dt,tol,Awave,BC);
%[C1,S1,D1] = AdaptiveDLAResdiual_TAN_RA_FE(x,v,k,C1,S1,D1,dt,tol,Awave,BC);


[C1,RC] = qr([C C1],0);
[D1,RD] = qr([D D1],0);
S_fill = zeros(size(S,1),size(S1,1));
Sf.mat = RC*([1/3*S,S_fill;S_fill',2/3*S1])*RD';
[Sf.U,Sf.S,Sf.V] = svd(Sf.mat);

r = sum(diag(Sf.S) > tol);

C1 = C1*Sf.U(:,1:r);
D1 = D1*Sf.V(:,1:r);
S1 = Sf.S(1:r,1:r);

end