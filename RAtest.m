n = 256;
k = 40;
p = 40;
pp = k+p;

R.S = diag([ones(k+p,1);zeros(n-(k+p),1)]);
load RA_dat.mat
R.U = UU;
R.V = VV;
R.mat = R.U*R.S*R.V;

err = zeros(1000,1);
for iter=1:1000
Q = rand(n,pp);
Q = R.mat*Q;
[Q,~] = qr(Q,0);
B = (R.mat'*Q)';
[Utmp,S1,D1] = svd(B,'econ');
C1 = Q*Utmp;
R_RA = C1*S1*D1';
err(iter) = norm(R.mat-R_RA,2);
end

bound = (1+11*sqrt(n*pp))*R.S(k+1,k+1);
fprintf('k = %2d, sig_{k+1}(R) = %f\n',k,R.S(k+1,k+1));
fprintf('error is %5.4e, bound is %e\n',norm(err,inf),bound);