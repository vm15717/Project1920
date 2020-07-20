%function eigLnmk1(N)
N=5;    
syms a v t;
[S_x,S_z] =hamiltonmat(N);
H0 = a*t*S_z;
V0 = 1i*v*S_x;
H1 = H0+V0;
[egvec, egval] = eig(H1);
aval = 1;
% vval = 1;
% tval = 0.5;
% eg1 = double(diag(subs(egval,[a,v,t],[aval,vval,tval])));
% H2 = H1([1,3,2,4,5],:);
% [egvec, egval] = eig(H2);
% eg2 = double(diag(subs(egval,[a,v,t],[aval,vval,tval])));
% H3 = H1([4,1,3,5,2],:);
% [egvec, egval] = eig(H3);
% eg3 = double(diag(subs(egval,[a,v,t],[aval,vval,tval])));
tval = -10:0.1:10;
egs = diag(subs(egval,[a,v],[aval,vval]));
egs_all = double(subs(egs,t,tval));
evs = (subs(egvec,[a,v],[aval,vval]));
evs_all = double(subs(evs,t,tval));
evs_all2 = evs_all./vecnorm(evs_all);
figure
hold on
plot(tval,real(egs_all))
plot(tval,imag(egs_all))
xlabel('Time parameter')
ylabel('Eigenvalues')
%end