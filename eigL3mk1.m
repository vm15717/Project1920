% 3 Level Landau Zener
% Eigenvalues and Eigenvectors of the Hamiltonian
% Hamiltonian H = [a*t v 0;v 0 v;0 v a*t]
syms a v t;
H1 = [a*t v 0;v 0 v;0 v -a*t];
[egvec, egval] = eig(H1);
aval = 1;
vval = 1;
tval = -10:0.1:10;
egs = diag(subs(egval,[a,v],[aval,vval]));
eg1 = egs(1);
eg2 = egs(2);
eg3 = egs(3);
eg1val = double(subs(eg1,t,tval));
eg2val = double(subs(eg2,t,tval));
eg3val = double(subs(eg3,t,tval));
evs1 = diag(subs(egvec(:,1),[a,v],[aval,vval]));
evs2 = diag(subs(egvec(:,2),[a,v],[aval,vval]));
evs3 = diag(subs(egvec(:,3),[a,v],[aval,vval]));
ev1val = double(subs(evs1,t,tval));
ev2val = double(subs(evs2,t,tval));
ev3val = double(subs(evs3,t,tval));
figure;
hold on
plot(tval, eg1val);
plot(tval, eg2val);
plot(tval, eg3val);
xlabel('Time parameter')
ylabel('Eigenvalues')