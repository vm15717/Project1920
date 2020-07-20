% 2 Level Landau Zener
% Eigenvalues and Eigenvectors of the Hamiltonian
% Hamiltonian H = [a*t v;v -a*t]
syms a v t;
H1 = [a*t v;v -a*t];
[egvec, egval] = eig(H1);
aval = 1;
vval = 1;
tval = -10:0.1:10;
egs = diag(subs(egval,[a,v],[aval,vval]));
eg1 = egs(1);
eg2 = egs(2);
eg1val = double(subs(eg1,t,tval));
eg2val = double(subs(eg2,t,tval));
evs1 = (subs(egvec(:,1),[a,v],[aval,vval]));
evs2 = (subs(egvec(:,2),[a,v],[aval,vval]));
evs1 = evs1/norm(evs1);
evs2 = evs2/norm(evs2);
ev1val = double(subs(evs1,t,tval));
ev2val = double(subs(evs2,t,tval));
innerprod = sum(ev1val.*ev2val,1);
figure;
hold on
plot(tval, eg1val);
plot(tval, eg2val);
xlabel('Time parameter')
ylabel('Eigenvalues')
figure;
plot(tval, innerprod);
xlabel('Time parameter')
ylabel('Innerproduct of eigenvectors')