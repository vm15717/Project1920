% 3 Level Landau Zener
% Eigenvalues and Eigenvectors of the Hamiltonian
% Hamiltonian H = [a*t v 0;v 0 v;0 v a*t]
syms a v t;
H1 = [a*t 1i*v 0;1i*v 0 1i*v;0 1i*v -a*t];
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
evs1 = (subs(egvec(:,1),[a,v],[aval,vval]));
evs2 = (subs(egvec(:,2),[a,v],[aval,vval]));
evs3 = (subs(egvec(:,3),[a,v],[aval,vval]));
evs1 = evs1/norm(evs1);
evs2 = evs2/norm(evs2);
evs3 = evs3/norm(evs3);
ev1val = double(subs(evs1,t,tval));
ev2val = double(subs(evs2,t,tval));
ev3val = double(subs(evs3,t,tval));
innerprod12 = sum((conj(ev1val).*ev2val),1);
innerprod23 = sum((conj(ev2val).*ev3val),1);
innerprod31 = sum((conj(ev3val).*ev1val),1);
figure;
hold on
plot(tval, real(eg1val),'r','DisplayName',strcat('\Re(E_{1})'));
plot(tval, real(eg2val),'g','DisplayName',strcat('\Re(E_{2})'));
plot(tval, real(eg3val),'b','DisplayName',strcat('\Re(E_{3})'));
plot(tval, imag(eg1val),'r','DisplayName',strcat('\Im(E_{1})'),'LineStyle','--');
plot(tval, imag(eg2val),'g','DisplayName',strcat('\Im(E_{2})'),'LineStyle','--');
plot(tval, imag(eg3val),'b','DisplayName',strcat('\Im(E_{3})'),'LineStyle','--');
legend()
xlabel('Time parameter')
ylabel('Eigenvalues')
figure;
hold on;
plot(tval, abs(innerprod12),'ro');
plot(tval, abs(innerprod23),'g');
plot(tval, abs(innerprod31),'b');
legend('Inner Product1','Inner Product 2', 'Inner Product 3');
xlabel('Time parameter')
ylabel('Innerproduct of eigenvectors')