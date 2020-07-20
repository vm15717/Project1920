% 2 Level Landau Zener
% Eigenvalues and Eigenvectors of the Hamiltonian
% Hamiltonian H = [a*t v;v -a*t]
syms a v t;
H1 = [a*t v*1i;v*1i -a*t];
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
innerprod = (sum(conj(ev1val).*ev2val,1));
figure;
hold on
plot(tval, real(eg1val),'r','DisplayName',strcat('\Re(E_{1})'));
plot(tval, real(eg2val),'k','DisplayName',strcat('\Re(E_{2})'));
plot(tval, imag(eg1val),'r','DisplayName',strcat('\Im(E_{1})'),'LineStyle','--');
plot(tval, imag(eg2val),'k','DisplayName',strcat('\Im(E_{2})'),'LineStyle','--');
xlabel('Time parameter')
ylabel('Eigenvalues')
legend()
figure;
hold on
%plot(tval, real(innerprod));
%plot(tval, imag(innerprod));
plot(tval, abs(innerprod));
xlabel('Time parameter')
ylabel('Innerproduct of eigenvectors')
%Start from initial state 
% evinit = [0 ;1];
% innerprod2 = (sum(conj(evinit).*ev2val,1));
% figure;
% hold on
% %plot(tval, real(innerprod2));
% %plot(tval, imag(innerprod2));
% plot(tval, abs(innerprod2));
% xlabel('Time parameter')
% ylabel('Innerproduct of eigenvectors')