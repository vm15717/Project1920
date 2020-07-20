f1 = @(x,t,T,y,a,v,h,eptol,b) rkhf_mk1(x,t,T,y,a,v,h,eptol,b) ;
f12 = @(x,t,T,y,a,v,h,eptol,b) runge_kutta_mk1(x,t,T,y,a,v,h,eptol,b) ;
f11 = @(x,t,T,y,a,v,h,eptol,b) rkck_mk1(x,t,T,y,a,v,h,eptol,b) ;
f2 = @(t,y,a,v) myODE2h(t,y,a,v) ;
f3 = @(t,y,a,v) myODE2nh(t,y,a,v) ;
f4 = @(t,y,a,v) myODE3h(t,y,a,v) ;
f5 = @(t,y,a,v) myODE3nh(t,y,a,v) ;
a=1;
v=1;
T = [-10,10];
h = (T(2)-T(1))/10000;
t=T(1);
y=[1;0];
b=0.9;
eptol =1e-6;
[ansdat, timedat]= f12(f2,t,T,y,a,v,h,eptol,b);
figure;
hold on
plot(timedat, abs(ansdat).^2);
prob = abs(ansdat).^2;
plot(timedat, sum(prob,2));
legend('Y1','Y2')
xlabel('Time')
ylabel('Y(t)')