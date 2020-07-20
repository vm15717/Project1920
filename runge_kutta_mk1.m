function [ansdat, timedat]=runge_kutta_mk1(solver1,t,T,y,a,v,h,eptol,b)
tic
%a=1;
%v=1;
%T = [-50,50];
%h = (T(2)-T(1))/10000;
%t=T(1);
%y=[1;0];
ansdat=[];
timedat=[];
while t<=T(end)
    k1 = h * solver1(t, y, a, v) ;
    k2 = h * solver1(t + 0.5 * h, y + 0.5 * k1,a,v) ;
    k3 = h * solver1(t + 0.5 * h, y + 0.5 * k2,a,v) ;
    k4 = h * solver1(t + h, y + k3,a,v) ;
    % Update next value of y 
    y = y + (1.0 / 6.0)*(k1 + 2 * k2 + 2 * k3 + k4) ;
    % Update next value of x 
    ansdat=vertcat(ansdat,y');
    timedat = vertcat(timedat,t);
    t = t + h ;
end
%plot(timedat, abs(ansdat));
toc