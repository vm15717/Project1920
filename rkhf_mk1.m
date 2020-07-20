function [ansdat, timedat]=rkhf_mk1(solver1,t,T,y,a,v,h,eptol,b)
tic
% a=1;
% v=1;
% T = [-50,50];
% h = (T(2)-T(1))/10000;
% t=T(1);
% y=[1;0;0];
ansdat=[];
timedat=[];
%eptol = 1e-6;
%b=0.9;
while t<=T(end)
    k1 = h * solver1(t, y, a, v) ;
    k2 = h * solver1(t + 0.25 * h, y + 0.25 * k1,a,v) ;
    k3 = h * solver1(t + (3/8) * h, y + (3/32) * k1 + (9/32) * k2,a,v) ;
    k4 = h * solver1(t + (12/13) * h, y + (1932/2197) * k1 + (-7200/2197) * k2 +...
        (7296/2197) * k3,a,v) ;
    k5 = h * solver1(t + 1 * h, y + (439/216) * k1 + (-8) * k2 +...
        (3680/513) * k3-(845/4104) * k4,a,v) ;
    k6 = h * solver1(t + (1/2) * h, y + (-8/27) * k1 + (2) * k2 +...
        (-3544/2565) * k3+(1859/4104) * k4-(11/40) * k5,a,v) ;
    %k4 = h * solver1(t + h, y + k3,a,v) ;
    % Update next value of y 
    y1 = y + (25/216)*(k1 )+ (1408/2565) * k3 + (2197/4104) * k4 + (-1/5)*k5 ;
    y2 = y + (16/135)*(k1 )+ (6656/12825) * k3 + (28561/56430) * k4 + (-9/50)*k5 +(2/55)*k6 ;
    %error
    ep1 = mean(abs(y1-y2));
    if ep1>=eptol
        h = b*h*(eptol/ep1)^(0.2);
    else
         h = b*h*(eptol/ep1)^(0.25);
    end
    % Update next value of x 
    ansdat=vertcat(ansdat,y1');
    y=y1;
    timedat = vertcat(timedat,t);
    t = t + h ;
end
%plot(timedat, abs(ansdat).^2);
toc