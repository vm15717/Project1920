function [ansdat, timedat]=rkck_mk1(solver1,t,T,y,a,v,h,eptol,b)
tic
%a=1;
%v=1;
%T = [-10,10];
%h = (T(2)-T(1))/10000;
%t=T(1);
%y=[0;1];
ansdat=[];
timedat=[];
%eptol = 1e-6;
%b=0.8;
while t<T(end)
    k1 = h * solver1(t, y, a, v) ;
    k2 = h * solver1(t + (1/5) * h, y + (1/5) * k1,a,v) ;
    k3 = h * solver1(t + (3/10) * h, y + (3/40) * k1 + (9/40) * k2,a,v) ;
    k4 = h * solver1(t + (3/5) * h, y + (3/10) * k1 + (-9/10) * k2 +...
        (6/5) * k3,a,v) ;
    k5 = h * solver1(t + 1 * h, y + (-11/54) * k1 + (5/2) * k2 +...
        (-70/27) * k3+(35/27) * k4,a,v) ;
    k6 = h * solver1(t + (7/8) * h, y + (1631/55296) * k1 + (175/512) * k2 +...
        (575/13824) * k3+(44275/110592) * k4+(253/4096) * k5,a,v) ;
    %k4 = h * solver1(t + h, y + k3,a,v) ;
    % Update next value of y 
    y1 = y + (37/378)*(k1 )+ (250/621) * k3 + (125/594) * k4 + (512/1771)*k5 ;
    y2 = y + (2825/27648)*(k1 )+ (18575/48384) * k3 + (13525/55296) * k4 + (277/14336)*k5 +...
    (1/4)*k6 ;
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