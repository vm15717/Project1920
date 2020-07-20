function dydt = myODE2nh(t, y,a,g)
%f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
%g = interp1(gt, g, t); % Interpolate the data set (gt, g) at times t
dydt = zeros(2,1);
%f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
dydt(1) = (-a*y(1)*t+y(2)*g*(1i))*(-1i); % Evalute ODE at times t
dydt(2) = (a*y(2)*t+y(1)*g*(1i))*(-1i); % Evalute ODE at times t
end