function dydt = myODE2h( t, y, a, v)
%f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
%g = interp1(gt, g, t); % Interpolate the data set (gt, g) at times t
dydt = zeros(2,1);
%f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
dydt(1) = (-a*y(1)*t+y(2)*v)*(-1i); % Evalute ODE at times t
dydt(2) = (a*y(2)*t+y(1)*v)*(-1i); % Evalute ODE at times t
end