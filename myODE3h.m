function dydt = myODE3h( t, y, a, v)
%f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
%g = interp1(gt, g, t); % Interpolate the data set (gt, g) at times t
dydt = zeros(3,1);
%f = interp1(ft, f, t); % Interpolate the data set (ft, f) at times t
dydt(1) = (-a*y(1)*t+y(2)*v)*(-1i); % Evalute ODE at times t
dydt(2) = (y(3)*v+y(1)*v)*(-1i); % Evalute ODE at times t
dydt(3) = (a*y(3)*t+y(2)*v)*(-1i); % Evalute ODE at times t
end