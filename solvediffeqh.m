%ft = linspace(-10, 100, 250); % Generate t for f 
%f = ft; % Generate f(t)
tic
a= 0.1:0.2:3;
g=1;
ansdat = [];
for i=1:length(a)
    TSPAN = [-50,50]; % Solve from t=1 to t=5
    IC = [0;1]; % y(t=0) = 1
    opts = odeset('RelTol',1e-16,'AbsTol',1e-16);
    [T Y] = ode45(@(t,y) myODE2nh(t, y, a(i), g), TSPAN, IC); % Solve ODE
    ansdat{i}=Y;
    timedat{i}=T;
end
for i=1:length(ansdat)
    figure;
    hold on;
    Y1 = ansdat{i}(:,1);
    Y2 = ansdat{i}(:,2);
    Ti=timedat{i};
    plot(Ti, abs(Y1).^(2));
    plot(Ti, abs(Y2).^(2));
    title('Plot of y as a function of time');
    legend('Y1','Y2');
    xlabel('Time'); ylabel('Y(t)');
end
toc