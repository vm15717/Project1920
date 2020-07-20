clear
syms v12 v23 v13
syms t real
syms a1 a2 a3 real
syms p1(t) 
e11 = ((-1j*v12^2*v23)-(1j*v12*v13*a3*t)+(1j*v12*v13*a2*t)+(1j*v13^2*conj(v23)));
p3num = 1i*v12*diff(p1(t), t, t)+ diff(p1(t),1)*(-v1e2*a1*t-v12*a2*t-v13*conj(v23))+...
    p1*(v12*a1-1i*v12*(conj(v13)*v13)-1i*v12*(conj(v12)*v12)-a1*t*(1i*v12*a2*t+1i*v13*conj(v23)));
p3 = p3num/e11;
p2 = (1i*diff(p1(t),t)-a1*t*p1-v13*p3)/v12;
p2 = simplify(expand(p2));
eq1 = 1i*diff(p3,t)-conj(v13)*p1-conj(v23)*p2-a3*t*p3;
eq2 = simplify(expand(eq1));
t1 = collect(eq2, diff(p1(t),t,t,t));
t2 = collect(eq2, diff(p1(t),t,t));
t3 = collect(eq2, diff(p1(t),t));
t4 = collect(eq2, p1(t));
eqrc = t1*(diff(p1(t),t,t,t))+t2*(diff(p1(t),t,t))+t3*(diff(p1(t),t))+t4*p1(t);
eqcons = simplify(expand(eq2-eqrc));