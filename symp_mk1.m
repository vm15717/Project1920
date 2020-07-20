syms p1 p2 p3
syms v12 v23 v13
syms t real
syms a1 a2 a3 real
e1 = (((-v12*a1*t)-(v12*a2*t)-(v13*conj(v23)))*...
    ((-1j*v12^2*v23)-(1j*v12*v13*a3*t)+(1j*v12*v13*a2*t)+(1j*v13^2*conj(v23))));
e2 = (v12^2*v13*a3)-(v12^2*v13*a2);
ea1 = simplify(expand(e1-e2));
e11 = ((-1j*v12^2*v23)-(1j*v12*v13*a3*t)+(1j*v12*v13*a2*t)+(1j*v13^2*conj(v23)));
e3 = -v12*a1-v12*a2+v12*a1-1i*v12*(v13*conj(v13))-1i*v12*(v12*conj(v12))-...
    a1*t*(1i*v12*a2*t+1i*v13*conj(v23));
eb1 = e3*e11;
e4 = -v12*a1*t-v12*a2*t-v13*conj(v23);
e41 = (-1i*v12*v13*a3+1i*v12*v13*a2);
eb2 = e4*e41;
eans2 =simplify(expand(eb1-eb2));
e5 = -2*a1*1i*v12*a2*t-a1*1i*v13*conj(v23);
et1 = e5*e11;
e51 = v12*a1-1i*v12*(conj(v13)*v13)-1i*v12*(conj(v12)*v12)-...
    a1*t*(1i*v12*a2*t+1i*v13*conj(v23));
et2 = e51*e41;
eans3 = simplify(expand(et1-et2));