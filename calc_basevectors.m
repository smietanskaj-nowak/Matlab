function [av, bv, cv]=calc_basevectors(a, b, c, alpha, beta, gamma)
av=a*[1 0 0];
bv=b*[cos(gamma/360*2*pi), sin(gamma/360*2*pi), 0];
cx=c*cos(beta/360*2*pi);
cy=c*(cos(alpha/360*2*pi)-cos(beta/360*2*pi)*cos(gamma/360*2*pi))/sin(gamma/360*2*pi);
cz=sqrt(c^2-cx^2-cy^2);
cv=[cx, cy, cz];
end