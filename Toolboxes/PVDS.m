function [v, P] = PVDS(solve)

t = solve.t;
z = solve.z;
Ci = solve.Ci;
c0 = solve.c0;
z0 = solve.z0;
z1 = solve.z1;
tvmin = solve.tvmin;
tvmax = solve.tvmax;
Pvmin = solve.Pvmin;
Pvmax = solve.Pvmax;

vmax= z1/tvmax ;
vmin = z0/tvmin;

a = z1/z0;
i = 2;
v(1) = vmax;

while 1
    vi = (a^i)*v(i-1);
    if vi < vmin
        break
    end
    v(i)=vi;
    i = i+1;
end

% v = linspace(vmin,vmax,100);


ti = z0./v;

Ci_i = interp1(t,Ci,ti);

try
    dddtCint = my_2FD_non_uniform(ti,Ci_i);
catch
    dddtCint = 0;
end

A = zeros(length(ti));
A = A + diag(ones(length(ti),1)) - diag(ones(length(ti)-1,1),-1)*a^2;
rhs = (dddtCint./v.^3)';
rhs(1) = Pvmax ; %vmax condition
rhs(end) = Pvmin ; %vmin condition
P = linsolve(A,rhs);

% P = rhs./(1-a);
P(P<0) = 0;
P = P/sum(P);
end

