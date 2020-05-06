function [v, P] = PVD_solve(solve)

t = solve.t;
z = solve.z;
Ci = solve.Ci;
c0 = solve.c0;
z0 = solve.z0;
z1 = solve.z1;
tvmin = solve.tvmin;
tvmax = solve.tvmax;
dv = solve.dv;
Pvmin = solve.Pvmin;
Pvmax = solve.Pvmax;

vmax= z1/tvmax ;
vmin = z0/tvmin;

v = vmax:dv:vmin;
ti = z0./v;

Ci_i = interp1(t,Ci,ti);

dddtCint = my_2FD_non_uniform(ti,Ci_i);

[v, P] = P_RHS_solve(z0,z1,ti,dddtCint,Pvmax,Pvmin);
end

