function [v, P] = PVD_solve(dat)
clearvars -except dat Disp

t = dat.t;
z = dat.z;
Ci = dat.Ci;
c0 = dat.c0;
z0 = dat.z0;
z1 = dat.z1;
tvmin = dat.tvmin;
tvmax = dat.tvmax;
dv = dat.dv;
Pvmin = dat.Pvmin;
Pvmax = dat.Pvmax;
%C_v = dat.C_v;
Disp = dat.Disp;

 if ~exist('Disp','var')
     % third parameter does not exist, so default it to something
      Disp = 0;
 end


vmax= z1/tvmax ;
vmin = z0/tvmin;

v = vmax:dv:vmin;
ti = z0./v;

Ci_i = interp1(t,Ci,ti);

dddtCint = my_2FD_non_uniform(ti,Ci_i);

[v, P] = P_RHS_solve(z0,z1,ti,dddtCint,Pvmax,Pvmin);
end

