function [vi,Pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax)

a = z0/z1;
vmin = z0/t(end);
vi(1) = vmin;
i=1;
while vi(i)<vmax
    vi(i+1) = vi(i)*a;
    i = i+1;
end

[~, z1id] = min( abs( z + z1) );
[~, z0id] = min( abs( z + z0) );

C = zeros(1,length(t)) ; 

% this section aims to find the spatial limits within which the relevant
% information is found.

for i = 1:length(t)
    C(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));
end


ti = z0./vi;

Ci = interp1(t,C,ti,'pchip')
ddtCi = my_2FD_non_uniform(ti,Ci);

C0 = Ci(1)/(z0-z1);

Pvmin = Ci(1)/(z0-z1)/C0/vmin;
N = length(vi);

A = full(gallery('tridiag',N,-a^(-2),1,0));

rhs = z0/C0*(ddtCi./vi.^3)';
rhs(1) = Pvmin ;
P = linsolve(A,rhs);
Pi = P;
dvi = [vmin diff(vi)]';
% figure
% plot(diff(vi));
vi(1) = 0;
Pi = P.*dvi;
Pi = Pi/sum(Pi);

end
