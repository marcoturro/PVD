
function [vi,Pi,pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax)

a = z0/z1 ;
vmin = z0/t(end) ;
vi(1) = vmin;
i=1 ; 

while vi(i)<vmax
    vi(i+1) = vi(i)*a ;
    i = i+1 ;
end


[~, z1id] = min( abs( z + z1) );
[~, z0id] = min( abs( z + z0) );

dz = z(z1id)-z(z0id);

C = zeros(1,length(t)) ; 

for i = 1:length(t)
    C(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i)) ;
end

ti = z0./vi ;
Ci = interp1(t,C,ti,'pchip') ;
order = 3;
framelen = 11;
Ci = sgolayfilt(Ci,order,framelen);
ddtCi = my_2FD_non_uniform(ti,Ci) ;

C0 = Ci(1)/dz;

Pvmin = Ci(1)/dz/C0/vmin;


N = length(vi) ;

A = full(gallery('tridiag',N,-a^(-2),1,0)) ;

rhs = z0/C0*(ddtCi./vi.^3)' ;


rhs(1) = Pvmin ;
pi = linsolve(A,rhs);
pi(pi<0) = 0;
dvi = [vmin diff(vi)]';

vi(1) = 0;
Pi = pi.*dvi;
Pi = Pi/sum(Pi);


end
