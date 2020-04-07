function Cmat = produce_data(P,t,z)


Nv = length(P);
vs = linspace(0,1,Nv);
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));

for k = 1:length(t)
    
    C = 0;
    for i = 1:Nv
        C = ci(P(i),vs(i),z,t(k)) + C;
    end
    Cmat(k,:) =C;
end
end