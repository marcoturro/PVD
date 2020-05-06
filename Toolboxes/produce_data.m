function Cmat = produce_data(P,vs,t,z)


Nv = length(P);
ci = @(ci0,vs,z,t)(ci0'.*(1-heaviside(z+vs'*t)));
Cmat = zeros(length(t),length(z));
for k = 1:length(t)
    Cmat(k,:) = sum(ci(P,vs,z,t(k)));
end

end