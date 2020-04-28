function [C,z,t] = data_norm(data,zL,vmin)
% Normalises the data to a standard shape in the z direct
% as well as adimensionalising t based on vmin and z_L
C = data.C;
z = data.z;
t = data.t;

sizeZ = 100;

if length(C(:,1)) ~= length(z)
    C = C';
end

C_zi = zeros(sizeZ,length(t));

for i=1:length(t)
    C_zi(:,i) = interp1(z,C(:,i),linspace(0,z(end),sizeZ))';
end

C_zi(C_zi<0) = 0;
z = linspace(0,z(end),sizeZ)/max(abs(z));
C = C_zi/max(max(C_zi));
t = t./(zL/vmin);


end