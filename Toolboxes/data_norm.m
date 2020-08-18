function [C,z,t,v,P] = data_norm(data,vmax)
% Normalises the data to a standard shape as well as
% adimensionalising the time vector based on the maximum 
% velocity (vmax) and the lowest space point (zL)

C = data.C;
z = data.z;
t = data.t;

sizeZ = 100;
sizeP = 350;

%%%%%%%%%%%%% Process of C matrix %%%%%%%%%%%%%%%

if length(C(:,1)) ~= length(z)
    C = C';
end

C = wdenoise(C,6);
C_zi = zeros(sizeZ,length(t));
zad = linspace(0,z(end),sizeZ);

for i=1:length(t)
    C_zi(:,i) = interp1(z,C(:,i),zad)';
end

if z(1) ~=0
    fnan = find(isnan(C_zi(:,1)),1,'last');
    for i=1:length(t)
        C_zi(1:fnan,i) = interp1([0 z(1)],[0 C_zi(fnan+1,i)],zad(1:fnan));
    end
end

C_zi(C_zi<0) = 0;
C = C_zi/max(max(C_zi));

%%%%%%%%%%%%% Process z vector %%%%%%%%%%%%%%%%%%

zL = max(abs(z));
z = zad/zL;

%%%%%%%%%%%%%% Process t vector %%%%%%%%%%%%%%%%%

t = t/(zL/vmax);

%%%%% Process, if existing v and P vectors %%%%%%

try
    v = data.v;
    P = data.P;
    v = v/vmax;
    P = interp1([v],[P],linspace(0,1,sizeP));
    v = linspace(0,1,sizeP);
    P = P/max(data.P);
catch
    v = zeros(1,sizeP);
    P = zeros(1,sizeP);
end

end