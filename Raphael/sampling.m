clear all;
close all;


vvmin = [1e-7,1e-6,1e-5];
vmax = 1e-3;

mrks = {'x','+','d'}

a = 1.05;

for j = 1:length(vvmin)
vmin = vvmin(j)
i = 1;
v = 0;
v(1) = vmin;
while v(i)<vmax
    v(i+1) = v(i)*a;
    i = i+1;
end

y = ones(1,length(v));
figure(1)
plot(v,v,mrks{j}); hold on;

end