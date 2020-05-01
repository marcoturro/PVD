function ddy = my_2FD_non_uniform(x,y)

N = length(x);
% h = x(2:end)-x(1:end-1);
% h1 = h(1:end-1);
% h2 = h(2:end);

for i = 2:N-1
    h1 = x(i)-x(i-1);

    h2 = x(i+1)-x(i);
    ddy(i) = (2* h2.*y(i-1) - 2 * (h1+h2).*y(i) + 2*h1.*y(i+1))./ (h1.*h2.*(h1+h2));
%ddy(2:N-1) = (2* h2.*y(1:end-2) - 2 * (h1+h2).*y(2:end-1) + 2*h1.*y(3:end))./ (h1.*h2.*(h1+h2));
end
ddy(1) = (2*y(1)-5*y(2)+4*y(3)-y(4))/(x(2)-x(1))^2;
ddy(N) = (2*y(N)-5*y(N-1)+4*y(N-2)-y(N-3))/(x(N)-x(N-1))^2;

end