m = 20+20i;
r = 0;
i = 0;
% for xr = -10:0.1:10
%     r = r+1;
%     for xi = -10:0.1:10
%         i = i + 1;
%         f(r,i) = -2*(xr^2-xi^2)+2*(real(m)*xr-imag(m)*xi);
%     end
% end

N = 40;
x=-N:0.1:N;  
y=-N:0.1:N;  

[xr,xi]=meshgrid(x,y);
f = -2*(xr.^2+xi.^2)+2*(real(m).*xr-imag(m).*xi);
mesh(xi,xr,f)