
f = @(x) x.^2.*cos(2*x) +1;
%n=14;
figure;
for n= 5  : 5 : 65
subplot(2,2,1);
a =-4;
b =-1;
fplot(f,[a,b]);

x_i = linspace(a,b,n+1);%для разбиения отрезка
y_i = f(x_i);
hold on; 
plot(x_i,y_i,'*r');

coef = polyfit(x_i,y_i,n);%для построения полинома
p = @(x) polyval(coef,x);
xx = linspace(a,b,1000);
p_xx = p(xx);
plot(xx,p_xx,'--g');
grid on;
legend('function','data','polynomial');
title(['Number of nodes = ', num2str(n+1)]);

hold off;
error = abs(p_xx - f(xx));
subplot(2,2,3);
semilogy(xx,error);
grid on;

error_max = max(error);
title(['Max error = ', num2str(error_max)]);

subplot(2,2,[2,4]);
semilogy(n+1, error_max,'or');
hold on;
grid on
xlabel('Number of nodes');
ylabel('Max error');
pause(1);
end