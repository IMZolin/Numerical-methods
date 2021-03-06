mkr1 = load('files\f_mkr1.txt');
mkr2 = load('files\f_mkr2.txt');
runge1 = load('files\f_runge1.txt');
runge2 = load('files\f_runge2.txt');
error_mkr1 = load('files\f_error_mkr1.txt');
error_mkr2 = load('files\f_error_mkr2.txt');
error_runge1 = load('files\f_error_runge1.txt');
error_runge2 = load('files\f_error_runge2.txt');
x1 = load('files\f_x1.txt');
x2 = load('files\f_x2.txt');
x = load('files\f_x.txt');
y = load('files\f_y.txt');
outrage_mkr = load('files\f_outrage_mkr.txt');
outrage_mkr_error = load('files\f_outrage_mkr_error.txt');
outrage_runge = load('files\f_outrage_mkr.txt');
outrage_runge_error = load('files\f_outrage_runge_error.txt');
steps = load('files\f_steps.txt');
y_h1 = load('files\f_y_h1.txt');
y_h2 = load('files\f_y_h2.txt');
steps_mkr = load('files\f_steps_mkr.txt');
steps_runge = load('files\f_steps_runge.txt');

figure;
subplot(1,2,1);
plot(x,y,'r');
hold on;
plot(x1,mkr1, 'm');
plot(x2,mkr2, 'g');
plot(x1, runge1,'c');
plot(x2, runge2,'y');
xlabel('x');
ylabel('y');
title('Graphs of the exact and obtained solutions on the segment');
legend('Exact fun', 'Boundary value problem, h = pi /22' ,'Boundary value problem, h = pi/44','Cauchy problem, h = pi /22', 'Cauchy problem, h = pi /44');

subplot(1,2,2);
plot(x1, error_mkr1,'m-*');
hold on;
plot(x2, error_mkr2,'g-*');
semilogy(x1, error_runge1,'c-*');
semilogy(x2, error_runge2,'y-*');
xlabel('x');
ylabel('Error');
title('Error graphs on a given segment');
legend('Boundary value problem, h = pi /22' ,'Boundary value problem, h = pi/44','Cauchy problem, h = pi /22', 'Cauchy problem, h = pi /44');

figure;
loglog(outrage_mkr(1:9), outrage_mkr_error(1:9), 'c');
hold on;
loglog(outrage_runge(1:9), outrage_runge_error(1:9), 'g');
loglog(outrage_runge(1:9),outrage_runge(1:9),'r');
xlabel('Perturbation');
ylabel('Error');
title('Dependence of the error rate on the magnitude of the disturbance at a fixed step');
legend('Boundary value problem','Cauchy problem','Bisector');

figure;
loglog(steps, y_h1,'m');
hold on;
loglog(steps, y_h2,'r');
loglog(steps, steps_mkr,'g');
loglog(steps, steps_runge,'c');
title('Dependence of the error rate on the magnitude of the disturbance at a fixed step');
xlabel('h');
ylabel('Error');
legend('h^2', 'h^3', 'MKR','Runge');