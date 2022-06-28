x = load('files\f_x.txt');
x1 = load('files\f_x1.txt');
x2 = load('files\f_x2.txt');
y = load('files\f_y.txt');
runge1 = load('files\f_runge1.txt');
runge2 = load('files\f_runge2.txt');
n = load('files\f_stepLength.txt');
error1 = load('files\f_error1.txt');
error2 = load('files\f_error2.txt');
iter = load('files\f_iter.txt');
eps = load('files\f_eps.txt');
max_error = load('files\f_maxError.txt');
outrage = load('files\f_outrage.txt');
outrage_error = load('files\f_outrage_error.txt');

figure;
plot(x, y, 'r');
hold on;
plot(x1(1:n(1)),runge1(1:n(1)),'b-o');
plot(x2(1:n(2)),runge2(1:n(2)),'g-*');
xlabel('x');
ylabel('y', 'Rotation', 0);
title('Graphs of exact and obtained solutions on a segment with two different steps');
legend('Function', 'Runge-Kutta1', 'Runge-Kutta2');
figure;
plot(x1, error1);
hold on;
plot(x2, error2);
legend('Step 1: 3\pi /46', 'Step 2: 3\pi/92');
title('Error graph on a given segment');

figure;
subplot(1,2,1);
loglog(eps, max_error);
hold on;
loglog(eps, eps, 'r');
loglog(eps, eps*30, 'b');
title('Dependence of the actual accuracy on the specified accuracy');
subplot(1,2,2);
semilogx(eps,iter);
title('Dependence of the number of iterations on the specified accuracy');

figure;
loglog(outrage,outrage_error);
title('Dependence of the error rate on the magnitude of the disturbance at a fixed accuracy');


