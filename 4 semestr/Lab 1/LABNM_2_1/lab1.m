x = load('f_x.txt');
y = load('f_y.txt');
x_h1 = load('f_x_h1.txt');
y_h1 = load('f_y_h1.txt');
x_h2 = load('f_x_h2.txt');
y_h2 = load('f_y_h2.txt');
x_h3 = load('f_x_h3.txt');
y_h3 = load('f_y_h3.txt');
polynom1 = load('f_polynom1.txt');
polynom2 = load('f_polynom2.txt');
polynom3 = load('f_polynom3.txt');
actual_error1 = load('f_actual_error1.txt');
actual_error2 = load('f_actual_error2.txt');
actual_error3 = load('f_actual_error3.txt');
theoritic_error1 = load('f_theoritic_error1.txt');
theoritic_error2 = load('f_theoritic_error2.txt');
theoritic_error3 = load('f_theoritic_error3.txt');
max_error = load('f_max_error.txt');
nodes = load('f_nodes.txt');

len1 = length(x);
len2 = length(polynom1);


figure;
grid on;
plot(x,y,'r');
hold on;
plot(x,polynom1,'b');
plot(x_h1,y_h1, 'm*');
xlim([0.01 1]);
ylim([-1 1]);
title('4 nodes');
xlabel('X');
ylabel('Polynom');

figure;
grid on;
plot(x,y,'r');
hold on;
plot(x,polynom2,'b');
plot(x_h2,y_h2, 'm*');
xlim([0.01 1]);
ylim([-1 1]);
title('5 nodes');
xlabel('X');
ylabel('Polynom');

figure;
grid on;
plot(x,y,'r');
hold on;
plot(x,polynom3,'b');
plot(x_h3,y_h3, 'm*');
xlim([0.01 1]);
ylim([-1 1]);
title('6 nodes');
xlabel('X');
ylabel('Polynom');


figure;
plot(x,actual_error1,'r');
hold on;
%grid on;
%plot(x_h1,a_e_n1, 'm*');
plot(x,actual_error2,'b');
plot(x,actual_error3, 'g');
plot(x,theoritic_error1,'m');
xlim([0.01 1]);
ylim([-1 1]);
legend('4 nodes', '5 nodes', '6 nodes','th_error 4 nodes');
title('Actual error of the 4,5,6 nodes');
xlabel('X');
ylabel('Actual error');

figure;
semilogy(x,theoritic_error1,'r');
hold on;
grid on;
semilogy(x,theoritic_error2,'b');
semilogy(x,theoritic_error3, 'g');
xlim([0.01 1]);
ylim([0 1]);
legend('4 nodes', '5 nodes', '6 nodes');
%loglog(x, actual_error,'b');
%semilogx(x, theoritic_error,'b');
title('Theoritic error of the 4,5,6 nodes');
xlabel('X');
ylabel('Actual error');


figure;
plot(nodes, max_error);
hold on;
grid on;
title('Dependence of the max error on the number of nodes');
xlabel('Nodes');
ylabel('Max error');