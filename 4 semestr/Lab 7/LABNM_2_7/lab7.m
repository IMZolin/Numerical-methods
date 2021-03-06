y1 = load('files\f_y1.txt.');
sol1 = load('files\f_sol1.txt');
x1 = load('files\f_x1.txt');
y2 = load('files\f_y1.txt.');
sol2 = load('files\f_sol2.txt');
x2 = load('files\f_x2.txt');
error1 = load('files\f_error1.txt');
error2 = load('files\f_error2.txt');
perturb_error = load('files\f_perturb_error.txt');
perturbation = load('files\f_perturbation.txt');
actual_error = load('files\f_actual_error.txt');
steps = load('files\f_steps.txt');

figure;
subplot(1,2,1);
plot(x1,y1,'r');
hold on;
plot(x1, sol1,'g-*');
plot(x2, sol2,'b-*');
title('Graphs of the exact and obtained solutions on the segment');
xlabel('x');
ylabel('y');

subplot(1,2,2);
plot(x1,error1,'g-*');
hold on;
plot(x2,error2, 'b-*');
title('Error graphs on a given segment');
xlabel('x');
ylabel('Error');

figure;
loglog(steps, actual_error,'b-*');
hold on;
loglog(steps, steps.^2,'r');
legend('actual error inf norm','h^2');
title('Graph of the dependence of the infinite norm of the actual accuracy on the magnitude of the step');
xlabel('h');
ylabel('Error');

figure;
loglog(perturbation, perturb_error,'b-*');
xlabel('Perturbation');
ylabel('Error');
title('Graph of the dependence of the error rate on the magnitude of the disturbance at a fixed step');

