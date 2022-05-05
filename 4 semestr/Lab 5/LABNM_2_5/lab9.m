error = load('f_maxError.txt');
eps = load('f_eps.txt');
x = load('f_x.txt');
y = load('f_y.txt');
euler = load('f_euler.txt');
stepLength = load('f_stepLength.txt');
iter = load('f_iter.txt');
wtfError = load('f_wtfError.txt');
wtfDer = load('f_wtfDer.txt');
eulerError = load('f_error.txt');

figure;
plot(x(1:30), y, 'r');
hold on;

finish = 30 + stepLength(1);
plot(x(31:finish), euler(1:stepLength), 'g-*');
%plot(x(1), euler(1), 'g*');
start = finish + 2;
start2 = 2 + stepLength(1);
plot(x(start:end), euler(start2:end), 'b--');

xlabel('x');
ylabel('y', 'Rotation', 0);
title('Graphs of exact and obtained solutions on a segment with two different steps');
legend('Function', 'Runge-Kutta1', 'Runge-Kutta2');

figure;
finish=30+stepLength(1);
plot(x(31:finish), eulerError(1:stepLength),'b');
hold on;
start = finish + 2;
start2 = 2 + stepLength(1);
plot(x(start:end), eulerError(start2:end), 'k');
legend('Step 3\pi /46', 'Step 3\pi/92');
title('Error graph on a given segment');
xlabel('x');
ylabel('Error', 'Rotation', 0);

maxError = zeros(1, length(eps));
maxWtfError = zeros(1, length(wtfDer));
for i=1:length(eps)
    start = 29 * (i - 1) + 1;
    finish = 29 * i;
    maxError(i) = max(error(start:finish));
    %maxWtfError(i) = max(wtfError(start:finish));
    if i <= length(wtfDer)
         maxWtfError(i) = max(wtfError(start:finish));
        %maxError(i) = max(error(start:finish));
    end
end
figure;
loglog(eps, eps, 'r');
hold on;
%loglog(eps, eps.^2, 'm');
loglog(eps, maxError, 'b');
loglog(eps, 29.*eps, 'm');
xlabel('Epsilon');
ylabel('Max error', 'Rotation', 0);
title('Graph of the dependence of the actual accuracy on the specified accuracy');

figure;
semilogx(eps, iter);
xlabel('Epsilon');
ylabel('Iterations', 'Rotation', 0);
title('Graph of the dependence of the number of iterations on the specified accuracy');

figure;
loglog(wtfDer,maxWtfError);
xlabel('Возмущение');
ylabel('Max error', 'Rotation', 0);
title('Graph of the dependence of the error rate on the magnitude of the disturbance at a fixed accuracy');