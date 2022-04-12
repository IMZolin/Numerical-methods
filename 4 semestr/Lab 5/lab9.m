error = load('f_error.txt');
eps = load('f_eps.txt');
x = load('f_x.txt');
y = load('f_y.txt');
euler = load('f_euler.txt');
stepLength = load('f_stepLength.txt');
iter = load('f_iter.txt');
wtfError = load('f_wtfError.txt');
wtfDer = load('f_wtfDer.txt');

figure;
plot(x(1:30), y, 'r');
hold on;

finish = 30 + stepLength(1);
plot(x(31:finish), euler(1:stepLength), 'b');
start = finish + 1;
start2 = 1 + stepLength(1);
plot(x(start:end), euler(start2:end), 'g');
xlabel('x');
ylabel('y', 'Rotation', 0);
legend('Функция', 'Эйлер1', 'Эйлер2');

maxError = zeros(1, length(eps));
maxWtfError = zeros(1, length(wtfDer));
for i=1:length(wtfDer)
    start = 29 * (i - 1) + 1;
    finish = 29 * i;
    maxWtfError(i) = max(wtfError(start:finish));
    if i <= length(eps)
        maxError(i) = max(error(start:finish));
    end
end
figure;
loglog(eps, eps, 'r');
hold on;
%loglog(eps, eps.^2, 'm');
loglog(eps, maxError, 'b');
xlabel('Точность');
ylabel('Максимальная ошибка', 'Rotation', 0);

figure;
semilogx(eps, iter);
xlabel('Точность');
ylabel('Итерации', 'Rotation', 0);

figure;
loglog(wtfDer, maxWtfError);
xlabel('Возсущение');
ylabel('Максимальна ошибка', 'Rotation', 0);