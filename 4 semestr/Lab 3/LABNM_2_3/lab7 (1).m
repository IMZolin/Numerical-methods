iter = load('f_iterrations.txt');
eps = load('f_epsilon.txt');
error = load('f_error.txt');

figure;
loglog(eps, error);
hold on;
loglog(eps, eps, 'r');
xlabel('Точность');
ylabel('Ошибка', 'Rotation', 0);
legend('График зависимости ошибки от точности', 'биссектриса','Location', 'North');

figure;
semilogx(eps, iter, 'o');
xlabel('Точность');
ylabel('Итерации', 'Rotation', 0);