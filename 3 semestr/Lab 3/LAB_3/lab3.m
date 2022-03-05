iterations = load('iterations.txt');
norm = load('actual_error.txt');
nevyaska = load('discrepancy.txt');
accuracy = load('epsilon.txt');

figure;
loglog(accuracy, norm, 'r');
hold on;
loglog(accuracy, nevyaska, 'b');
legend('Фактичекая ошибка','Невязка');
xlabel('Эпсилон');

figure
semilogx(accuracy, iterations, 'bo');
ylabel('Итерации');
xlabel('Эпсилон');