
condition_number = load('condition.txt');
accuracy = load('actual_error.txt');
nevyazka = load('discrepancy.txt');

figure;
loglog(condition_number, accuracy, 'r');
hold on;
loglog(condition_number, nevyazka, 'b');
legend('Фактическая точность', 'Невязка', 'Location', 'NorthWest');
title('График зависимости фактической точности и невязки от числа обусловленности');
xlabel('x');
ylabel('y', 'Rotation',0);