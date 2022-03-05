error = load('f_err.txt');
ermit = load('f_p.txt');
x = load('f_x.txt');
y = load('f_y.txt');
gridN = load('f_gridN.txt');
gridX = load('f_gridX.txt');
gridY = load('f_gridY.txt');
theoryError = load('f_theory.txt');
n = x(1);
x = x(2:end);
start = 1;
ending = n;
gridStart = 1;
gridEnd = gridN(1);
color = ['r', 'm', 'b'];
figure;
plot(x, y, 'm');
title('График функции');
grid on;
for i=1:3
figure;
plot(x, ermit(start:ending), color(i));
hold on;
start = ending + 1;
ending = ending + n;
titl = strcat('График полинома с ', num2str(gridN(i)),' узлами');
title(titl);
if (i ~= 1)
gridStart = gridEnd + 1;
gridEnd = gridEnd + gridN(i);
end
plot(gridX(gridStart:gridEnd), gridY(gridStart:gridEnd), 'k*');
for j=gridStart:gridEnd
plot([gridX(j), gridX(j)], [gridY(j), gridY(j) + 3], 'k--', 'LineWidth', 0.5);
text(gridX(j) - 0.1, gridY(j) + 3, sprintf('(%.2f, %.2f)', gridX(j), gridY(j)));
end
xlabel('x');
ylabel('y', 'Rotation',0);
grid on;
end
start = 1;
ending = 30;
figure;
title('График фактической ошибки от числа узлов');
for i=1:3
semilogy(x, error(start: ending), strcat(color(i), 'o'));
hold on;
start = ending + 1;
ending = ending + 30;
end
loglog(x, theoryError, 'go');
xlabel('x');
ylabel('Ошибка', 'Rotation',0);
legend('4 узла', '6 узлов', '8 узлов', 'теоретическая ошибка для 4 узлов', 'Location', 'NorthEastOutside');
maxError = zeros(1, length(gridN));
start = 1;
ending = 30;
for i=1:length(gridN)
maxError(i) = max(error(start:ending));
start = ending + 1;
ending = ending + 30;
end
figure;
semilogy(gridN, maxError);
xlabel('узлы');
ylabel('ошибка', 'Rotation',0);
title('График зависимости макс ошибки от числа узлов');