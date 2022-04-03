iter = load('f_iterrations.txt');
eps = load('f_epsilon.txt');
error = load('f_error.txt');

figure;
loglog(eps, error);
hold on;
loglog(eps, eps, 'r');
xlabel('��������');
ylabel('������', 'Rotation', 0);
legend('������ ����������� ������ �� ��������', '�����������','Location', 'North');

figure;
semilogx(eps, iter, 'o');
xlabel('��������');
ylabel('��������', 'Rotation', 0);