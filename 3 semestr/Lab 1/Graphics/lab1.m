%алгебраическая функция
%2*x^4-x^2-10
f = @(x) 2*x.^4-x.^2-10;
%трансцедентная функция
%x+lg(1+x)-1.5
g = @(x) x+log10(1+x)-1.5;

%Итерации
f_bisectional_iterations = [3, 6, 9, 13, 16, 19];
f_fixedpoint_iterations = [4, 5, 6, 7, 7, 8];
g_bisectional_iterations = [3, 6, 9, 13, 16, 19];
g_fixedpoint_iterations = [45, 70, 95, 120, 146, 171];

%Погрешность
errors = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6];

f_x_roots = zeros(6, 1);
f_y_roots = zeros(6, 1);
g_x_roots = zeros(6, 1);
g_y_roots = zeros(6, 1);
f_exitflag = zeros(6, 1);
g_exitflag = zeros(6, 1);
f_output = struct();
g_output = struct();
fzero_f_iterations = zeros(6, 1);
fzero_g_iterations = zeros(6, 1);

for i=1:6
    options = optimset('TolX', errors(i));
    [f_x_roots(i), f_y_roots(i), f_exitflag(i), f_output] = fzero(f, [1, 2], options);
    [g_x_roots(i), g_y_roots(i), g_exitflag(i), g_output] = fzero(g, [0.7, 1.7], options);
    fzero_f_iterations(i) = f_output.iterations;
    fzero_g_iterations(i) = g_output.iterations;
end

f_y_roots = abs(f_y_roots);
g_y_roots = abs(g_y_roots);

fx_approximations = [1.500000000000000  1.75000000000000, 1.625000000000000, 1.562500000000000];
fy_approximations = f(fx_approximations);
gx_approximations = [1.200000000000000, 0.950000000000000, 1.075000000000000, 1.137500000000000];
gy_approximations = g(gx_approximations);

fx_bisectional_roots = [1.562500000000000, 1.585937500000000, 1.581054687500000, 1.581115722656250, 1.581138610839844, 1.581139564514160];
fx_fixedpoint_roots = [1.580497205008714, 1.581105506152811, 1.581137118936483, 1.581138742270902, 1.581138742270902, 1.581138825577894];
gx_bisectional_roots = [1.137500000000000, 1.160937500000000, 1.163867187500000, 1.164660644531250, 1.164622497558594, 1.164617729187012];
gx_fixedpoint_roots = [1.172551052819379, 1.163785762849198, 1.164677445442137, 1.164586702905622, 1.164594305740714, 1.164595163340671];

fy_fixedpoint_roots = abs(f(fx_fixedpoint_roots));
fy_bisectional_roots = abs(f(fx_bisectional_roots));
gy_fixedpoint_roots = abs(g(gx_fixedpoint_roots));
gy_bisectional_roots = abs(g(gx_bisectional_roots));

%Графики функций f(x) и g(x)
figure;
x=-2:0.0001:2;%диапазон и шаг для х
plot(x,f(x),'b'); %строю функцию, синей сплошной линией
grid on;
xlabel('x'); %оси
ylabel('y', 'Rotation',0);
title('График алгебраичсекой функции');
figure;
x=-2:0.0001:3;%диапазон и шаг для х
plot(x, g(x), 'k');
grid on;
xlabel('x'); %оси
ylabel('y', 'Rotation',0);
title('График трансцедннтной функции');

graphic_interpretation(1.4, 1.8, f, fx_approximations, fy_approximations,'метода половинного деления');
graphic_interpretation(0.9, 1.3, g, gx_approximations, gy_approximations, 'метода половинного деления');

iterations_dependences(errors, f_fixedpoint_iterations, f_bisectional_iterations, fzero_f_iterations,'2*x^4-x^2-10');
iterations_dependences(errors, g_fixedpoint_iterations, g_bisectional_iterations, fzero_f_iterations, 'x+log10(1+x)-1.5');

function_value_dependences(errors, fy_fixedpoint_roots, fy_bisectional_roots, f_y_roots, '2*x^4-x^2-10');
function_value_dependences(errors, gy_fixedpoint_roots, gy_bisectional_roots, g_y_roots, 'x+log10(1+x)-1.5');

%Графическая интерпритация метода половинного деления
function graphic_interpretation(left, right, func, roots, froots, method_name)
    x = left:0.1:right;
    figure;
    grid on;
    hold all;
    plot(x, func(x), 'm');
    xlabel('x'); %оси
    ylabel('y', 'Rotation',0);
    title(sprintf('Графическая интерпретация %s', method_name));
    plot([left, right], [0, 0], 'k');
    for i=1:4
        plot(roots(i), froots(i), 'r*');  
    end
    plot([roots(1), roots(1)],[froots(1), 0], 'b', 'LineWidth', 0.2);
    plot([roots(2), roots(2)],[froots(2), 0], 'b', 'LineWidth', 0.2);
    plot([roots(3), roots(3)],[froots(3), 0], 'b', 'LineWidth', 0.2);
    plot([roots(4), roots(4)],[froots(4), 0], 'b', 'LineWidth', 0.2);
end

%Зависимость количества итераций функций от погрешности
function iterations_dependences(errors, fixedpoint_iterations, bisectional_iterations, fzero_iterations, func_name)
    figure;
    semilogx(errors, fixedpoint_iterations, 'r');
    hold on;
    semilogx(errors, bisectional_iterations, 'm');
    semilogx(errors, fzero_iterations, 'b');
    plot(errors, fixedpoint_iterations, 'ro');
    plot(errors, bisectional_iterations, 'mo');
    plot(errors, fzero_iterations, 'bo');
    legend('Метод простых итераций', 'Метод половинного деления', 'Метод fzero');
    xlabel('Погрешность');
    ylabel('Число итераций');
    title(sprintf('Графики зависимости количества итераций %s от погрешности', func_name));
end

%Зависимость модуля значения функций в найденном корне от заданной точности
function function_value_dependences(errors, fixedpoint_iteration_values, bisectional_values, fzero_values, func_name)
    figure;
    loglog(errors, fzero_values, 'm');
    hold on;
    semilogx(errors, bisectional_values, 'b');
    semilogx(errors, fixedpoint_iteration_values, 'r');
    plot(errors, errors, 'k');
    legend('Метод матлаба fzero', 'Метод половинного деления', 'Метод простых итераций', 'Точность', 'Location', 'NorthWest');
    xlabel('Погрешность')
    ylabel('Модуль значения функции в корне');
    title(sprintf('Графики зависимости модуля значения %s в найденном корне от заданной точности', func_name));
end
