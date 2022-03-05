f = @(x) log10(x) + 7/(2*x+6);
%x=load('x.txt');
%lagrange =load('coefficients.txt');
x10 = [0.9944210634216271 0.9505337813261752 0.8666548714109941 0.7502298445041029 0.6115931672713687 0.4630509233581410 0.3177884639640546 0.1887000102397788 0.0872440978404414 0.0224264600916150];
lagrange10 =[0.8737924116564679 0.8639237821111662 0.8430213575293520 0.8084704676978695 0.7555640263085135 0.6762982703828422 0.5570577039394093 0.3733979119905708 0.0744332595448101 -0.4912292938159535];
x7 =[0.9890848353201114 0.9046020904940162 0.7502298445041029 0.5526338005665857 0.3459459843610037 0.1658689014122528 0.0435084154665060 ];
lagrange7 =[0.8726277689126107 0.8528357374649322 0.8084704676978695 0.7276220131143313 0.5850502850277922 0.3253066406782242 -0.2114381270120695];
x3 =[0.9505337813261752 0.6115931672713688 0.1887000102397788];
lagrange3 =[0.8639237821111662 0.7555640263085135 0.3733979119905708];

p10 = @(x) polyval(lagrange10,x);
p_x10 = p10(x10);
figure;
a=0;
b=1;
fplot(f,[a,b],'b')
hold on;
loglog(x10,lagrange10);
loglog(x7,lagrange7);
loglog(x3,lagrange3,'g');
legend('function','Lagrange: nodes 10','Lagrange: nodes 7','Lagrange: nodes 3');