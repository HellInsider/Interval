%% init
Tol = @(A,b,x) min(rad(b) - mag(mid(b) - A * x));

A = [infsup(0, 2), infsup(1, 3); 1, -infsup(2, 4); infsup(1, 3), 0];
b = [infsup(3, 7); 0; infsup(5, 7)];

[tolMax,argMax,envs,ccode] = tolsolvty(inf(A), sup(A), inf(b), sup(b));

n = 100;
levels = 30;
draw(A,b,n,levels)

%% %   vector b correction   %%%
e = [infsup(-1, 1); infsup(-1, 1); infsup(-1, 1)];
C = 1.5 * abs(tolMax);
b1 = b + C * e;
draw(A,b1,n,levels)

% variability
[tolMax1,argMax1,~,~] = tolsolvty(inf(A), sup(A), inf(b1), sup(b1));
ive1 = ive(A, b1);
rve1 = rve(A, tolMax1);
iveBox = [midrad(argMax1(1), ive1);midrad(argMax1(2), ive1)];
rveBox = [midrad(argMax1(1), rve1);midrad(argMax1(2), rve1)];
plotintval([iveBox, rveBox], 'n');

%% %   matrix A correction   %%%
koef = 1.5;
b2 = [infsup(-1, 4); infsup(-2, 2); infsup(2, 7)];
[tolMax2,argMax2,~,~] = tolsolvty(inf(A), sup(A), inf(b2), sup(b2));
E = 1.036 * [argMax2(1) * 0.23 0.96455;0 0.96455;argMax2(1) * 0.23 0];
A1 = infsup(inf(A) + E, sup(A) - E);
draw(A1,b2,n,levels);
% variability
[tolMax2,argMax2,~,~] = tolsolvty(inf(A1), sup(A1), inf(b2), sup(b2));
ive2 = ive(A1, b2);
rve2 = rve(A1, tolMax2);
iveBox2 = [midrad(argMax2(1), ive2);midrad(argMax2(2), ive2)];
rveBox2 = [midrad(argMax2(1), rve2);midrad(argMax2(2), rve2)];
plotintval([iveBox2, rveBox2], 'n');


%% УПРАВЛЕНИЕ ПОЛОЖЕНИЕМ МАКСИМУМА

iterations = 10;
figure
A2 = A;
for i = 1:iterations
    A2 = A2 ./ 2;
    [~,argMax,~,~] = tolsolvty(inf(A2), sup(A2), inf(b), sup(b));
    plot(argMax(1), argMax(2), 'k+');
    hold on
end
grid on

figure 
A2 = A;
line1 = [1, 1; 0, 0; 0, 0]
line2 = [0, 0; 1, 1; 0, 0]
line3 = [0, 0; 0, 0; 1, 1]
for i = 1:iterations
    C = rad(A2)./2 .*line1;
    A2 = infsup(inf(A2)+C, sup(A2)-C);
    [~,argMax,~,~] = tolsolvty(inf(A2), sup(A2), inf(b), sup(b));
    plot(argMax(1), argMax(2), 'k+');
    hold on

    C = rad(A2)./2 .*line2;
    A2 = infsup(inf(A2)+C, sup(A2)-C);
    [~,argMax,~,~] = tolsolvty(inf(A2), sup(A2), inf(b), sup(b));
    plot(argMax(1), argMax(2), 'k+');
    hold on

    C = rad(A2)./2 .*line3;
    A2 = infsup(inf(A2)+C, sup(A2)-C);
    [~,argMax,~,~] = tolsolvty(inf(A2), sup(A2), inf(b), sup(b));
    plot(argMax(1), argMax(2), 'k+');
    hold on
    A2
end
grid on