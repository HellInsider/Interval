A = kinterval([3 5;-1 -3], [4 6;1 1]);
b1 = kinterval([-3;-1], [3;2]);
b2 = kinterval([-3;-1],[4;2]);

[x1, opt1, history1] = subdiff(A, b1);
opt1.tau = 1;
opt1.max_iterations_num = 300;
[x2, opt2, history2] = subdiff(A, b2, opt1);

tol1 = [0 0.5 0.75 0.8889 0.5 0.2 0 -0.75 -0.8235 0;0.3333 0.1667 0 -0.1111 -0.5 -0.6 -0.5 0 0.0588 0.3333];
tol2 = [0 1 0.5 0.2 0 -0.75 -0.8236 0;0.3333 0 -0.5 -0.6 -0.5 0 0.0588 0.3333];

figure
for i=1:size(history1, 2) - 1
    x = history1(:, i);
    plot([inf(x(1)), inf(x(1)), sup(x(1)), sup(x(1)), inf(x(1))], [inf(x(2)), sup(x(2)), sup(x(2)), inf(x(2)), inf(x(2))], ':')
    hold on
    plot(mid(x(1)), mid(x(2)))
    hold on
endfor
i = size(history1, 2);
x = history1(:, i);
plot([inf(x(1)), inf(x(1)), sup(x(1)), sup(x(1)), inf(x(1))], [inf(x(2)), sup(x(2)), sup(x(2)), inf(x(2)), inf(x(2))], 'r')
hold on
plot(mid(x(1)), mid(x(2)))
hold on
plot(tol1(1,:), tol1(2,:), '--')

figure
for i=1:size(history2, 2) - 1
    x = history2(:, i);
    plot([inf(x(1)), inf(x(1)), sup(x(1)), sup(x(1)), inf(x(1))], [inf(x(2)), sup(x(2)), sup(x(2)), inf(x(2)), inf(x(2))], ':')
    hold on
    plot(mid(x(1)), mid(x(2)))
    hold on
endfor
i = size(history2, 2);
x = history2(:, i);
plot([inf(x(1)), inf(x(1)), sup(x(1)), sup(x(1)), inf(x(1))], [inf(x(2)), sup(x(2)), sup(x(2)), inf(x(2)), inf(x(2))], 'r')
hold on
plot(mid(x(1)), mid(x(2)))
hold on
plot(tol2(1,:), tol2(2,:), '--')