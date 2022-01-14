function [] = draw(A,b,n,levels)
Tol = @(A1,b1,x) min(rad(b1) - mag(mid(b1) - A1 * x));
x = linspace(-5, 5, n);
y = linspace(-5, 5, n);
[xx, yy] = meshgrid(x, y);
zz = zeros([size(xx, 1), size(xx, 2)]);
for i=1:size(xx, 1)
    for j=1:size(yy,1)
        X = [xx(i, j); yy(i, j)];
        zz(i, j) = Tol(A, b, X);
    end
end
figure
[~,~,~,~,~] = EqnTolR2(inf(A), sup(A), inf(b), sup(b));
hold on
contour(xx, yy, zz, levels)
hold on
%%%   Tolmaximizing   %%%
[tolMax,argMax,~,~] = tolsolvty(inf(A), sup(A), inf(b), sup(b));
tolMax
plot(argMax(1), argMax(2), 'k+');
hold on 
end