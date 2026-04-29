function plotter(X, Y, Z, title_str)
if nargin < 4
    title_str = '';
end
figure()
surf(X, Y, Z)
shading interp
colorbar
xlabel('e')
ylabel('r_a')
zlabel('\Delta V')
title(title_str)
end