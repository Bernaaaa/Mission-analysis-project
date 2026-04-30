function plotter(X, Y, Z, title_str, xlabel_str, ylabel_str, zlabel_str)
if nargin < 4
    title_str = '';
    xlabel_str = '';
    ylabel_str = '';
    zlabel_str = '';   
end
figure()
surf(X, Y, Z)
shading interp
colorbar
xlabel(xlabel_str)
ylabel(ylabel_str)
zlabel(zlabel_str)
title(title_str)
end