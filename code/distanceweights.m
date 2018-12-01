close all;
Nx = 64;
Ny = 40;

[X Y]= meshgrid(linspace(-1.6,1.6,Nx), linspace(-1,1,Ny));

Z= 1./(sqrt((X).^2+(Y).^2))+0.1;
%Z= exp(-sqrt(X.^2+Y.^2).^2/0.1)+0.1;

surf(X, Y, Z);
shading('flat');
colormap('parula');