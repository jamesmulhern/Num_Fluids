function [] = PlotVectors(x,y,u,v)
%PLOTVECTORS Creates a quiver plot of the velocity field
%   at the centers of the cells for pressure. A linear interpolation is
%   employed.

Nplot = 20;

[nx,~] = size(u);
[~,ny] = size(v);
Nplot = min([Nplot,nx,ny]);
nskipx = round(nx/Nplot);
nskipy = round(ny/Nplot);

% Interpolate to the centers of cells for pressure
iuy = 2:nskipy:ny-1;
ivx = 2:nskipx:nx-1;
ucenter = ( u(iuy,1:nskipx:end-1) + u(iuy,2:nskipx:end) ) / 2 ;
vcenter = ( v(1:nskipy:end-1,ivx) + v(2:nskipy:end,ivx) ) / 2 ;

% Plot velocity vectors
quiver(x(ivx),y(iuy),ucenter,vcenter);

% Axis labels
axis([0 1 0 1]);
xlabel('$x$',"Interpreter","latex");
ylabel('$y$',"Interpreter","latex");



end

