function [] = PlotProfiles(xcorners,ycorners,u,v,BC)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


xcenter = ( xcorners(1:end-1) + xcorners(2:end) ) / 2 ;
ycenter = ( ycorners(1:end-1) + ycorners(2:end) ) / 2 ;

Nx = length(xcenter); % number of cells
Ny = length(ycenter);

% Create 2 overlay axes
figure(); % create new figure window
t = tiledlayout(1,1); % since matlab version R2019b


%% Vertical profile of horizontal velocity

color = 'r';
ax1 = axes(t);
hold on;

j = find(xcorners(2:end-1)==0.5);
k = (1:Ny)';
m = k + (j-1)*Ny;

uplot = [ BC.u.bottom(j); u(m) ; BC.u.top(j) ] ;
yplot = [ ycorners(1); ycenter; ycorners(end) ];
plot(ax1, uplot, yplot, color);

ax1.XColor = color;
ax1.YColor = color;
ax1.Box = 'off';
ax1.XAxisLocation = 'top';
xlabel(ax1,'$u(x=0.5,y)$',"Interpreter","latex");
ylabel(ax1,'$y$',"Interpreter","latex");



%% Horizontal profile of vertical velocity

ax2 = axes(t);
hold on;
color = 'k';

k = find(ycorners(2:end-1)==0.5);
j = 1:Nx;
m = k + (j-1)*(Ny-1);

vplot = [ BC.v.left(k), v(m), BC.v.right(k) ] ;
xplot = [ xcorners(1), xcenter, xcorners(end) ] ;
plot(ax2,xplot,vplot,color);

ax2.Color = 'none';
ax2.Box   = 'off';
ax2.YAxisLocation = 'right';
xlabel(ax2,'$x$',"Interpreter","latex");
ylabel(ax2,'$v(x,y=0.5)$',"Interpreter","latex");



%% Extrema

% u_max
[ymax,umax] = get_maximum(yplot,uplot);
plot(ax1,umax,ymax,'*');

% v_min
[xmin,vmin] = get_maximum(xplot,-vplot);
plot(ax2,xmin,-vmin,'ko');

% v_max
[xmax,vmax] = get_maximum(xplot,vplot);
plot(ax2,xmax,vmax,'kd');

% O. Botella and R. Peyret.  Comput. Fluids, 27:421-433, 1998. 
% doi: 10.1016/S0045-7930(98)00002-4
ymax(2) = 0.4581;
umax(2) = 0.2140424;
plot(ax1,umax(2),ymax(2),'b*');

xmin(2) = 0.1896;
vmin(2) = -0.2538030;
plot(ax2,xmin(2),vmin(2),'bo');

xmax(2) = 0.7630;
vmax(2) = 0.1795728;
plot(ax2,xmax(2),vmax(2),'bd');


% Print table to command window
Grid = [sprintf("%u x %u",Nx,Ny); "49 x 49"];
T = table(Grid,umax',ymax',vmax',xmax',vmin',xmin',...
    'VariableNames',["Grid","u_max","y_max","v_max","x_max","v_min","x_min"],...
    'RowName',["present","Botella and Peyret (1998)"]);
disp(T);


end

function [xmax,ymax] = get_maximum(x,y)

    [~,imax] = max(y);
    imax = imax-1:imax+1;
    xmax = x(imax);
    ymax = y(imax);
    p = polyfit(xmax,ymax,2);
    dydx = polyder(p);
    xmax = roots(dydx);
    ymax = polyval(p,xmax);

end