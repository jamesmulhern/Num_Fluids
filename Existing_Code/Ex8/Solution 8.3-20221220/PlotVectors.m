function [] = PlotVectors(x,y,u,v,BC)
%PLOTVECTORS Creates a quiver plot of the velocity field
%   at the centers of the cells for pressure. A linear interpolation is
%   employed.

% Append boundary values needed for interpolation
u_with_bc = [ BC.u.left, u, BC.u.right ];
v_with_bc = [ BC.v.bottom; v; BC.v.top ];

% Interpolate to the centers of cells for pressure
ucenter = ( u_with_bc(:,1:end-1) + u_with_bc(:,2:end) ) / 2 ;
vcenter = ( v_with_bc(1:end-1,:) + v_with_bc(2:end,:) ) / 2 ;

% Plot velocity vectors
quiver(x,y,ucenter,vcenter);

xlabel('$x$',"Interpreter","latex");
ylabel('$y$',"Interpreter","latex");



end

