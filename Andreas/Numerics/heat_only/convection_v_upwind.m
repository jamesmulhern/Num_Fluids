function [c_v] = convection_v_upwind(u,v,BC,dx,dy)
%CONVECTION_V Computes the convection flux of y-momentum
%
%   using the Finite Volume Method (FVM) on a staggered grid. Linear
%   interpolation is employed where needed.
%
%   Inputs:
%    * u (Ny  ,Nx-1), 
%    * v (Ny-1,Nx  ) - velocity components on their internal nodes
%    * BC - structure containing boundary values of u and v:
%       |
%       |__                 shape  
%       |   u
%       |   |__ left        (Ny,1   ) 
%       |   |__ right       (Ny,1   )
%       |   |__ bottom      (1 ,Nx-1)
%       |   |__ top         (1 ,Ny-1)
%       |__ 
%           v
%           |__ left        (Ny-1,1 )
%           |__ right       (Ny-1,1 )
%           |__ bottom      (1   ,Nx)
%           |__ top         (1   ,Nx)
%
%    * dx,dy - grid spacings in x and y
%
%   Output: vector of the total convective flux for each cell


% First compute the velocities at the centers of edges of each cell.
% We have (Ny-1)*Nx cells for v.

% Initiate arrays
[ny,Nx] = size(v); % number of cells for v
u_e = NaN(ny,Nx);   
u_w = NaN(ny,Nx);

% Nodal values
v_E = [ v(:,2:end) , 2*BC.v.right - v(:,end) ];
v_W = [ 2*BC.v.left - v(:,1) , v(:,1:end-1)  ];
v_N = [ v(2:end,:) ; BC.v.top     ];
v_S = [ BC.v.bottom; v(1:end-1,:) ];

% Interpolate u and v to the centers of cell faces
u_e(:,1:end-1) = ( u(2:end,:) + u(1:end-1,:) ) / 2; 
u_w(:,2:end  ) = u_e(:,1:end-1);
v_n            = ( v_N + v ) / 2;
v_s            = ( v_S + v ) / 2;

% Boundary values from boundary conditions
u_e(:,end) = ( BC.u.right(2:end) + BC.u.right(1:end-1) ) / 2 ;
u_w(:,1  ) = ( BC.u.left (2:end) + BC.u.left (1:end-1) ) / 2 ;

% Fluxes with upwind scheme
uv_e = u_e .* ( v  .*(u_e>0) + v_E.*(u_e<0) );
uv_w = u_w .* ( v_W.*(u_w>0) + v  .*(u_w<0) );
vv_n = v_n .* ( v  .*(v_n>0) + v_N.*(v_n<0) );
vv_s = v_s .* ( v_S.*(v_s>0) + v  .*(v_s<0) );

% Finally, sum-up the convective fluxes through all faces
c_v = (uv_e-uv_w)/dx + (vv_n-vv_s)/dy;


end

