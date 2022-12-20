function [c_u] = convection_u(u,v,BC,dx,dy)
%CONVECTION_U Computes the convection flux of x-momentum
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
% We have Ny*(Nx-1) cells for u.

% Nodal values
u_E = [ u(:,2:end) , BC.u.right   ];
u_W = [ BC.u.left  , u(:,1:end-1) ];
u_N = [ u(2:end,:) ; 2*BC.u.top - u(end,:)   ];
u_S = [ 2*BC.u.bottom - u(1,:); u(1:end-1,:) ];
v_N = [ v(2:end,:) ; BC.v.top     ];
v_S = [ BC.v.bottom; v(1:end-1,:) ];

% Interpolate u and v to the centers of cell faces
u_e = ( u_E + u ) / 2 ; 
u_w = ( u_W + u ) / 2 ;
v_n = ( v_N + v ) / 2 ;
v_s = ( v_S + v ) / 2 ;

% Upwind scheme
uu_e = u_e .* ( u  .*(u_e>0) + u_E.*(u_e<0) ) ;
uu_w = u_w .* ( u_W.*(u_w>0) + u  .*(u_w<0) ) ;
uv_n = v_n .* ( u  .*(v_n>0) + u_N.*(v_n<0) ) ;
uv_s = v_s .* ( u_S.*(v_s>0) + u  .*(v_s<0) ) ;

% Finally, sum-up the convective fluxes through all faces
c_u = (uu_e-uu_w)/dx + (uv_n-uv_s)/dy;


end

