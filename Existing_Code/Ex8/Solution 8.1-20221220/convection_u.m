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

% Get the number of finite volumes for u
[Ny,nx] = size(u);

% Initialize the variables of correct sizes
u_e = NaN(Ny,nx);   u_w = NaN(Ny,nx);
u_n = NaN(Ny,nx);   u_s = NaN(Ny,nx);
v_n = NaN(Ny,nx);   v_s = NaN(Ny,nx);

% Interpolate u and v to the centers of cell faces
u_e(:,1:end-1) = ( u(:,2:end) + u(:,1:end-1) ) / 2; 
u_w(:,2:end  ) = u_e(:,1:end-1);
u_n(1:end-1,:) = ( u(2:end,:) + u(1:end-1,:) ) / 2;
u_s(2:end  ,:) = u_n(1:end-1,:);
v_n(1:end-1,:) = ( v(:,2:end) + v(:,1:end-1) ) / 2;
v_s(2:end  ,:) = v_n(1:end-1,:);

% Boundary values from boundary conditions
u_e(:,end) = ( BC.u.right + u(:,end) ) / 2 ;
u_w(:,1  ) = ( BC.u.left  + u(:,1  ) ) / 2 ;
u_n(end,:) = BC.u.top;
u_s(1  ,:) = BC.u.bottom;
v_n(end,:) = ( BC.v.top   (2:end) + BC.v.top   (1:end-1) ) / 2;
v_s(1  ,:) = ( BC.v.bottom(2:end) + BC.v.bottom(1:end-1) ) / 2;

% Finally, sum-up the convective fluxes through all faces
c_u = (u_e.^2-u_w.^2)/dx + (u_n.*v_n-u_s.*v_s)/dy;

end

