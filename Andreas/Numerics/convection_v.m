function [c_v] = convection_v(u,v,BC,dx,dy)
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

% Get the number of finite volumes for v
[ny,Nx] = size(v);

% Initialize the variables of correct sizes
u_e = NaN(ny,Nx);   u_w = NaN(ny,Nx);
v_e = NaN(ny,Nx);   v_w = NaN(ny,Nx);
v_n = NaN(ny,Nx);   v_s = NaN(ny,Nx);

% Interpolate u and v to the centers of cell faces
u_e(:,1:end-1) = ( u(2:end,:) + u(1:end-1,:) ) / 2; 
u_w(:,2:end  ) = u_e(:,1:end-1);
v_e(:,1:end-1) = ( v(:,2:end) + v(:,1:end-1) ) / 2;
v_w(:,2:end  ) = v_e(:,1:end-1);
v_n(1:end-1,:) = ( v(2:end,:) + v(1:end-1,:) ) / 2;
v_s(2:end  ,:) = v_n(1:end-1,:);

% Boundary values from boundary conditions
u_e(:,end) = ( BC.u.right(2:end) + BC.u.right(1:end-1) ) / 2 ;
u_w(:,1  ) = ( BC.u.left (2:end) + BC.u.left (1:end-1) ) / 2 ;
v_e(:,end) =   BC.v.right;
v_w(:,1  ) =   BC.v.left ;
v_n(end,:) = ( BC.v.top    + v(end,:) ) / 2;
v_s(1  ,:) = ( BC.v.bottom + v(1  ,:) ) / 2;

% Finally, sum-up the convective fluxes through all faces
c_v = (u_e.*v_e-u_w.*v_w)/dx + (v_n.^2-v_s.^2)/dy;

end

