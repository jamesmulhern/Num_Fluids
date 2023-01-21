function [c_T] = convection_T_upwind_new(u,v,BC,dx,dy)
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

% Interpolate u and v to the centers of cell faces
u_e = ( u_E + u ) / 2 ; 
u_w = ( u_W + u ) / 2 ;
v_n = [ (    v ( : ,  2:end) +    v( : ,   1:end-1) ) / 2 ;
        ( BC.v.top   (2:end) + BC.v.top   (1:end-1) ) / 2 ];
v_s = [ ( BC.v.bottom(2:end) + BC.v.bottom(1:end-1) ) / 2 ;
        (    v ( : ,  2:end) +    v ( : ,  1:end-1) ) / 2 ];

% Nodal Temperature
T_E = [ TM(:,2:end) , BC.T.right   ];
T_W = [ BC.T.left  , TM(:,1:end-1) ];
T_N = [ TM(2:end,:) ; 2*BC.T.top - TM(end,:)   ];
T_S = [ 2*BC.T.bottom - TM(1,:); TM(1:end-1,:) ];

% Interpolate Temperature
T_e = ( T_E + TM ) / 2 ;
T_w = ( T_W + TM ) / 2 ;
T_n = ( T_N + TM ) / 2 ;
T_s = ( T_S + TM ) / 2 ;


% Upwind scheme
u_eT_e = u_e .* ( TM  .*(u_e>0) + T_E.*(u_e<0) ) ;
u_wT_w = u_w .* ( T_W.*(u_w>0) + TM  .*(u_w<0) ) ;
v_nT_n = v_n .* ( TM  .*(v_n>0) + T_N.*(v_n<0) ) ;
v_sT_s = v_s .* ( T_S.*(v_s>0) + TM  .*(v_s<0) ) ;

% Finally, sum-up the convective fluxes through all faces
c_T = (u_eT_e-u_wT_w)/dx + (v_nT_n-v_sT_s)/dy;


end

