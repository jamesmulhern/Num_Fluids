% Driver script for the computation of an incompressible flow in a  
% lid-driven cavity with the Finite Volume Method (FVM)
%
%
%% Clear workspace

clearvars;
close all;
clc


%% Parameters

% Problem formulation
Re = 50;

% Spatial discretization
Lx = 1;  % width
Ly = 1;  % height
Nx = 20; % number of cells in x-direction
Ny = 20; % number of cells in y-direction

% Time discretization
dt = 0.01;
t_end = 10;

% Visualization
t_plot = 0.1;


%% Mesh

% Dependent parameters
Nx_lines = Nx+1; % number of lines
Ny_lines = Ny+1;
dx = Lx/Nx; % size of the finite volumes
dy = Ly/Ny;

% Cell edges
xcorners =  0:dx:Lx  ;
ycorners = (0:dy:Ly)';

% Cell centers
xcenter = (xcorners(1:end-1)+xcorners(2:end)) / 2;
ycenter = (ycorners(1:end-1)+ycorners(2:end)) / 2;


%% Initial condition

% Let us test the convection flux in a simple shear flow
u = zeros(Ny,Nx-1);
v = zeros(Ny-1,Nx);


%% Boundary conditions

% u                             v
BC.u.left   = zeros(Ny,1   );   BC.v.left   = zeros(Ny-1,1 );
BC.u.right  = zeros(Ny,1   );   BC.v.right  = zeros(Ny-1,1 );
BC.u.top    = ones (1 ,Nx-1);   BC.v.top    = zeros(1   ,Nx);
BC.u.bottom = zeros(1 ,Nx-1);   BC.v.bottom = zeros(1   ,Nx);


%% LHS matrices

[Lap_u,b_u] = Laplacian_u(Nx,Ny,BC,dx,dy);
[Lap_v,b_v] = Laplacian_v(Nx,Ny,BC,dx,dy);
[Lap_p    ] = Laplacian_p(Nx,Ny,   dx,dy); 

Nu = (Nx-1)*Ny; % number of cells for u
Nv = Nx*(Ny-1); % number of cells for v

% Step 1: transport of momentum
LHS_u = speye(Nu,Nu) - dt/Re *Lap_u ;
LHS_v = speye(Nv,Nv) - dt/Re *Lap_v ;

% Pre-compute matrix decomposition
dLHS_u = decomposition(LHS_u);
dLHS_v = decomposition(LHS_v);
dLHS_p = decomposition(Lap_p);


%% Time-stepping loop

for t = 0:dt:t_end
    
    % Step 1: Transport of momentum
    c_u = convection_u(u,v,BC,dx,dy);
    c_v = convection_v(u,v,BC,dx,dy);
    
    rhs_u = u(:) - dt*c_u(:) - dt/Re*b_u; 
    rhs_v = v(:) - dt*c_v(:) - dt/Re*b_v; 
    
    u = dLHS_u \ rhs_u ;
    v = dLHS_v \ rhs_v ;
    
    u = reshape(u,Ny,Nx-1); % reshape the vectors of nodal values to the
    v = reshape(v,Ny-1,Nx); % shape of the mesh. Useful for the next step.
    
    % Step 2: Pressure
    u_with_bc = [ BC.u.left, u, BC.u.right ] ; % use boundary conditions 
    v_with_bc = [ BC.v.bottom; v; BC.v.top ] ; % for computation of deriv-
                                               % atives
    dudx = diff(u_with_bc, 1, 2) / dx ; % x-derivative of u
    dvdy = diff(v_with_bc, 1, 1) / dy ; % y-derivative of v
                                 % __ ->    
    rhs_p = (dudx + dvdy) / dt ; % \/ u  : /_\ t
    p     = dLHS_p \ rhs_p(:) ;
    p     = reshape(p, Nx, Ny);
    
    % Step 3: Correct velocity field
    u = u - dt * diff(p,1,2)/dx ;
    v = v - dt * diff(p,1,1)/dy ;
    
    % Visualize velocity field
    if mod(t,t_plot)==0
        u_with_bc = [ BC.u.left, u, BC.u.right ];
        v_with_bc = [ BC.v.bottom; v; BC.v.top ];
        ucenter = ( u_with_bc(:,1:end-1) + u_with_bc(:,2:end) ) / 2 ;
        vcenter = ( v_with_bc(1:end-1,:) + v_with_bc(2:end,:) ) / 2 ;
        quiver(xcenter,ycenter,ucenter,vcenter);
        xlabel('x');
        ylabel('y');
        title(sprintf('t = %g',t));
        drawnow;
    end
    
end