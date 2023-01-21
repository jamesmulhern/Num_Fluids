% Driver script for the computation of an incompressible flow in a  
% lid-driven cavity with the Finite Volume Method (FVM)
%
%
%% Clear workspace
clearvars;
close all;
clc


%% Parameters
Re = 50;                % Reynolds number     

% Spatial discretization
Lx = 1;  % width
Ly = 1;  % height
Nx = 20; % number of cells in x-direction
Ny = 20; % number of cells in y-direction


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
xcenter = 1/2*(xcorners(1:end-1)+xcorners(2:end));
ycenter = 1/2*(ycorners(1:end-1)+ycorners(2:end));


%% Initial condition
% Let us test the convection flux in a simple shear flow
u = ones(Ny,Nx-1).*ycenter;
v = zeros(Ny-1,Nx);


%% Boundary conditions
% u                             v
BC.u.left   = ycenter       ;   BC.v.left   = zeros(Ny-1,1 );
BC.u.right  = ycenter       ;   BC.v.right  = zeros(Ny-1,1 );
BC.u.top    = ones (1 ,Nx-1);   BC.v.top    = zeros(1   ,Nx);
BC.u.bottom = zeros(1 ,Nx-1);   BC.v.bottom = zeros(1   ,Nx);


%% Test convective flux
c_u = convection_u(u,v,BC,dx,dy);

% Test the output
if all(size(c_u)==size(u))
    disp('Good job, the size of the output seems correct!');
else
    error('The size of your output does not match the number of cells for u');
end
if max(abs(c_u))==0
    disp('The values seem right. The convective flux in a shear flow is zero.');
else
    fprintf('The maximum convective flux is: %g',max(abs(c_u)));
    error('Hm, the values seem wrong. The convective flux shuld be zero everywhere.')
end