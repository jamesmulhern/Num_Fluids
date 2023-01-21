% 
%  Finite Volume method for the transport of a scalar in a given 
%  velocity field
%
%  convection and diffusion
%
%
% definition of the control volumes
%

clearvars;
close all;
clc

% domain size
Lx = 1;
Ly = 1;
%
% number of cells
Nx = 100;
Ny = 100;
% number of lines
Nx_lines = Nx+1;
Ny_lines = Ny+1;
%
dx = Lx/Nx;
dy = Ly/Ny;
%
% cell corners
%
xcorners = 0:dx:Lx;
ycorners = 0:dy:Ly;
%
% cell centres
xcenter = 1/2*(xcorners(1:end-1)+xcorners(2:end));
ycenter = 1/2*(ycorners(1:end-1)+ycorners(2:end));
%
% velocity field at the corners
% ux = x
% uy = -y
%
ux = xcorners * 0;
uy = -ycorners * 0;
%
% diffusivity coefficient
%
kappa = 0.01;
%
% left boundary condition
%
Tw = 1 - ycenter ;
%
% we are going to solve the problem in the form of a linear system
%
%      A T = B
%
% matrices initialization
%
M = Nx*Ny;
A = zeros(M,M);
B = zeros(M,1);
%

% these are help coefficients 
%
AP = zeros(M,1);
AE = AP;
AW = AP;
AS = AP;
AN = AP;

%
%%
%====================================
% Coefficients for interior cells
%====================================

for k = 1:Ny  
    for j = 1:Nx
        m = k + (j-1)*Ny;
        AW(m) = -0.5*dy*ux(j)-kappa*dy/dx;
        AE(m) = 0.5*dy*ux(j+1) -kappa*dy/dx;
        AS(m) = -0.5*dx*uy(k)-kappa*dx/dy;
        AN(m) = 0.5*dx*uy(k+1) -kappa*dx/dy;
        AP(m) = 2*kappa*(dx/dy+dy/dx);
    end
end
%%
% assign coefficients
% interior cells
for k = 2:Ny-1  
    for j = 2:Nx-1
        m = k + (j-1)*Ny;      
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1 ) = AS(m);
        A(m,m+1 ) = AN(m);
%         B(m)      = 0;
    end
end

%%
% fill in the coefficients for the boundary conditions
% south boundary
k = 1;
for j = 2:Nx-1
        m = k + (j-1)*Ny;  
        AP(m)     = kappa*(2*dy/dx+dx/dy);
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m+Ny) = AE(m);
        A(m,m+1)  = AN(m);
%         B(m)      = 0;
end
%
%%
% west boundary
j = 1;
for k = 2:Ny-1
        m = k + (j-1)*Ny;
        AP(m)     = kappa*(3*dy/dx+2*dx/dy);
        A(m,m)    = AP(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1)  = AS(m);
        A(m,m+1)  = AN(m);
        B(m)      = 2*kappa*dy/dx*Tw(k);
end
%%
% north boundary
k = Ny;
for j = 2:Nx-1
        m = k + (j-1)*Ny;
        AP(m)     = -0.5*dx*uy(k+1)+kappa*(2*dy/dx+3*dx/dy);
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m+Ny) = AE(m);
        A(m,m-1)  = AS(m);
        B(m)      = 0;
end
%%
% east boundary
j = Nx;
for k = 2:Ny-1
        m = k + (j-1)*Ny;
        AP(m)     = 0.5*dy*ux(j+1) + kappa*(dy/dx+2*dx/dy);
        A(m,m)    = AP(m);
        A(m,m-Ny) = AW(m);
        A(m,m-1)  = AS(m);
        A(m,m+1)  = AN(m);
        B(m)      = 0;
end
%%
%south-west corner
j = 1;
k = 1;
m = k + (j-1)*Ny;
AP(m)     = kappa*(3*dy/dx+dx/dy);
A(m,m)    = AP(m);
A(m,m+Ny) = AE(m);
A(m,m+1)  = AN(m);
B(m)      = 2*kappa*dy/dx*Tw(k);

%%
%south-east corner
j = Nx;
k = 1;
m = k + (j-1)*Ny;
AP(m)     = kappa*(dy/dx+dx/dy) + 0.5*dy*ux(j+1) ;
A(m,m)    = AP(m);
A(m,m-Ny) = AW(m);
A(m,m+1 ) = AN(m);
B(m)      = 0;

%%
%north-east corner
j = Nx;
k = Ny;
m = k + (j-1)*Ny;
AP(m)     = kappa*(dy/dx+3*dx/dy) + 0.5*(dy*ux(j+1) -dx*uy(k+1));
A(m,m)    = AP(m);
A(m,m-Ny) = AW(m);
A(m,m-1 ) = AS(m);
B(m)      = 0;

%%
%north-west corner
j = 1;
k = Ny;
m = k + (j-1)*Ny;
AP(m)     = 3*kappa*(dy/dx+dx/dy) -0.5*dx*uy(k+1);
A(m,m)    = AP(m);
A(m,m+Ny) = AE(m);
A(m,m-1 ) = AS(m);
B(m)      = 2*kappa*dy/dx*Tw(k);

%%
% 
% now me make A a sparse matrix
A = sparse(A);
%

%% we finally solve the linear system using backslash operator
T = A\B;
%

%% we reshape
TM = reshape(T,Ny,Nx);

%% and plot
[c, h] = contourm(xcenter,ycenter,TM);
clegendm(c,h,4);


