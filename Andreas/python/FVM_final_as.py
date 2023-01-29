import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse as sp
import scipy.sparse.linalg as la

def T_BC_gradient(x):
    return Tc + (Th-Tc)/Lx *x

def Laplacian_T(Nx,Ny,dx,dy):
    # 1. set coefficents for A
    A_E = 1/(dx**2)
    A_W = 1/(dx**2)
    A_S = 1/(dy**2)
    A_N = 1/(dy**2)
    A_P = -(2/(dx**2) + 2/(dy**2))


    
    #Set up matrix Lt
    zeros = np.linspace(Ny-1,Nx*Ny-1,Nx)[:-1].astype("int")
    A_N_vec = A_N * np.ones(Ny*Nx-1)
    A_N_vec[zeros] = 0

    A_S_vec = A_S * np.ones(Ny*Nx-1)
    A_S_vec[zeros] = 0

    Lt = sp.diags([A_P,A_E,A_W,A_N_vec,A_S_vec],offsets = [0,Ny,-Ny,1,-1],shape=[Nx*Ny,Nx*Ny],format="csc")
    #plt.spy(Lt)
    
    #Applying BC to A_P
    A_P_side_indecies = list(np.arange(0,Ny,1)) + list(np.arange(Nx*Ny-Ny,Nx*Ny,1))
    A_P_topbottom_indecies = list(np.linspace(0,Nx*Ny-Ny,Nx).astype("int")) + list(np.linspace(Ny-1,Nx*Ny-1,Nx).astype("int"))
    
    Lt[A_P_side_indecies,A_P_side_indecies] = Lt[A_P_side_indecies,A_P_side_indecies] - 1/(dx**2)
    Lt[A_P_topbottom_indecies,A_P_topbottom_indecies] = Lt[A_P_topbottom_indecies,A_P_topbottom_indecies] -1 /(dy**2) 
    
    #Setting up bt
    bt = np.zeros(Nx*Ny)

    left_BC = - 2/(dx**2)*Tc
    right_BC = - 2/(dx**2)*Th

    bt[0:Ny] = bt[0:Ny] + left_BC
    bt[-Ny:] = bt[-Ny:] + right_BC

    x_p = np.linspace(dx/2,Lx -dx/2,Nx )

    T_BC_grad = T_BC_gradient(x_p)

    top_BC = -2/(dy**2) * T_BC_grad  # top_BC in bottom_BC are identical!
    bottom_BC = -2/(dy**2) * T_BC_grad

    bt[np.linspace(Ny-1,Nx*Ny-1,Nx).astype("int")] = bt[np.linspace(Ny-1,Nx*Ny-1,Nx).astype("int")] + top_BC
    bt[np.linspace(0,Nx*Ny-Ny,Nx).astype("int")] = bt[np.linspace(0,Nx*Ny-Ny,Nx).astype("int")] + bottom_BC
    
    return Lt,bt

#4.6
def tempconv(T_n,u_n,v_n,dx,dp):
    Mat_u_w = sp.diags([1],offsets = [0],shape=[Nx*Ny,(Nx+1)*Ny],format="csc")
    Mat_u_e = sp.diags([1],offsets = [Ny],shape=[Nx*Ny,(Nx+1)*Ny],format="csc")

    Mat_v_s = sp.kron(sp.eye(Nx),sp.diags([1],offsets = [0],shape=[Ny,Ny+1],format="csc"))
    Mat_v_nn = sp.kron(sp.eye(Nx),sp.diags([1],offsets = [1],shape=[Ny,Ny+1],format="csc"))
    
    T_1 = T_n
    T_2 = np.roll(T_n,-Ny)
    T_2[-Ny:] = 2*Th-T_n[-Ny:] #ghost points
    T_e = np.where(Mat_u_e*u_n>=0,T_1,T_2) 

    T_3 = np.roll(T_n,Ny)
    T_3[:Ny] = 2*Tc-T_n[:Ny] #ghost points
    T_4 = T_n 
    T_w = np.where(Mat_u_w*u_n>=0,T_3,T_4) 

    x_p = np.linspace(dx/2,Lx -dx/2,Nx )

    T_BC_grad = T_BC_gradient(x_p)

    T_5 = np.roll(T_n,1)
    T_5[0::Ny] = 2*T_BC_grad -T_n[0::Ny] #ghost points
    T_6 = T_n
    T_s = np.where(Mat_v_s*v_n>=0,T_5,T_6) 

    T_7 = T_n
    T_8 = np.roll(T_n,-1)
    T_8[Ny-1::Ny]= 2*T_BC_grad-T_n[Ny-1::Ny] #ghost points
    T_nn = np.where(Mat_v_nn*v_n>=0,T_7,T_8)  #nn becuase we dont want to confuse nord n with n as n as a time step
    
    return Mat_u_e*u_n * T_e * 1/dx -Mat_u_w*u_n * T_w *1/dx + Mat_v_nn *v_n* T_nn *1/dy - Mat_v_s* v_n * T_s *1/dy


#Functions for u
def Laplacian_u(Nx,Ny,dx,dy):
    #Laplace u
    # u size is ny * (nx-1)

    # 1. set coefficents for A
    A_E = 1/(dx**2)
    A_W = 1/(dx**2)
    A_S = 1/(dy**2)
    A_N = 1/(dy**2)
    A_P = -(2/(dx**2) + 2/(dy**2))



    #Set up matrix Lt
    zeros = np.linspace(Ny-1,(Nx-1)*Ny-1,Nx-1)[:-1].astype("int")
    A_N_vec = A_N * np.ones(Ny*(Nx-1)-1)
    A_N_vec[zeros] = 0

    A_S_vec = A_S * np.ones(Ny*(Nx-1)-1)
    A_S_vec[zeros] = 0

    Lu = sp.diags([A_P,A_E,A_W,A_N_vec,A_S_vec],offsets = [0,Ny,-Ny,1,-1],shape=[(Nx-1)*Ny,(Nx-1)*Ny],format="csc")
    #plt.spy(Lt)

    #Applying BC to A_P
    A_P_topbottom_indecies = list(np.linspace(0,(Nx-1)*Ny-Ny,Nx-1).astype("int")) + list(np.linspace(Ny-1,(Nx-1)*Ny-1,Nx-1).astype("int"))

    Lu[A_P_topbottom_indecies,A_P_topbottom_indecies] = Lu[A_P_topbottom_indecies,A_P_topbottom_indecies] -1 /(dy**2) 

    #Setting up bt <- it is zero
    bu = np.zeros((Nx-1)*Ny)

    left_BC = - 1/(dx**2)*0
    right_BC = - 1/(dx**2)*0

    bu[0:Ny] = bu[0:Ny] + left_BC
    bu[-Ny:] = bu[-Ny:] + right_BC



    top_BC = -2/(dy**2) * 0  
    bottom_BC = -2/(dy**2) * 0

    bu[np.linspace(Ny-1,(Nx-1)*Ny-1,Nx-1).astype("int")] = bu[np.linspace(Ny-1,(Nx-1)*Ny-1,Nx-1).astype("int")] + top_BC
    bu[np.linspace(0,(Nx-1)*Ny-Ny,Nx-1).astype("int")] = bu[np.linspace(0,(Nx-1)*Ny-Ny,Nx-1).astype("int")] + bottom_BC

    return Lu, bu

def u_conv(u_n,v_n,dx,dy):

    # u convection 
    u_n_lim = u_n[Ny:-Ny] # this part takes out the boundary u components from u_n , the first and last Ny

    # u_e^2
    u_P = u_n_lim
    u_E = np.roll(u_n_lim,-Ny)
    u_E[-Ny:]= 0 # right BC

    u_e_2 = np.where((u_E + u_P)/2 >=0, (u_E+u_P)/2 *u_P,(u_E+u_P)/2 * u_E )

    # u_w^2
    u_P = u_n_lim
    u_W = np.roll(u_n_lim,Ny)
    u_W[:Ny] = 0 # left BC

    u_w_2 = np.where((u_P + u_W)/2 >=0, (u_P+u_W)/2 *u_W,(u_P+u_W)/2 * u_P )

    # u_n *v_n
    u_P = u_n_lim
    u_N = np.roll(u_n_lim,-1)
    u_N[Ny-1::Ny] = - u_n_lim[Ny-1::Ny] # upper BC

    v_nw = np.delete(v_n[:-(Ny+1)], np.arange(0,len(v_n)-(Ny+1),1)[0::(Ny+1)])
    v_ne = np.delete(v_n[(Ny+1):], np.arange(0,len(v_n)-(Ny+1),1)[0::(Ny+1)])

    u_n_v_n = np.where((v_ne+v_nw)/2 >= 0,(v_ne+v_nw)/2 * u_P,(v_ne+v_nw)/2 * u_N)

    # u_s * v_s
    u_P = u_n_lim
    u_S = np.roll(u_n_lim,1)
    u_S[0::Ny] = - u_n_lim[0::Ny] # bottom BC

    v_sw = np.delete(v_n[:-(Ny+1)], np.arange(0,len(v_n)-(Ny+1),1)[(Ny+1)-1::(Ny+1)])
    v_se = np.delete(v_n[(Ny+1):], np.arange(0,len(v_n)-(Ny+1),1)[(Ny+1)-1::(Ny+1)])

    u_s_v_s = np.where((v_se+v_sw)/2 >= 0,(v_se+v_sw)/2 * u_S,(v_se+v_sw)/2 * u_P)
    return 1/dx * u_e_2 - 1/dx * u_w_2 + 1/dy * u_n_v_n - 1/dy * u_s_v_s

#Functions for v
def Laplacian_v(Nx,Ny,dx,dy):
    #Laplace v
    # v size is (ny-1) *nx
    # 1. set coefficents for A
    A_E = 1/(dx**2)
    A_W = 1/(dx**2)
    A_S = 1/(dy**2)
    A_N = 1/(dy**2)
    A_P = -(2/(dx**2) + 2/(dy**2))



    #Set up matrix Lt
    zeros = np.linspace(Ny-1-1,Nx*(Ny-1)-1,Nx)[:-1].astype("int")
    A_N_vec = A_N * np.ones((Ny-1)*Nx-1)
    A_N_vec[zeros] = 0

    A_S_vec = A_S * np.ones((Ny-1)*Nx-1)
    A_S_vec[zeros] = 0

    Lv = sp.diags([A_P,A_E,A_W,A_N_vec,A_S_vec],offsets = [0,(Ny-1),-(Ny-1),1,-1],shape=[Nx*(Ny-1),Nx*(Ny-1)],format="csc")
    #plt.spy(Lv)

    #Applying BC to A_P
    A_P_side_indecies = list(np.arange(0,Ny-1,1)) + list(np.arange(Nx*(Ny-1)-(Ny-1),Nx*(Ny-1),1))

    Lv[A_P_side_indecies,A_P_side_indecies] = Lv[A_P_side_indecies,A_P_side_indecies] - 1/(dx**2)

    #Setting up bt <- it is zero
    bv = np.zeros(Nx*(Ny-1))

    left_BC = - 2/(dx**2)*0
    right_BC = - 2/(dx**2)*0

    bv[0:(Ny-1)] = bv[0:(Ny-1)] + left_BC
    bv[-(Ny-1):] = bv[-(Ny-1):] + right_BC



    top_BC = -1/(dy**2) * 0  
    bottom_BC = -1/(dy**2) * 0

    bv[np.linspace(Ny-2,(Nx-1)*Ny-1,Nx).astype("int")] = bv[np.linspace(Ny-2,(Nx-1)*Ny-1,Nx).astype("int")] + top_BC
    bv[np.linspace(0,Nx*(Ny-1)-(Ny-1),Nx).astype("int")] = bv[np.linspace(0,Nx*(Ny-1)-(Ny-1),Nx).astype("int")] + bottom_BC

    return Lv, bv

def v_conv(u_n,v_n,dx,dy):
    # v convection 
    v_n_lim = np.delete(v_n, list(np.linspace(Ny,len(v_n)-1,Nx).astype("int"))+list( np.linspace(0,len(v_n)-(Ny+1),Nx).astype("int"))) # this part takes out the boundary v components from v_n

    # v_n^2
    v_P = v_n_lim
    v_N = np.roll(v_n_lim,-1)
    v_N[Ny-2::(Ny-1)] = 0 # upper BC
    v_n_2 = np.where((v_N+v_P)/2 >=0,(v_N+v_P)/2 *v_P, (v_N+v_P)/2 *v_N)

    # v_s^2
    v_P = v_n_lim
    v_S = np.roll(v_n_lim,1)
    v_S[0:-(Ny-2):(Ny-1)] = 0 # upper BC
    v_s_2 = np.where((v_S+v_P)/2 >=0,(v_S+v_P)/2 *v_S, (v_S+v_P)/2 *v_P)

    # v_e * u_e 

    v_P = v_n_lim
    v_E =  np.roll(v_n_lim,-(Ny-1))
    v_E[-(Ny-1):] = 0 #BC

    u_en = np.delete(u_n[Ny:],np.linspace(0,Ny*Nx-Ny,Nx).astype("int"))
    u_es = np.delete(u_n[Ny:],np.linspace(Ny-1,Ny*Nx-1,Nx).astype("int"))

    v_e_u_e = np.where((u_en+u_es)/2 >=0,(u_en+u_es)/2 *v_P, (u_en+u_es)/2 *v_E)

    # v_w * u_w

    v_P = v_n_lim
    v_W =  np.roll(v_n_lim,(Ny-1))
    v_W[:(Ny-1)] = 0 #BC

    u_wn = np.delete(u_n[:-Ny],np.linspace(0,Ny*Nx-Ny,Nx).astype("int"))
    u_ws = np.delete(u_n[:-Ny],np.linspace(Ny-1,Ny*Nx-1,Nx).astype("int"))

    v_w_u_w = np.where((u_wn+u_ws)/2 >=0,(u_wn+u_ws)/2 *v_W, (u_wn+u_ws)/2 *v_P)

    return 1/dy * v_n_2 - 1/dy * v_s_2 + 1/dx * v_e_u_e - 1/dx *v_w_u_w

def Laplacian_p(Nx,Ny,dx,dy):
    # 1. set coefficents for A
    A_E = 1/(dx**2)
    A_W = 1/(dx**2)
    A_S = 1/(dy**2)
    A_N = 1/(dy**2)
    A_P = -(2/(dx**2) + 2/(dy**2))


    
    #Set up matrix Lt
    zeros = np.linspace(Ny-1,Nx*Ny-1,Nx)[:-1].astype("int")
    A_N_vec = A_N * np.ones(Ny*Nx-1)
    A_N_vec[zeros] = 0

    A_S_vec = A_S * np.ones(Ny*Nx-1)
    A_S_vec[zeros] = 0

    Lp = sp.diags([A_P,A_E,A_W,A_N_vec,A_S_vec],offsets = [0,Ny,-Ny,1,-1],shape=[Nx*Ny,Nx*Ny],format="csc")
    #plt.spy(Lt)
    
    #Applying BC to A_P
    A_P_side_indecies = list(np.arange(0,Ny,1)) + list(np.arange(Nx*Ny-Ny,Nx*Ny,1))
    A_P_topbottom_indecies = list(np.linspace(0,Nx*Ny-Ny,Nx).astype("int")) + list(np.linspace(Ny-1,Nx*Ny-1,Nx).astype("int"))
    
    Lp[A_P_side_indecies,A_P_side_indecies] = Lp[A_P_side_indecies,A_P_side_indecies] + 1/(dx**2)
    Lp[A_P_topbottom_indecies,A_P_topbottom_indecies] = Lp[A_P_topbottom_indecies,A_P_topbottom_indecies] +1 /(dy**2) 
    
    #THIS IS TO FIX THE BOTOM MOST RIGHT EDGE TO zero Pa for
    Lp[0,0]=Lp[0,0] -2* 1/(dy**2) # -1 for Ps = 0 => PS= -Pa and -1 for what was added just before by normal BC
 
    return Lp



#Code

#data from 5.1.
Lx = 1
Ly = 1
Nx = 50
Ny = 50

dt = 0.01
t_end = 150
Pr = 0.71
Gr = 3.1*10**6

#read from instructions
Tc = 0 
Th = Lx #if we decided to scale to e.g. 2 instead of 1

dx = Lx/Nx
dy = Ly/Ny

x_p = np.linspace(dx/2,Lx -dx/2,Nx )

T_BC_grad = T_BC_gradient(x_p)
T_IC = np.repeat(T_BC_grad,Ny) #put T_BC_gradient over the whole T grid points

#Inital Conditions
u_n = np.zeros((Nx+1)*Ny)
v_n = np.zeros(Nx*(Ny+1))
T_n = T_IC

#### Contruction of the LHS ####

# T
L_T, b_T = Laplacian_T(Nx,Ny,dx,dy)
I_T = sp.eye(Nx*Ny) #Identity matrix
LHS_T = 1/dt * I_T - 1/(Pr*np.sqrt(Gr)) *L_T

# u* = u_star
L_u_star, b_u_star = Laplacian_u(Nx,Ny,dx,dy)
I_u_star = sp.eye((Nx-1)*Ny)
LHS_u_star = 1/dt * I_u_star - 1/(np.sqrt(Gr))*L_u_star

#v* = v_star
L_v_star, b_v_star = Laplacian_v(Nx,Ny,dx,dy)
I_v_star = sp.eye(Nx*(Ny-1))
LHS_v_star = 1/dt * I_v_star - 1/(np.sqrt(Gr))*L_v_star

#p_np1
LHS_p_np1 = Laplacian_p(Nx,Ny,dx,dy)

#Factorizing 
dLHS_T = sp.linalg.splu(sp.csc_matrix(LHS_T))
dLHS_u_star = sp.linalg.splu(sp.csc_matrix(LHS_u_star))
dLHS_v_star = sp.linalg.splu(sp.csc_matrix(LHS_v_star))
dLHS_p_np1 = sp.linalg.splu(sp.csc_matrix(LHS_p_np1))

#### Stepping through time
print("Running...")
for t in np.arange(0,t_end,dt):
    #RHS temperature
    T_convective_term = tempconv(T_n,u_n,v_n,dx,dy)
    RHS_T = - T_convective_term -  1/(Pr*np.sqrt(Gr))*b_T + 1/dt * T_n
    #solve temperature 
    T_np1 = dLHS_T.solve(RHS_T)

    #RHS u_star
    u_star_convective_term = u_conv(u_n,v_n,dx,dy)
    u_n_lim = u_n[Ny:-Ny]  #removing BC from v so that the size of u_n matches u_star 
    RHS_u_star = - 1/(np.sqrt(Gr)) * b_u_star + 1/dt * u_n_lim - u_star_convective_term  
    #solve u_star
    u_star = dLHS_u_star.solve(RHS_u_star)

    #RHS v_star
    v_star_convective_term = v_conv(u_n,v_n,dx,dy)
    v_n_lim = np.delete(v_n, list(np.linspace(Ny,len(v_n)-1,Nx).astype("int"))+list( np.linspace(0,len(v_n)-(Ny+1),Nx).astype("int"))) #removing BC from v so that the size of v_n matches v_star 
    v_T_term = (np.delete(T_n,np.linspace(0,Nx*Ny-Ny,Nx).astype(int))+ np.delete(T_n,np.linspace(Ny-1,Nx*Ny-1,Nx).astype(int)))/2 #check sign
    RHS_v_star = - 1/(np.sqrt(Gr)) * b_v_star + 1/dt * v_n_lim - v_star_convective_term + v_T_term
    #solve v_star
    v_star = dLHS_v_star.solve(RHS_v_star)

    #RHS pressure
    u_e = u_n[Ny:]
    u_w = u_n[:-Ny]
    v_s = np.delete(v_n, list(np.linspace(Ny,len(v_n)-1,Nx).astype("int")))
    v_nn = np.delete(v_n, list(np.linspace(0,len(v_n)-(Ny+1),Nx).astype("int")))
    div_vec_u_star = 1/dx *u_e - 1/dx *u_w + 1/dy * v_nn -1/dy*v_s
    RHS_p_np1 = 1/dt * div_vec_u_star
    #solve pressure
    p_np1 = dLHS_p_np1.solve(RHS_p_np1)

    # Doing some reshaping to use the diff function of numpy
    p_np1_mat = p_np1.reshape(Nx,Ny)
    p_np1_mat = p_np1_mat.transpose()

    dp_dx = np.diff(p_np1_mat, axis = 1) / dx
    dp_dy = np.diff(p_np1_mat, axis = 0) / dy

    dp_dx = dp_dx.transpose()
    dp_dy = dp_dy.transpose()
    dp_dx = dp_dx.flatten()
    dp_dy = dp_dy.flatten()


    #computing corrections u'=u_line and v'=v_line
    # p_np1_E = p_np1[Ny:]
    # p_np1_W = p_np1[:-Ny]
    # u_line =- dt/dx*(p_np1_E-p_np1_W) #maybe missing a rho?
    u_line =- dt * dp_dx #maybe missing a rho?

    # p_np1_N = np.delete(p_np1, list(np.linspace(0,len(p_np1)-Ny,Nx).astype("int")))
    # p_np1_S = np.delete(p_np1, list(np.linspace(Ny-1,len(p_np1)-1,Nx).astype("int")))
    # v_line = - dt/dy*(p_np1_N-p_np1_S)
    v_line = - dt * dp_dy

    # computing un+1 and vn+1
    u_np1 = u_star + u_line
    v_np1 = v_star + v_line

    #adding back to known BC to u and v
    u_np1 = np.pad(u_np1, (Ny, Ny), 'constant', constant_values=(0, 0)) # scale back u* to the size of u_n by adding left and right BC, since they are 0, we can simply use the pad function to add zeros left and right
    v_np1 = np.pad(np.reshape(v_np1,(Nx,Ny-1)),((0,0),(1,1)), 'constant').flatten() #Adding back top and bottom BC
    
    #preparing for next itteration
    T_n = T_np1
    u_n = u_np1
    v_n = v_np1
    p_n = p_np1

print("Done")

#Displaying things

x_p = np.linspace(dx/2,Lx -dx/2,Nx )
y_p = np.linspace(dy/2,Ly -dy/2,Ny )
yv, xv  = np.meshgrid(y_p,x_p) # creating a mesh

#"D matrix-reshaping vectors and calculating required parameters for visualisation
T_n_vis  = T_n.reshape(Nx,Ny)
p_n_vis  = p_n.reshape(Nx,Ny)

u_e = u_n[Ny:]
u_w = u_n[:-Ny]
u_centerpoint = (u_e + u_w)/2
u_centerpoint = u_centerpoint.reshape(Nx,Ny)

v_nn = np.delete(v_n, np.linspace(0,len(v_n)-(Ny+1),Nx).astype("int"))
v_s = np.delete(v_n, np.linspace(Ny,len(v_n)-1,Nx).astype("int"))
v_centerpoint = (v_nn+v_s)/2
v_centerpoint = v_centerpoint.reshape(Nx,Ny)

#%%

# Plotting
#---------------------------

plt.contourf(xv,yv,T_n_vis,cmap="turbo")
# plt.contour(xv,yv,T_n_vis,cmap="turbo")
plt.colorbar()
# plt.quiver(xv,yv,u_centerpoint,v_centerpoint)
plt.title('Temperature (Dimensionless)')
plt.xlabel('Length in x direction (dim.less)')
plt.ylabel('Length in y direction (dim.less)')
plt.savefig('temperature_fc_python.pdf')
plt.show()


#%%
plt.quiver(xv,yv,u_centerpoint,v_centerpoint)
plt.show()


#%%



#%%

# Plot streamlines




def Laplacian_psi(Nx, Ny, dx, dy):
    # 1. set coefficents for A
    A_E = 1 / (dx ** 2)
    A_W = 1 / (dx ** 2)
    A_S = 1 / (dy ** 2)
    A_N = 1 / (dy ** 2)
    A_P = -(2 / (dx ** 2) + 2 / (dy ** 2))

    # Set up matrix Lt
    zeros = np.linspace(Ny - 1, Nx * Ny - 1, Nx)[:-1].astype("int")
    A_N_vec = A_N * np.ones(Ny * Nx - 1)
    A_N_vec[zeros] = 0

    A_S_vec = A_S * np.ones(Ny * Nx - 1)
    A_S_vec[zeros] = 0

    Lp = sp.diags([A_P, A_E, A_W, A_N_vec, A_S_vec], offsets=[0, Ny, -Ny, 1, -1], shape=[Nx * Ny, Nx * Ny],
                  format="csc")
    # plt.spy(Lt)

    # Applying BC to A_P
    A_P_side_indecies = list(np.arange(0, Ny, 1)) + list(np.arange(Nx * Ny - Ny, Nx * Ny, 1))
    A_P_topbottom_indecies = list(np.linspace(0, Nx * Ny - Ny, Nx).astype("int")) + list(
        np.linspace(Ny - 1, Nx * Ny - 1, Nx).astype("int"))

    Lp[A_P_side_indecies, A_P_side_indecies] = Lp[A_P_side_indecies, A_P_side_indecies] + 1 / (dx ** 2)
    Lp[A_P_topbottom_indecies, A_P_topbottom_indecies] = Lp[A_P_topbottom_indecies, A_P_topbottom_indecies] + 1 / (
                dy ** 2)

    # THIS IS TO FIX THE BOTOM MOST RIGHT EDGE TO zero Pa for
    Lp[0, 0] = Lp[0, 0] - 2 * 1 / (
                dy ** 2)  # -1 for Ps = 0 => PS= -Pa and -1 for what was added just before by normal BC

    return Lp

#%%

x_corners = np.arange(0, Lx, dx)
y_corners = np.arange(0, Ly, dy)

nx = np.size(x_corners) - 1
ny = np.size(y_corners) - 1

dx = x_corners[1] - x_corners[0]
dy = y_corners[1] - y_corners[0]

u_n_mat = u_n[0:Nx * Ny].reshape(Nx, Ny)
u_n_mat = u_n_mat.transpose()
u_n_mat = u_n_mat[:, 1:]

v_n_mat = v_n[0:Nx * Ny].reshape(Nx, Ny)
v_n_mat = v_n_mat.transpose()
v_n_mat = v_n_mat[1:, :]

rhs = np.diff(u_n_mat, axis=0) / dy - np.diff(v_n_mat, axis = 1)

LHS_psi = Laplacian_psi(nx, ny, dx, dy)
dLHS_psi = sp.linalg.splu(sp.csc_matrix(LHS_psi))

psi = dLHS_psi.solve(rhs.flatten())




#%%

plt.contour(xv, yv, psi)
plt.show()

