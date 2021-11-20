from math import *
from numpy import array
from numpy import linalg 
from numpy import arange 
from matplotlib.pyplot import * 
import matplotlib.pyplot as plt
from numpy import *
import numpy as np
from numpy.typing import _32Bit
from fd_rk45 import *
from lqr2 import *
import icon as i

global K
global cm_b_d_alpha
global az_turbulence
global ax_wind
global az_wind

global c_lin
global c_mass
global c_vel
global c_pow
global n  # dimenstion of x
global m  # dimension of u

n=6
m=2

problem = 2
problem_txt = 'GLIDSLOPE FOLLWOING'

def declare_matrix( m, n ):
  return array( [ [ 0.0 for j in range(0,n)] for i in range(0,m)] )

def declare_vector( n  ):
  return array( [ 0.0 for i in range(0,n) ] )

def RHS( x , t , u_control):
  
  deg2rad = pi/180.0
  g = 9.81

  S = 19.185 * c_lin**2  
  c =  1.607 * c_lin 
  b = 11.85 * c_lin 
  Sh = 3.59  * c_lin**2 
  ch = 0.871 * c_lin 
  Lh = 4.856 * c_lin 
  
  m  = 2070 * c_mass 
  Iy = 4500 * c_mass*(c_lin**2) 

  vx = x[0] + ax_wind
  vz = x[1] + az_wind

  lambd = b*b/S
  kappa = (Sh*Lh)/(S*c)

  alpha = atan( vz/vx )
  V = sqrt( vz**2 + vx**2 )

  a = 4.92#4.513
  alpha_Cz_0 = -3.3*deg2rad #1.8
  CL = a * ( alpha - alpha_Cz_0 )

  e = 0.8 
  CD_0 = 0.0419 #0.023
  CD = CD_0 + CL**2/(pi*lambd*e)

  ro_0 = 1.225
  ro = ro_0 * ( 1.0 - abs(x[4])/44300.0 )**4.256

  Q_dyn = 0.5*ro*(V**2)

  L = Q_dyn * S * CL
  D = Q_dyn * S * CD
  G = m * g

  Power_0_HP = 2*220 * c_pow 
  Power = 735.5 * Power_0_HP * ro/ro_0
  Thrust = Power/V * u_control[1]

  cm_ac = -0.095#-0.05 
  cm_w = -0.110#-0.065 
  cm_b_d_alpha = 0.1982#0.547
  cm_wb = cm_w + cm_b_d_alpha*( alpha - alpha_Cz_0 ) #Cmbu

  alpha_h_0 = -2.0*deg2rad 
  d_eps_d_alfa = 0.4815#0.365 
  a_1 = 3.4#4.06 
  a_2 = 2.513#2.96 
  alpha_h = alpha_h_0 + ( 1.0 - d_eps_d_alfa )*( alpha - alpha_Cz_0 )
  cm_h = a_1*alpha_h + a_2*u_control[0]

  cm_q = -20.14#-7.62 
  om_y_v = x[2] / V #predkosc kątowa pochylania
  
  dx_dt = declare_vector (n)

  dx_dt[0] = ( -D*cos(alpha) + L*sin(alpha) - G*sin(x[5]) + Thrust ) / m  -  x[2]*vz
  #dx_dt(2) = ( -D*sin(alpha) - L*cos(alpha) + G*cos(x(6))          ) / m  +  x(3)*vx
  dx_dt[1] = ( -D*sin(alpha) - L*cos(alpha) + G*cos(x[5])          ) / m  +  x[2]*vx  + az_turbulence
  dx_dt[2] = ( Q_dyn*S*c*( cm_ac + cm_wb - kappa*cm_h + cm_q*om_y_v ) ) / Iy
  dx_dt[3] =  cos(x[5])*vx + sin(x[5])*vz
  dx_dt[4] = -sin(x[5])*vx + cos(x[5])*vz
  dx_dt[5] = x[2]

  return dx_dt

def Jacob_AB( RHS, y, t , u_control , n , m ):
    
  A = declare_matrix( n , n )
  B = declare_matrix( n , m )

  dy = 1.0e-6
  f0 = RHS( y , t , u_control )
  
  for i in range( 0 , n ) :
    yp = array( y )
    yp[i] += dy; 
    f = RHS( yp , t , u_control )
    for j in range( 0 , n ) :
      A[j,i]= ( f[j] - f0[j] ) / dy; 
    
  for i in range( 0 , m ) :
    up = array([u_control[0], u_control[1]])
    up[i] += dy
    f = RHS( y , t , up )
    for j in range( 0 , n ):
      B[j,i]= ( f[j] - f0[j] ) / dy; 
   
  return A , B

plt.figure( figsize = ( 15 , 7 ) )
plt.ion()
ax = plt.gca()
plt.axis( [ -10 , 140 , 0 , 200 ] )

liney, = ax.plot( [] , [] , 'r' )
lineu, = ax.plot( [] , [] , 'g' )
lineg, = ax.plot( [] , [] , 'b' )

pntg, = ax.plot( [] , [] , 'ob' ,marker = i.plane_marker, markersize=40 )

handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels)

c_lin = 17.44/19.185
c_mass = 3018.0/2070.0
c_vel = 255.0/250.0
c_pow = 989.0/(2*220.0*0.7355)

z0 = 120.0
c_turb = 300.0
c_turb = 0.0
X_turb_1 = 500.0
X_turb_2 = 700.0

az_turbulence = 0.0
ax_wind = 0.0
az_wind = 0.0

rad2deg = 180/pi
deg2rad = pi/180

ils_slope = atan(-3.0*deg2rad)

Vx = c_vel * 250.0/3.6
Vz = 0.0
 
x = declare_vector( n ) 
x = [Vx, Vz, 0.0, 0.0, -z0, 0.0]

u_control = declare_vector( m )

tp  = np.zeros(3000)
yp = np.zeros( (len(x),3000) )
up  = np.zeros( (len(u_control),3000) )
gp = np.zeros(3000)
zp = np.zeros(3000)

dt = 0.01
t = 0.0

R = np.eye(m, m)
R[0,0] = 0.435
R[1,1] = 0.004

Q = np.eye(n,n)
Q[0,0] = 0.865
Q[1,1] = 1.555
Q[2,2] = 0.235
Q[3,3] = 0.000
Q[4,4] = 0.245
Q[5,5] = 0.000

i = 0

for i in range (1, 3000):

  X = x[3]
  Z0 = 100.0
  z_ref = Z0 + ils_slope*X
  z_terr = 0.0
  Vz = Vx * ils_slope

  tp[i] = t
  yp[:, i] = x
  up[:, i] = u_control
  gp[i] = z_terr
  zp[i] = z_ref

  e = np.zeros( n )

  e[0] = x[0] - ( cos(x[5])*Vx - sin(x[5])*Vz )
  e[1] = x[1] + ( sin(x[5])*Vx + cos(x[5])*Vz )
  #e(1) = ( cos(x(6))*x(1) + sin(x(6))*x(2)) - Vx
  #e(2) = ( sin(x(6))*x(1) - cos(x(6))*x(2)) - Vz
  e[2] = x[2] - 0.0
  e[3] = 0.0
  e[4] = x[4] - (-z_ref)
  e[5] = 0.0

  [A, B] = Jacob_AB(RHS, x, t, u_control, n, m)
  K = lqr2(A, B, Q, R)[0]

  u_control = matmul(-K, e).A1
  
  u_control[0] = max( -15.0*deg2rad , min( u_control[0] , 15.0*deg2rad ) )
  u_control[1] = max(  0.1 , min( u_control[1] , 1.0 ))

  az_turbulence = 0.0
  ax_wind = 0.0
  az_wind = 0.0
  
  if X > X_turb_1 and X < X_turb_2:
    az_turbulence = c_turb*( 1.0 - 2.0*np.random.rand())
    ax_wind = 0.0
    az_wind = 0.0;  #15.5 + 3.0*( 1.0 - 2.0*rand() )

  
  liney.set_data( yp[3, :i] , zp[:i] )
  lineu.set_data( yp[3, :i] , gp[:i] )
  lineg.set_data( yp[3, :i] , -yp[4, :i] )
  
  pntg.set_data( x[3] , -x[4] )
  
  x = fd_rk45( RHS , x , t , dt, u_control )

  i += 1
  
  if i % 10 == 0: 

    xl = max( 0 , X - 500 )
    xp = X + 50
    zl = -10
    zu = z0 + 40

    v_z = sin(x[5])*x[0] - cos(x[5])*x[1]
    V = sqrt( x[0]**2 + x[1]**2 ) * 3.6
    e_v = Vx - sqrt( ( cos(x[5])*x[0] + sin(x[5])*x[1] )**2 )
    teta = x[5]*rad2deg
    alt = -x[4]
    thr = u_control[1]
    elev = u_control[0]*rad2deg
    gamma_ = atan( ( sin(x[5])*x[0] - cos(x[5])*x[1])  / ( cos(x[5])*x[0] + sin(x[5])*x[1]) )*rad2deg

    txt = '{}:\n  t = {:.9f}   V = {:.9f} km/h   gamma = {:.9f}   v_z = {:.9f} alt = {:.9f}\n teta={:.9f} thr = {:.9f} elev = {:.9f} e_v = {:.9f} e(q) = {:.9f} e(z) = {:.9f}'.format(problem_txt, t , V , gamma_ , v_z, alt, teta, thr, elev, e_v, e[2], e[4])
    
    ax.relim()
    ax.autoscale_view( False, True, False )  
    ax.grid( True )
    #legend([ liney , lineu , lineg, linez ], ["z_ref", "z_terr" , "z", ""] , loc = 2 )
    plt.subplots_adjust( left = 0.07 , right = 0.95 , bottom = 0.1 , top = 0.87 )
    plt.xlabel( '' , size = 15 )
    plt.ylabel( '' , size = 15 )
    plt.axis( [ 0 , xp , zl , zu ] )
    plt.title( txt )
    plt.draw()
    plt.pause(0.01)

  if x[4] > 0:
    plt.show(block=True)
    break

  t += dt
  
