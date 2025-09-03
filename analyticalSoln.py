import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

#set initial conditons and duration of sim
t = .1
x0 = 0

rho_l = 1
P_l = 1
u_l = 0

rho_r = .125
P_r = .1
u_r = 0

gamma = 1.4
mu = np.sqrt((gamma - 1)/(gamma + 1))

#calculate speed of sound
c_l = ((gamma * P_l)/rho_l) ** .5
c_r = ((gamma * P_r)/rho_r) ** .5

#calculate valuea 
def func(P):
    return (P - P_r)*((((1 - mu ** 2) ** 2)*((rho_r * (P + mu * mu * P_r)) ** - 1) ) ** .5) - 2 * (np.sqrt(gamma)/(gamma - 1))*(1 - P ** ((gamma - 1)/(2*gamma)))

P_post = fsolve(func, 0)
u_post = 2 * np.sqrt(gamma)/(gamma -1) * (1- P_post ** ((gamma-1)/(2*gamma)))
rho_post = rho_r * ((P_post/P_r) + mu ** 2)/(1 + mu * mu * (P_post/P_r))
u_shock = u_post * ((rho_post/rho_r)/((rho_post/rho_r)-1))
rho_mid = rho_l * (P_post/P_l) ** (1/gamma)

x1 = x0 - c_l * t
x3 = x0 + u_post * t
x4 = x0 + u_shock * t
c_2 = c_l - ((gamma - 1)/2) * u_post
x2 = x0 + (u_post - c_2) * t

n = 1000
x_min = -.5
x_max = .5
x = np.linspace(x_min, x_max, n)

rho_data = np.zeros(n)
P_data = np.zeros(n)
u_data = np.zeros(n)

for i in range(len(x)):
    #region 1
    if x[i] < x1:
        rho_data[i] = rho_l
        P_data[i] = P_l
        u_data[i] = u_l
        
    #region 2
    elif x1 <= x[i] and x[i] <= x2:
        c = mu * mu * ((x0 - x[i])/t) + (1 - mu * mu) * c_l
        rho_data[i] = rho_l * (c/c_l) ** (2/(gamma - 1))
        P_data[i] = P_l * (rho_data[i]/rho_l) ** gamma
        u_data[i] = (1 - mu * mu) * ((-(x0 - x[i])/t) + c_l)
        
    #region 3
    elif x2 <= x[i] and x[i] <= x3:
        rho_data[i] = rho_mid
        P_data[i] = P_post
        u_data[i] = u_post
        
    #region 4
    elif x3 <= x[i] and x[i] <= x4:
        rho_data[i] = rho_post
        P_data[i] = P_post
        u_data[i] = u_post
    
    #region 5    
    elif x4 < x[i]:
        rho_data[i] = rho_r
        P_data[i] = P_r
        u_data[i] = u_r         
                

#plot data vs. x
f1 = plt.figure()
plt.plot(x, rho_data)
plt.xlabel('x')
plt.ylabel('Density')
plt.title('Density vs. Position')
plt.grid()
plt.xlim(-.5,.5)

f2 = plt.figure()
plt.plot(x, P_data)
plt.xlabel('x')
plt.ylabel('Pressure')
plt.title('Pressure vs. Position')
plt.grid()
plt.xlim(-.5,.5)

f3 = plt.figure()
plt.plot(x, u_data)
plt.xlabel('x')
plt.ylabel('Velocity')
plt.title('Velocity vs. Position')
plt.grid()
plt.xlim(-.5,.5)

plt.show()

        

      