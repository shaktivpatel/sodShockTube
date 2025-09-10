#sod shock tube numerical solution using MacCormack (2 step Lax-Wendroff Scheme) scheme

#set up computational domain


#set initial conditons of shock tube problem and variables
t = .1
x0 = 0

rho_l = 1
P_l = 1
u_l = 0

rho_r = .125
P_r = .1
u_r = 0

gamma = 1.4
mu = sqrt((gamma - 1)/(gamma + 1))

