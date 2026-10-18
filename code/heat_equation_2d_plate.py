# -*- coding: utf-8 -*-
"""
Created on Thu Feb 15 12:27:05 2024
Nom : BOUTABIA                    Prénom :Yanis                 Classe : 3PB1
@author: ipsa
"""
###################################### Etape 1 #############################
# numpy : Gère les opérations numériques, comme les manipulations de tableaux.
# matplotlib.pyplot : Crée des graphiques statiques et visualise les données.
# matplotlib.animation : Crée des animations pour visualiser l'évolution de la température au fil du temps.
import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import FuncAnimation

print("2D heat equation solver,finite difference method")

###################################### Etape 2 #############################

plate_length = 50
max_iter_time = 100   # 100

alpha = 1     
delta_x = 1   

delta_t = (delta_x ** 2)/(4 * alpha)
gamma = (alpha * delta_t) / (delta_x ** 2)
# solution will be numerically stable whenever delta_T < or = to delta _x^2/4*alpha
print("delta_t=",delta_t,"delta_x=",delta_x,"gamma=",gamma," conditions de stabilité delta_t <= 1/gamma ?")

###################################### Etape 3 #############################
# Commentez 
# Initialize solution: the grid of T(k, i, j)
T = np.empty((max_iter_time, plate_length, plate_length))
###################################### Etape 4    #############################
# Initial condition everywhere inside the grid
T_initial = 400

# Boundary conditions : T_top, T_bottom, T_left, T_right : Temperatures sur les quatres cotés.
T_top = 800     # 500 °C    #   500  #  500  #   800  #   100 # 20
T_left = 250     # 20  °C    #   20  #  500  #   250  #  200 # 250
T_bottom = 40   # 20   °C   #   500 #  500  #   40  #  500   # 40
T_right = 800    # 20  °C    #   20   # 20    #   800 #   200  # 600
print('T_top=', T_top,'°C,','T_left=',T_left,'°C,','T_bottom=',T_bottom ,'°C,','T_right =',T_right,'°C,')
print()

###################################### Partie 5 #############################
# Set the initial condition
T.fill(T_initial)
#print('T.fill=',T.fill(T_initial))

# Set the boundary conditions
T[:, (plate_length-1):, :] = T_top
T[:, :, :1] = T_left
T[:, :1, 1:] = T_bottom
T[:, :, (plate_length-1):] = T_right
###################################### Partie 6 ###################################################################################
# Comment 

#print(T[:, (plate_length-1):, :])
#print(T[:, :, :1])

def calculate(T):
    for k in range(0, max_iter_time-1, 1):
        for i in range(1, plate_length-1, delta_x):
            for j in range(1, plate_length-1, delta_x):
                T[k + 1, i, j] = gamma * (T[k][i+1][j] + T[k][i-1][j] + T[k][i][j+1] + T[k][i][j-1] - 4*T[k][i][j]) + T[k][i][j]
 #               print(k,i,j)
    return T
###################################### Etape 7  ################################################################################
# Comment
def plotheatmap(T_k, k):
    # Clear the current plot figure
    plt.clf()
    plt.title(f"Temperature at t = {k*delta_t:.3f} unit time")
    #plt.grid()
   #plt.title(f"T_top=T_top=', ,'°C,)
  
    plt.xlabel("x")
  
 
    plt.ylabel("y")
  
    # This is to plot T_k (T at time-step k)
    # color bar and value of max and min of temperature
    plt.pcolormesh(T_k, cmap=plt.cm.jet, vmin=0, vmax=800)
    plt.colorbar()
    
    return plt
    
###################################### Etape 8  #####################################################################################
# Comment
# Do the calculation here
T = calculate(T)
###################################### Etape neuf  #####################################################################################
# comment 
def animate(k):
    plotheatmap(T[k], k)

anim = animation.FuncAnimation(plt.figure(), animate, interval=1, frames=max_iter_time, repeat=False)
#

#anim.save("fdm_2d_heat_equation_solution.gif")
#print("Done!")
#print("**********************************************************************")
#print("2D heat equation solver,finite difference method")
#print("**********************************************************************")
#print('T initial =', T_initial,'°C')
#print("**********************************************************************")
#print('T_top=', T_top,'°C,','T_left=',T_left,'°C,','T_bottom=',T_bottom ,'°C,','T_right =',T_right,'°C,')
#print("**********************************************************************")
#print("delta_t=",delta_t,"delta_x=",delta_x,"gamma=",gamma,)

anim.save("fdm_2d_heat_equation_solution.html")
#anim.save("heat_equation_solution.gif")
#plt
#anim.save("heat_equation_solution.gif")

k = 1
plotheatmap(T[k], k)
plt.show()

k = 10
plotheatmap(T[k], k)
plt.show()

k= 99
plotheatmap(T[k], k)
plt.show()


# Look at your folder  "fdm_2d_heat_equation_solution.html"  + "fdm_2d_heat_equation_solution_frames"




