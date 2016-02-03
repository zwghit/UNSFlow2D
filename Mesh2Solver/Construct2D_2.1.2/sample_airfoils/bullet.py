import numpy as np
from matplotlib import pyplot as plt
from math import pi, sin, cos

def main():

  # Length and width

  l = 1.0
  w = 2./3.

  # Number of points on curve

  numcrv = 101

  # Allocate x and y

  npoint = numcrv+4
  x = np.zeros((npoint))
  y = np.zeros((npoint))

  # Back/top points

  x[0] = l
  y[0] = 0.
  x[1] = l
  y[1] = w/2.
  
  # Semicircle points

  rad = w/2.
  cenx = rad
  ceny = 0.
  dtheta = pi/float(numcrv-1)
  theta0 = pi/2.
  for i in range(0,numcrv):
    x[i+2] = cenx + rad*cos(theta0 + i*dtheta)
    y[i+2] = ceny + rad*sin(theta0 + i*dtheta)

  # Back/bottom points
  
  x[npoint-2] = l
  y[npoint-2] = -w/2.
  x[npoint-1] = l
  y[npoint-1] = 0.

  # Plot

  ax = plt.subplot(111)
  ax.set_aspect('equal','datalim')
  ax.plot(x,y,'-x')
  plt.show()

  # Write to file

  f = open('bullet.dat','w')
  f.write('Bullet\n')
  for i in range(0,npoint):
    f.write(str(x[i]) + '  ' + str(y[i]) + '\n')
  f.close()

if __name__ == "__main__":
  main()
