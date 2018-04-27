#Because we don't have a display, we need to put this lines else we'll get errors.

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
plt.ioff() #http://matplotlib.org/faq/usage_faq.html (interactive mode)

import numpy as np

def readfile(file):
  xtemp = []
  ytemp = []
  vtemp = []
  psitemp = []
  with  open(file) as f:
    for line in f:
      if len(line)>1:    
        data = line.split()
        xtemp.append(float(data[0]))
        ytemp.append(float(data[1]))
        psitemp.append(float(data[2]))
        vtemp.append(float(data[3]))
        
    return xtemp,ytemp,psitemp,vtemp

files = ["Lf1.txt", "Lf2.txt","Lf5.txt"]
# files = ["steer5.txt", "steer6.txt","steer7.txt"]
# style = ['bo', 'r--', 'g^']
x = []
y = []
v = []
psi = []

#Create lists to plot data
for file in files:
    xi,yi,psii,vi = readfile("../logs/"+file)
    x.append(xi)
    y.append(yi)
    psi.append(psii)
    v.append(vi)


#Plot
plt.plot(x[0], y[0], label='Lf=1', marker = '*' )
plt.plot(x[1], y[1], label='Lf=2', marker = '*' )
plt.plot(x[2], y[2], label='Lf=5', marker = '*' )

#Plot configuration
minx = min(x[0]+x[1]+x[2])
maxx = max(x[0]+x[1]+x[2])
miny = min(y[0]+y[1]+y[2])
maxy = max(y[0]+y[1]+y[2])
extra_pad = 5

plt.xticks(np.arange(minx, maxx, step=10))
plt.xlim((minx-extra_pad, maxx+extra_pad))   # set the ylim to ymin, ymax
plt.yticks(np.arange(miny, maxy, step=10))
plt.ylim((miny-extra_pad, maxy+extra_pad))   # set the ylim to ymin, ymax
plt.xlabel('x')
plt.ylabel('y')
plt.title("Different Lf values")
plt.legend()

#Save plot to file
plt.savefig('../plots/Lf.png') 
plt.show()
 