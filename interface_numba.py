import math
import numpy as np
import time
import sys
from numba import njit, prange, set_num_threads

args = sys.argv

start_step = int(args[1])
dump_step = int(args[2])
file_num = int(args[3])
threads_num = int(args[4])

time_start = time.time()

set_num_threads(int(threads_num))

dx = 1.0
K = 10.0

grid_x = 260
grid_y = 265
grid_z = 240

x_center = 100.0
y_center = 100.0

x0 = x_center - (grid_x/2)*dx
y0 = y_center - (grid_y/2)*dx
z0 = 0.0

x1 = x0 - dx * K
x2 = x_center + (grid_x/2)*dx + dx * K
y1 = y0 - dx * K
y2 = y_center + (grid_y/2)*dx + dx * K
z1 = z0 - dx * K
z2 = z0 + grid_z*dx + dx * K

count = 0

@njit('f8[:, :, :](f8, f8, f8)', parallel=True, nogil=True)
def rho(x, y, z):
    C_1 = (2.0*K*dx)**-3.0
    C_2 = math.cos(math.pi / (K * dx)
    
	temp = np.zeros((grid_x, grid_y, grid_z))
	for i in prange(grid_x):
    C_x1 = (i*dx+x0)-x
    C_x2 = C_x1**2.0
    
		for j in range(grid_y):
        C_y1 = (j*dx+y0)-y
        C_y2 = C_y1**2.0
        
			for k in range(grid_z):
            C_z1 = (k*dx+z0)-z
            C_z2 = C_z1**2.0
            
				if math.sqrt(C_x1+C_y1+C_z1) < K * dx:
					temp[i,j,k] = C_1 * (1.0 + C_2 * C_x1) * (1.0 + C_2 * C_y1) * (1.0 + C_2 * C_z1)
	return temp

for num in range(int(file_num)):
	step = int(start_step) + int(num) * int(dump_step)
	lineCount = 0

	rho_D = np.zeros((grid_x, grid_y, grid_z))

	print("\nNow Step:"+str(step))

	for line in open(str(step) +".dmp","r"):
		lineCount += 1
		if lineCount < 10:
			continue
		#print("Line:"+str(lineCount))		
		data = line.split()
		
		X = float(data[3])
		Y = float(data[4])
		Z = float(data[5])

		if int(data[1]) == 1:
			if z2 > Z > z1 and x2 > X > x1 and y2 > Y > y1:
				count += 1
				rho_D += rho(X, Y, Z)

	with open("rho_D_"+str(step)+".txt",mode="w") as f:
		f.write("step:"+str(step)+" atoms:"+str(count)+" dx:"+str(dx)+" K:"+str(K))
		for i in range(grid_x):
			for j in range(grid_y):
				for k in range(grid_z):
					f.write("\n"+str(i*dx+x0)+" "+str(j*dx+y0)+" "+str(k*dx+z0)+" "+str(rho_D[i,j,k]))
	count = 0

	time_end = time.time()
	tim = (time_end - time_start) / 60.0

	est_time = tim / (num + 1.0) * (float(file_num) - num - 1.0)

	with open("Calculation_time.txt",mode="w") as f2:
		f2.write("\n"+str(step)+" "+str(tim))
		
	print("Calculation Time:"+str(tim)+" (min)")
	print("Remaining Time:"+str(est_time)+" (min)")
