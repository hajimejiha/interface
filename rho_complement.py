import numpy as np
import math

start_step = input("Start Step:")
dump_step = input("Dump Timestep:")
file_num = input("Number of Dumpfile:")

grid_x = 250
grid_y = 250
grid_z = 250

div = 3
length_threshold = 5.0

interface_z = np.zeros((3, grid_x*grid_y))

def length(a, b, c, d, e, f):
	length = math.sqrt((a-b)**2+(c-d)**2+(e-f)**2)
	return length
	
for num in range(int(file_num)):
	step = int(start_step) + int(num) * int(dump_step)
	lineCount = 0
	
	with open("add_z_" + str(step) + ".xyz", mode="w") as f0:
		f0.write("type x y z")
	
	print("Now Step:"+str(step))

	for line in open("interface_Z_" + str(step) + ".xyz",mode="r"):
		lineCount += 1
		data = line.split()
		
		if lineCount < 3:
			point_num = data[0]
			continue

		interface_z[0, lineCount-3] = data[1]
		interface_z[1, lineCount-3] = data[2]
		interface_z[2, lineCount-3] = data[3]
	
	
	for i in range(grid_x-1):
		for j in range(grid_y):
			point_length = length(interface_z[0, i*grid_y+j], \
					  interface_z[0, (i+1)*grid_y+j], \
					  interface_z[1, i*grid_y+j], \
					  interface_z[1, (i+1)*grid_y+j], \
					  interface_z[2, i*grid_y+j], \
					  interface_z[2, (i+1)*grid_y+j])

			if point_length > length_threshold:
				for k in range(div-1):
					x = (interface_z[0, i*grid_y+j] + interface_z[0, (i+1)*grid_y+j]) / div * (k+1)
					y = (interface_z[1, i*grid_y+j] + interface_z[1, (i+1)*grid_y+j]) / div * (k+1)
					z = (interface_z[2, i*grid_y+j] + interface_z[2, (i+1)*grid_y+j]) / div * (k+1)
					with open("add_z_" + str(step) + ".xyz", mode="a") as f0:
						f0.write("\n"+"A "+str(x)+" "+str(y)+" "+str(z))
	
	for i in range(grid_x):
		for j in range(grid_y-1):
			point_length2 = length(interface_z[0, i*grid_y+j], \
					  interface_z[0, i*grid_y+(j+1)], \
					  interface_z[1, i*grid_y+j], \
					  interface_z[1, i*grid_y+(j+1)], \
					  interface_z[2, i*grid_y+j], \
					  interface_z[2, i*grid_y+(j+1)])

			if point_length2 > length_threshold:
				for k in range(div-1):
					x = (interface_z[0, i*grid_y+j] + interface_z[0, i*grid_y+(j+1)]) / div * (k+1)
					y = (interface_z[1, i*grid_y+j] + interface_z[1, i*grid_y+(j+1)]) / div * (k+1)
					z = (interface_z[2, i*grid_y+j] + interface_z[2, i*grid_y+(j+1)]) / div * (k+1)
					with open("add_z_" + str(step) + ".xyz", mode="a") as f0:
						f0.write("\n"+"A "+str(x)+" "+str(y)+" "+str(z))
