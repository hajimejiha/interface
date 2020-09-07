import numpy as np
import math
import sys

args = sys.argv

start_step = int(args[1])
dump_step = int(args[2])
file_num = int(args[3])
rho_center = float(args[4])

grid_x = 260
grid_y = 265
grid_z = 240

x_center = (-2.6094400009644420e-01 + 2.9680005700009639e+02)/2.0
y_center = (-1.2665453803565583e+00 + 3.0217111738035612e+02)/2.0

dx = 1.0

x0 = x_center - (grid_x/2)*dx
y0 = y_center - (grid_y/2)*dx
z0 = 20.0

interface_z = np.zeros((3, grid_x*grid_y))



def length(a, b, c, d, e, f):
	length = math.sqrt((a-b)**2+(c-d)**2+(e-f)**2)
	return length



def interface(step, rho_center):
    itr = 0
    lineCount = 0
    rho_z = 0.0
    rho_z_old = 0.0
    
    for line in open("rho_D." + str(step) +".txt","r"):
        lineCount += 1
        if lineCount < 2:
            continue

        if (lineCount - 2) % grid_z == 0:
            check = False
            
            itr += 1
            rho_z = 0.0
            rho_z_old = 0.0

        if check == False:
            data = line.split()
            rho_z_old = rho_z
            rho_z = float(data[3])

            if (lineCount - 2) % grid_z == 1:
                rho_z_old = rho_z

            if rho_z < rho_center < rho_z_old:
                check = True
                
                interface_z[0, itr] = float(data[0])
                interface_z[1, itr] = float(data[1])
                interface_z[2, itr] = float(data[2]) + dx *  (rho_center - rho_z) / (rho_z_old - rho_z)


    with open("interface_z." + str(step) + ".xyz",mode="w") as f:
        f.write(str(grid_x*grid_y) + "\nInterface")
        for i in range(grid_x*grid_y):
            f.write("\n"+"X "+str(interface_z[0, i])+" "+str(interface_z[1, i])+" "+str(interface_z[2, i]))



def rho_comp(step):
    count = 0
    lineCount = 0
    div = 3
    length_threshold = 5.0
        
    for i in range(grid_x-1):
        for j in range(grid_y):
            point_length = (length(interface_z[0, i*grid_y+j],
                                   interface_z[0, (i+1)*grid_y+j],
                                   interface_z[1, i*grid_y+j],
                                   interface_z[1, (i+1)*grid_y+j],
                                   interface_z[2, i*grid_y+j],
                                   interface_z[2, (i+1)*grid_y+j]))

            if point_length > length_threshold:
                for k in range(div-1):
                    count += 1
                    x = (interface_z[0, i*grid_y+j] + interface_z[0, (i+1)*grid_y+j]) / div * (k+1)
                    y = (interface_z[1, i*grid_y+j] + interface_z[1, (i+1)*grid_y+j]) / div * (k+1)
                    z = (interface_z[2, i*grid_y+j] + interface_z[2, (i+1)*grid_y+j]) / div * (k+1)
                    with open("interface_z_comp." + str(step) + ".xyz", mode="a") as f:
                        f.write("\n"+"A "+str(x)+" "+str(y)+" "+str(z))
    
    for i in range(grid_x):
        for j in range(grid_y-1):
            point_length2 = (length(interface_z[0, i*grid_y+j],
                                    interface_z[0, i*grid_y+(j+1)],
                                    interface_z[1, i*grid_y+j],
                                    interface_z[1, i*grid_y+(j+1)],
                                    interface_z[2, i*grid_y+j],
                                    interface_z[2, i*grid_y+(j+1)]))

            if point_length2 > length_threshold:
                for k in range(div-1):
                    count += 1
                    x = (interface_z[0, i*grid_y+j] + interface_z[0, i*grid_y+(j+1)]) / div * (k+1)
                    y = (interface_z[1, i*grid_y+j] + interface_z[1, i*grid_y+(j+1)]) / div * (k+1)
                    z = (interface_z[2, i*grid_y+j] + interface_z[2, i*grid_y+(j+1)]) / div * (k+1)
                    with open("interface_z_comp." + str(step) + ".xyz", mode="a") as f:
                        f.write("\n"+"A "+str(x)+" "+str(y)+" "+str(z))
                        


if __name__ == '__main__':

    for num in range(file_num):
        step = start_step + num * dump_step
        print("Now Step:"+str(step))
        
        print("Interface")
        interface(step, rho_center)
        
        print("Complement Point\n")
        rho_comp(step)
