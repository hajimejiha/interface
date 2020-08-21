import numpy as np
import math
import sys

args = sys.argv

start_step = int(args[1])
dump_step = int(args[2])
file_num = int(args[3])

grid_x = 260
grid_y = 265
grid_z = 240

x0 = x_center - (grid_x/2)*dx
y0 = y_center - (grid_y/2)*dx
z0 = 20.0

dx = 1.0

interface_z = np.zeros((3, grid_x*grid_y))



def length(a, b, c, d, e, f):
	length = math.sqrt((a-b)**2+(c-d)**2+(e-f)**2)
	return length

    

def histgram():
    bin_num = 1000
    hist_ave = np.zeros(bin_num+1)
    
    for num in range(int(file_num)):
        step = int(start_step) + int(num) * int(dump_step)
        lineCount = 0
        
        rho = np.zeros(grid_x*grid_y*grid_z)
        


        for line in open("rho_D." + str(step) +".txt","r"):
            lineCount += 1
            if lineCount < 3:
                continue		
            data = line.split()

            rho[lineCount - 3] = data[3]
            
        if num == 0:
            min=np.min(rho)
            max=np.max(rho)		

        hist, bins = np.histogram(rho, bins=bin_num, range=(min, max))

        for i in range(bin_num):
            hist_ave[i] += hist[i]/float(file_num)

        with open("histgram_"+str(step)+"_rho_D.txt",mode="w") as f:
            for i in range(bin_num):
                f.write("\n"+str(hist[i])+" "+str(bins[i]))
            f.write("\n0 "+str(bins[bin_num]))

    hist_ave2 = np.delete(hist_ave, 0)
    index_liq = np.where(hist_ave2 == np.max(hist_ave2))

    rho_center = (bins[int(index_liq[0] // 2)] + bins[int(-(-index_liq[0] // 2))]) / 2.0	

    with open("histgram_average_rho_D.txt",mode="w") as f:
        f.write("rho_center:"+str(rho_center))
        for i in range(bin_num+1):
            f.write("\n"+str(hist_ave[i])+" "+str(bins[i]))



def interface(step):
    itr = 0
    lineCount = 0
    rho_z = 0.0
    rho_z_old = 0.0
    chack = 0
    
    for line in open("rho_D." + str(step) +".txt","r"):
        lineCount += 1
        if lineCount < 3:
            continue

        if (lineCount - 2) % grid_z == 0:
            itr += 1
            rho_z = 0.0
            rho_z_old = 0.0
            chack = 0

        if chack == 0:
            data = line.split()

            rho_z_old = rho_z
            rho_z = float(data[3])

            if (lineCount - 2) % grid_z == 1:
                rho_z_old = rho_z

            if rho_z < rho_center < rho_z_old:
                interface_z[0, itr] = float(data[0])
                interface_z[1, itr] = float(data[1])
                interface_z[2, itr] = float(data[2]) + dx *  (rho_center - rho_z) / (rho_z_old - rho_z)
                chack = 1

    with open("interface_Z." + str(step) + ".xyz",mode="w") as f:
        f.write(str(grid_x*grid_y) + "\nInterface")
        for i in range(grid_x*grid_y):
            f.write("\n"+"X "+str(interface_z[0, i])+" "+str(interface_z[1, i])+" "+str(interface_z[2, i]))



def rho_comp(step):
    div = 3
    length_threshold = 5.0
        
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
                    with open("interface_z." + str(step) + ".xyz", mode="a") as f:
                        f.write("\n"+"A "+str(x)+" "+str(y)+" "+str(z))
    
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
                    with open("interface_z." + str(step) + ".xyz", mode="a") as f:
                        f.write("\n"+"A "+str(x)+" "+str(y)+" "+str(z))



if __name__ == '__main__':
    print("Histgram")
    histgram()
    
    for num in range(int(file_num)):
        step = int(start_step) + int(num) * int(dump_step)
        print("Now Step:"+str(step))
        
        print("Interface")
        interface(step)
        
        print("Complement Point")
        rho_comp(step)
