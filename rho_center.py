import numpy as np

start_step = int(input("Start Step:"))
dump_step = int(input("Dump Timestep:"))
file_num = int(input("Number of Dumpfile:"))

grid_x = 260
grid_y = 265
grid_z = 240

bin_num = 1000
hist_ave = np.zeros(bin_num+1)

for num in range(file_num):
    step = start_step + num * dump_step
    print("Now Step:"+str(step))
    
    lineCount = 0
    
    rho = np.zeros(grid_x*grid_y*grid_z)
    
    for line in open("rho_D." + str(step) +".txt","r"):
        lineCount += 1
        if lineCount < 2:
            continue		
        data = line.split()

        rho[lineCount - 2] = data[3]
        
    if num == 0:
        min=np.min(rho)
        max=np.max(rho)		

    hist, bins = np.histogram(rho, bins=bin_num, range=(min, max))

    for i in range(bin_num):
        hist_ave[i] += hist[i]/float(file_num)
        
"""
    with open("histgram."+str(step)+".rho_D.txt",mode="w") as f:
        for i in range(bin_num):
            f.write("\n"+str(hist[i])+" "+str(bins[i]))
        f.write("\n0 "+str(bins[bin_num]))
"""

hist_ave2 = np.delete(hist_ave, 0)
index_liq = np.where(hist_ave2 == np.max(hist_ave2))

rho_center = (bins[int(index_liq[0] // 2)] + bins[int(-(-index_liq[0] // 2))]) / 2.0	

with open("histgram_average_rho_D.txt",mode="w") as f:
    f.write("rho_center:"+str(rho_center))
    for i in range(bin_num+1):
        f.write("\n"+str(hist_ave[i])+" "+str(bins[i]))