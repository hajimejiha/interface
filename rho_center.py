import numpy as np

start_step = int(input("Start Step:"))
dump_step = int(input("Dump Timestep:"))
file_num = int(input("Number of Dumpfile:"))

grid_x = 260
grid_y = 265
grid_z = 240

bin_num = 1000
bins2 = np.zeros(bin_num)
hist_ave = np.zeros(bin_num)

rho_max = 0
rho_min = float("inf")



print("\n--- rho_max rho_min ---")
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

    if np.max(rho) > rho_max:
        rho_max = np.max(rho) 
    if np.min(rho) < rho_min:
        rho_min = np.min(rho) 



print("\n--- rho_center ---")
for num in range(file_num):
    step = start_step + num * dump_step
    print("Now Step:"+str(step))
    
    lineCount = 0
    
    
    for line in open("rho_D." + str(step) +".txt","r"):
        lineCount += 1

        if lineCount < 2:
            continue

        data = line.split()

        rho[lineCount - 2] = data[3]

    hist, bins = np.histogram(rho, bins=bin_num, range=(rho_min, rho_max))

    if num == 0:
        for i in range(bin_num):
            bins2[i] = (bins[i] + bins[i+1]) / 2

    for i in range(bin_num):
        hist_ave[i] += hist[i]
        
"""
    with open("histgram."+str(step)+".rho_D.txt",mode="w") as f:
        for i in range(bin_num):
            f.write("\n"+str(hist[i])+" "+str(bins[i]))
        f.write("\n0 "+str(bins[bin_num]))
"""

for i in range(bin_num):
    hist_ave[i] /= file_num

hist_ave2 = np.delete(hist_ave, 0)
index_liq = np.where(hist_ave2 == np.max(hist_ave2))

rho_center = (bins2[int(index_liq[0] // 2)] + bins2[int(-(-index_liq[0] // 2))]) / 2

with open("histgram_average_rho_D.txt",mode="w") as f:
    f.write("rho_center:%f rho_max:%f rho_min:%f" % (rho_center, rho_max, rho_min))
    for i in range(bin_num):
        f.write("\n%f %.3f" % (bins2[i], hist_ave[i]))
