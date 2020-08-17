import numpy as np

start_step = input("Start Step:")
dump_step = input("Dump Timestep:")
file_num = input("Number of Dumpfile:")

grid_x = 100
grid_y = 100
grid_z = 70

x0 = 100.0
y0 = 100.0
z0 = 180.0

dx = 1.0

bin_num = 1000

hist_ave = np.zeros(bin_num+1)

interface_z = np.zeros((3, grid_x*grid_y))

print("Histgram")

for num in range(int(file_num)):
	step = int(start_step) + int(num) * int(dump_step)
	lineCount = 0
	
	rho = np.zeros(grid_x*grid_y*grid_z)
	
	print("Now Step:"+str(step))

	for line in open("rho_D_" + str(step) +".txt","r"):
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

print("Interface")

for num in range(int(file_num)):
	step = int(start_step) + int(num) * int(dump_step)
	itr = 0
	lineCount = 0
	rho_z = 0.0
	rho_z_old = 0.0
	chack = 0
	
	print("Now Step:"+str(step))

	for line in open("rho_D_" + str(step) +".txt","r"):
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

	with open("interface_Z_" + str(step) + ".xyz",mode="w") as f:
		f.write(str(grid_x*grid_y) + "\nInterface")
		for i in range(grid_x*grid_y):
			f.write("\n"+"X "+str(interface_z[0, i])+" "+str(interface_z[1, i])+" "+str(interface_z[2, i]))
