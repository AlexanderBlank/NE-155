from diffusion_solver import *

M = 111 # y cell count == number of rows
N = 128 # x cell count == number of columns

width = [.02 if i in range(50, 90) else 0.123 for i in range(N)] #width of each column at each mesh x index
height = [.05 if i in range(0, 60) else .1 for i in range(M)]

mesh = [[PhysicsVals() for i in range(N)] for j in range(M)]

def integral(arr):
	result = [0]
	for x in arr:
		result.append(result[-1] + x)
	return result

xVals = integral(width)
yVals = integral(height)

xCellCenters = [(xVals[i]+xVals[i+1])/2 for i in range(N)]
yCellCenters = [(yVals[i]+yVals[i+1])/2 for i in range(M)]

for y_index in range(M):
	for x_index in range(N):
		
		c = mesh[y_index][x_index]
		
		c.D, c.Sigma_a, c.S = 1, 0.1, 1
		
		if (yCellCenters[y_index] // 1) % 2 != 0:
			c.D, c.Sigma_a, c.S = 0.1, 0.1, 0
		
		if (xCellCenters[x_index] - 6.5)**2 + (yCellCenters[y_index] - 4)**2 < 2**2:
			c.D, c.Sigma_a, c.S = 0.5, 0.4, 4
		
		if abs(xCellCenters[x_index] - 9) + (yCellCenters[y_index] - 3) < 3:
			c.D, c.Sigma_a, c.S = .2, 1.6, 3
		
	

mesh[30][50] = PhysicsVals(.7, 0.2, 1600)

import matplotlib.pyplot as plt
plt.set_cmap('gray')

D = [[val.D for val in row] for row in mesh]
plt.pcolormesh(xVals, yVals, D)
plt.colorbar()
plt.show()

outFile = open('demoOutput.txt', 'a')

solver = diffusion_solver(mesh, width, height, leftBoundary=0, rightBoundary='ref',\
	topBoundary=5, bottomBoundary='ref', outputStream=outFile)

solver.printGeneralInfo()
solver.printInput()

import time
start_time = time.time()
solver.solve()
#print(solver.flux)
print("%s seconds" % (time.time() - start_time))

solver.printFlux()
solver.printCellAvgFlux()

plt.pcolormesh(xVals, yVals, solver.cellAvgFlux())
plt.colorbar()
plt.show()
