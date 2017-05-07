from diffusion_solver import *
from random import uniform as randRange

def compareAnalytical1DFixed():
	"""
	checks the solver output against the analytical solution to the diffusion equation
	for vacuum boundary conditions with uniform source and physics values in one dimension
	"""
	with open('compareAnalytical1DFixed.txt', 'a') as out:
		#use randomly selected fixed values
		D_0 = randRange(0.1, 20)
		Sigma_a_0 = randRange(0.1, 20)
		S_0 = randRange(0, 20)
		
		M = 128 # y cell count == number of rows
		N = 128 # x cell count == number of columns
		
		cellWidth = randRange(0.01, .2)
		
		width = [cellWidth for i in range(N)] #width of each column at each mesh x index
		height = [0.5 for i in range(M)]
		
		mesh = [[PhysicsVals(D_0, Sigma_a_0, S_0) for i in range(N)] for j in range(M)]
		
		def integral(arr):
			result = [0] * (len(arr)+1)
			for i in range(len(arr)):
				result[i+1] = result[i] + arr[i]
			return result
		
		xVals = integral(width)
		yVals = integral(height)
		
		solver = diffusion_solver(mesh, width, height,\
			leftBoundary=0, rightBoundary=0, topBoundary='ref', bottomBoundary='ref',\
			outputStream=out)
		
		solver.printInput()
		
		import time
		start_time = time.time()
		solver.solve()
		print("%s seconds" % (time.time() - start_time), file=out)
		
		solver.printFlux()
		solver.printCellAvgFlux()
		
		cellAvgFlux = solver.cellAvgFlux()
		
		from math import cosh, sqrt
		
		#from homework 6
		def fixedExact(x):
			return S_0/Sigma_a_0*(1 - cosh(sqrt(Sigma_a_0/D_0)*x) / cosh(sqrt(Sigma_a_0/D_0)*(xVals[-1]/2)))
		
		
		fluxAnalytical = [[fixedExact(x-(xVals[-1]/2)) for x in xVals] for y in yVals]
		
		cellAvgFluxAnalytical = [[sum(fluxAnalytical[l][k] for k in [i, i+1] for l in [j, j+1])/4 for i in range(N)] for j in range(M)]
		
		diff = [[solver.flux[j][i]-fluxAnalytical[j][i] for i in range(N)] for j in range(M)]
		
		err = (sum(sum(d**2 for d in row) for row in diff))**0.5 / (M*N)
		
		print('root-mean-square error from analytical solution:', err, file=out)
		
		
		import matplotlib.pyplot as plt
		plt.set_cmap('gray')
		
		plt.pcolormesh(xVals, yVals, cellAvgFlux)
		plt.colorbar()
		plt.show()
		
		plt.pcolormesh(xVals, yVals, cellAvgFluxAnalytical)
		plt.colorbar()
		plt.show()
		
		return err
	

def compareAnalytical1DCos():
	"""
	checks the solver output against the analytical solution to the diffusion equation
	for vacuum boundary conditions with cos(x) source and physics values in one dimension
	"""
	with open('compareAnalytical1DCos.txt', 'a') as out:
		#use randomly selected fixed values
		D_0 = randRange(0.1, 20)
		Sigma_a_0 = randRange(0.1, 20)
		S_0 = randRange(0, 20)
		
		M = 128 # y cell count == number of rows
		N = 128 # x cell count == number of columns
		
		cellWidth = randRange(0.005, 0.015) #limited so that source is never negative
		
		width = [cellWidth for i in range(N)] #width of each column at each mesh x index
		height = [0.5 for i in range(M)]
		
		def integral(arr):
			result = [0] * (len(arr)+1)
			for i in range(len(arr)):
				result[i+1] = result[i] + arr[i]
			return result
		
		xVals = integral(width)
		yVals = integral(height)
		
		xCellCenters = [(xVals[i]+xVals[i+1])/2 for i in range(N)]
		yCellCenters = [(yVals[i]+yVals[i+1])/2 for i in range(M)]
		
		from math import cos
		
		mesh = [[PhysicsVals(D_0, Sigma_a_0, S_0*cos((xCellCenters[i] - xCellCenters[-1]/2))) for i in range(N)] for j in range(M)]
		
		solver = diffusion_solver(mesh, width, height,\
			leftBoundary=0, rightBoundary=0, topBoundary='ref', bottomBoundary='ref',\
			outputStream=out)
		
		solver.printInput()
		
		import time
		start_time = time.time()
		solver.solve()
		print("%s seconds" % (time.time() - start_time), file=out)
		
		solver.printFlux()
		solver.printCellAvgFlux()
		
		cellAvgFlux = solver.cellAvgFlux()
		
		from math import cosh, sqrt
		
		#from homework 6
		def fixedExact(x):
			return S_0/(D_0 + Sigma_a_0)*(cos(x) - cos((xVals[-1]/2))*cosh(sqrt(Sigma_a_0/D_0)*x) / cosh(sqrt(Sigma_a_0/D_0)*(xVals[-1]/2)))
		
		
		fluxAnalytical = [[fixedExact(x-(xVals[-1]/2)) for x in xVals] for y in yVals]
		
		cellAvgFluxAnalytical = [[sum(fluxAnalytical[l][k] for k in [i, i+1] for l in [j, j+1])/4 for i in range(N)] for j in range(M)]
		
		diff = [[solver.flux[j][i]-fluxAnalytical[j][i] for i in range(N)] for j in range(M)]
		
		err = (sum(sum(d**2 for d in row) for row in diff))**0.5 / (M*N)
		
		print('root-mean-square error from analytical solution:', err, file=out)
		
		
		import matplotlib.pyplot as plt
		plt.set_cmap('gray')
		
		plt.pcolormesh(xVals, yVals, cellAvgFlux)
		plt.colorbar()
		plt.show()
		
		plt.pcolormesh(xVals, yVals, cellAvgFluxAnalytical)
		plt.colorbar()
		plt.show()
		
		return err
	

def compareAnalytical2DFixed():
	"""
	checks the solver output against the analytical solution to the diffusion equation
	for vacuum boundary conditions with uniform source and physics values in two dimensions
	"""
	with open('compareAnalytical2DFixed.txt', 'a') as out:
		#use randomly selected fixed values
		D_0 = randRange(0.1, 20)
		Sigma_a_0 = randRange(0.1, 20)
		S_0 = randRange(0, 20)
		
		M = 128 # y cell count == number of rows
		N = 128 # x cell count == number of columns
		
		cellWidth = randRange(0.01, .2)
		
		width = [cellWidth for i in range(N)] #width of each column at each mesh x index
		height = [0.5 for i in range(M)]
		
		mesh = [[PhysicsVals(D_0, Sigma_a_0, S_0) for i in range(N)] for j in range(M)]
		
		def integral(arr):
			result = [0] * (len(arr)+1)
			for i in range(len(arr)):
				result[i+1] = result[i] + arr[i]
			return result
		
		xVals = integral(width)
		yVals = integral(height)
		
		solver = diffusion_solver(mesh, width, height,\
			leftBoundary=0, rightBoundary=0, topBoundary=0, bottomBoundary=0,\
			outputStream=out)
		
		solver.printInput()
		
		import time
		start_time = time.time()
		solver.solve()
		print("%s seconds" % (time.time() - start_time), file=out)
		
		solver.printFlux()
		solver.printCellAvgFlux()
		
		cellAvgFlux = solver.cellAvgFlux()
		
		from math import cosh, sqrt
		
		#based on separation of variables
		def fixedExact_x(x):
			return sqrt(S_0/Sigma_a_0)*(1 - cosh(sqrt(Sigma_a_0/D_0)*x) / cosh(sqrt(Sigma_a_0/D_0)*(xVals[-1]/2)))
		
		def fixedExact_y(y):
			return sqrt(S_0/Sigma_a_0)*(1 - cosh(sqrt(Sigma_a_0/D_0)*y) / cosh(sqrt(Sigma_a_0/D_0)*(yVals[-1]/2)))
		
		fluxAnalytical = [[fixedExact_x(x-(xVals[-1]/2))*fixedExact_y(y-(yVals[-1]/2)) for x in xVals] for y in yVals]
		
		cellAvgFluxAnalytical = [[sum(fluxAnalytical[l][k] for k in [i, i+1] for l in [j, j+1])/4 for i in range(N)] for j in range(M)]
		
		diff = [[solver.flux[j][i]-fluxAnalytical[j][i] for i in range(N)] for j in range(M)]
		
		err = (sum(sum(d**2 for d in row) for row in diff))**0.5 / (M*N)
		
		print('root-mean-square error from analytical solution:', err, file=out)
		
		
		import matplotlib.pyplot as plt
		plt.set_cmap('gray')
		
		plt.pcolormesh(xVals, yVals, cellAvgFlux)
		plt.colorbar()
		plt.show()
		
		plt.pcolormesh(xVals, yVals, cellAvgFluxAnalytical)
		plt.colorbar()
		plt.show()
		
		return err
	


def verifyReflection():
	"""
	runs the solver on the same problem, once using reflecting boundary conditions,
	and once by building a mesh four times as large with fixed boundary conditions
	and compares the results from the two
	"""
	with open('verifyReflection.txt', 'a') as out:
		#start with bottom left quadrant with reflecting boundaries on the top and right
		
		M = int(randRange(32, 64)) # y cell count == number of rows
		N = int(randRange(32, 64)) # x cell count == number of columns
		
		width = [randRange(0.01, 4) for i in range(N)] #width of each column at each mesh x index
		height = [randRange(0.01, 4) for i in range(M)]
		
		#populate with random values
		mesh = [[PhysicsVals(randRange(0, 20), randRange(0, 20), randRange(0, 20)) for i in range(N)] for j in range(M)]
		
		solver = diffusion_solver(mesh, width, height,\
			leftBoundary=randRange(0, 5), rightBoundary='ref', topBoundary='ref', bottomBoundary=randRange(0, 5),\
			outputStream=out)
		
		solver.printInput()
		
		import time
		start_time = time.time()
		solver.solve()
		print("%s seconds" % (time.time() - start_time), file=out)
		
		solver.printFlux()
		solver.printCellAvgFlux()
		
		#now for the full solution, without using reflection
		
		M_full = M*2
		N_full = N*2
		
		width_full = width + list(reversed(width))
		height_full = height + list(reversed(height))
		
		mesh2 = [row + list(reversed(row)) for row in mesh]
		
		mesh_full = mesh2 + list(reversed(mesh2))
		
		
		solver_full = diffusion_solver(mesh_full, width_full, height_full,\
			leftBoundary=solver.leftBoundary, rightBoundary=solver.leftBoundary,\
			topBoundary=solver.bottomBoundary, bottomBoundary=solver.bottomBoundary,\
			outputStream=out)
		
		start_time = time.time()
		solver_full.solve()
		print("%s seconds" % (time.time() - start_time), file=out)
		
		solver_full.printFlux()
		solver_full.printCellAvgFlux()
		
		flux2 = [row + list(reversed(row))[1:] for row in solver.flux]
		flux_expanded = flux2 + list(reversed(flux2))[1:]
		
		flux_full = solver_full.flux
		
		diff = [[flux_expanded[j][i]-flux_full[j][i] for i in range(N*2)] for j in range(M*2)]
		
		err = (sum(sum(d**2 for d in row) for row in diff))**0.5 / (M*N*4)
		
		print('root-mean-square error between full solution and expanded quarter-solution:', err, file=out)
		
		def integral(arr):
			result = [0] * (len(arr)+1)
			for i in range(len(arr)):
				result[i+1] = result[i] + arr[i]
			return result
		
		xVals = integral(width_full)
		yVals = integral(height_full)
		
		import matplotlib.pyplot as plt
		plt.set_cmap('gray')
		
		plt.pcolormesh(xVals, yVals, solver_full.cellAvgFlux())
		plt.colorbar()
		plt.show()
		
		cellAvgFlux2 = [row + list(reversed(row)) for row in solver.cellAvgFlux()]
		cellAvgFlux_expanded = cellAvgFlux2 + list(reversed(cellAvgFlux2))
		
		plt.pcolormesh(xVals, yVals, cellAvgFlux_expanded)
		plt.colorbar()
		plt.show()
		
		return err
	

def main():
	print("running comparison to analytical solution for fixed source in one dimension...")
	err = compareAnalytical1DFixed()
	print('root-mean-square error from analytical solution:', err)
	print('images should be identical')
	print('error should be as small as possible, expected around', 1e-6)
	
	print("running comparison to analytical solution for cos(x) source in one dimension...")
	err = compareAnalytical1DCos()
	print('root-mean-square error from analytical solution:', err)
	print('images should be identical')
	print('error should be as small as possible, expected around', 1e-6)
	
	print("running comparison to analytical solution for fixed source in two dimensions...")
	err = compareAnalytical2DFixed()
	print('root-mean-square error from analytical solution:', err)
	print('images should be identical')
	print('error should be as small as possible, expected around', 1e-6)
	
	print("verfifying reflecting boundary conditions...")
	err = verifyReflection()
	print('root-mean-square error between full solution and expanded quarter-solution:', err)
	print('images should be identical')
	print('error should be as small as possible, expected around', 1e-11)
	
	
	print("tests complete")
	


if __name__ == "__main__":
	main()