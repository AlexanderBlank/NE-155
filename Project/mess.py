from scipy.sparse import linalg
from scipy.sparse import dia_matrix
import numpy as np


class PhysicsVals:
	"""
	D: diffusion coefficient [cm]
	Sigma_a: macroscopic absorption cross section [cm^-1]
	S: fixed neutron source [neutrons / (cm^3 s)]
	"""
	pass


M = 32 # y cell count == number of rows
N = 16 # x cell count == number of columns

width = [.2 for i in range(N)] #width of each column at each mesh x index
height = [.3 for i in range(M)]

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


#for a fixed source boundary, assign a number; assing 0 for vacuum
#for reflecting, assign 'reflecting'
#periodic not implemented
#other boundary conditions not implemented
leftBoundary = 0
rightBoundary = 'reflecting'
topBoundary = 0
bottomBoundary = 'reflecting'
#implementing periodic boundaries would require a 9-banded matrix
#so don't do anything that's too specific to 5-banded stuff
#TODO verify that periodic boundaries are paired

for x_index in range(N):
	for y_index in range(M):
		
		c = mesh[y_index][x_index]
		
		c.S = 1
		
		if (xCellCenters[x_index] - 1)**2 + (yCellCenters[y_index] - 4)**2 < 1.5**2:
			c.D, c.Sigma_a = .5, 0.05
		else:
			c.D, c.Sigma_a = .2, 0.05
		
	


#TODO verify that all mesh cells have positive D, Sigma_a, S, W, and H



BOTTOM, LEFT, CENTER, RIGHT, TOP = range(5)

dirs = [(0, -1), (-1, 0), (0, 0), (1, 0), (0, 1)]

def avgSigma(i, j):
	"""
	avergaes values of Sigma_a from cells adjacent to vertex (i, j)
	"""
	
	'''
	Sigma_a will end up being half as large on the edges, which makes sense
	for a vacumm boundary,
	but with a reflected boundary there's the same Sigma_a on the other side
	so average Sigma_a shouldn't be cut in half?
	k = min(max(0, k), M-1)
	to count local values as values across reflecting boundary
	
	
	ask Prof. Slaybaugh about this
	
	the implementation below seems to give the right answer and the alternative
	above gives the wrong answer
	'''
	
	result = 0
	
	
	for k in [i-1, i]:
		for l in [j-1, j]:
			if leftBoundary == 'periodic':
				k %= N
			if topBoundary == 'periodic':
				l %= M
			if k >= 0 and k < N and l >= 0 and l < M:
				result += mesh[l][k].Sigma_a * width[k]*height[l]/4
	
	if i == 0 and isNumeric(leftBoundary):
		result = 1
	if i == N and isNumeric(rightBoundary):
		result = 1
	if j == 0 and isNumeric(bottomBoundary):
		result = 1
	if j == M and isNumeric(topBoundary):
		result = 1
	
	return result

def isNumeric(x):
	return isinstance(x, (int, float))

def avgSource(i, j):
	"""
	averages values of source from cells adjacent to vertex (i, j)
	"""
	
	#at ambiguous corners, use the average of the source term from each surface
	
	edgeConditions = []
	
	if i == 0 and isNumeric(leftBoundary):
		edgeConditions.append(leftBoundary)
	if i == N and isNumeric(rightBoundary):
		edgeConditions.append(rightBoundary)
	if j == 0 and isNumeric(bottomBoundary):
		edgeConditions.append(bottomBoundary)
	if j == M and isNumeric(topBoundary):
		edgeConditions.append(topBoundary)
	
	if len(edgeConditions) > 0:
		return sum(edgeConditions) / len(edgeConditions)
	
	#if we're not at a fixed source edge, sum over the adjacent cells
	
	result = 0
	
	for k in [i-1, i]:
		for l in [j-1, j]:
			if leftBoundary == 'periodic':
				k %= N
			if topBoundary == 'periodic':
				l %= M
			if k >= 0 and k < N and l >= 0 and l < M:
				result += mesh[l][k].S * width[k]*height[l]/4
	
	return result
	

def matrixElement(i, j, direction):
	"""
	returns matrix element corresponding to the effect on the flux at (i, j)
	due to the flux in the adjacent cell in given direction
	"""
	
	
	#correct for off-by-one between edge-centered and cell-centered coordinates
	i_next = i
	j_next = j
	i = i-1
	j = j-1
	
	if i == -1 and isNumeric(leftBoundary):
		return 0
	if i_next == N and isNumeric(rightBoundary):
		return 0
	if j == -1 and isNumeric(bottomBoundary):
		return 0
	if j_next == M and isNumeric(topBoundary):
		return 0
	
	if i == -1:
		if leftBoundary == 'periodic':
			i = N-1 # wrap around to the right edge (negative list indexing in python would have had the same effect)
		elif leftBoundary == 'reflecting':
			if (direction == BOTTOM):
				return -(mesh[j][i_next].D*width[i_next]) / (2*height[j])
			if (direction == TOP):
				return -(mesh[j_next][i_next].D*width[i_next]) / (2*height[j_next])
	
	if j == -1:
		if bottomBoundary == 'periodic':
			i = M-1
		elif bottomBoundary == 'reflecting':
			if (direction == LEFT):
				return -(mesh[j_next][i].D*height[j_next]) / (2*width[i])
			if (direction == RIGHT):
				return -(mesh[j_next][i_next].D*height[j_next]) / (2*width[i_next])
	
	if i_next == N:
		if rightBoundary == 'periodic':
			i_next = 0 # wrap back to the left edge
		elif rightBoundary == 'reflecting':
			if (direction == BOTTOM):
				return -(mesh[j][i].D*width[i]) / (2*height[j])
			if (direction == TOP):
				return -(mesh[j_next][i].D*width[i]) / (2*height[j_next])
	
	if j_next == M:
		if topBoundary == 'periodic':
			j_next = 0 # wrap back to the bottom edge
		elif topBoundary == 'reflecting':
			if (direction == LEFT):
				return -(mesh[j][i].D*height[j]) / (2*width[i])
			if (direction == RIGHT):
				return -(mesh[j][i_next].D*height[j]) / (2*width[i_next])
	
	
	#copied from lecture notes
	if (direction == LEFT):
		return -(mesh[j][i].D*height[j] + mesh[j_next][i].D*height[j_next]) / (2*width[i])
	if (direction == RIGHT):
		return -(mesh[j][i_next].D*height[j] + mesh[j_next][i_next].D*height[j_next]) / (2*width[i_next])
	if (direction == BOTTOM):
		return -(mesh[j][i].D*width[i] + mesh[j][i_next].D*width[i_next]) / (2*height[j])
	if (direction == TOP):
		return -(mesh[j_next][i].D*width[i] + mesh[j_next][i_next].D*width[i_next]) / (2*height[j_next])
	
	print('a bad direction was passed to a(i, j, direction)')
	print('i:', i, 'j:', j, 'direction:', direction)


diagonals = [[] for d in range(5)]

#spacing for when the diagonals are built into a diagonal matrix
diagonals[RIGHT] += [0]
diagonals[TOP] += [0] * (N+1)

for j in range(M+1):
	diagonals[LEFT].append(0)
	for i in range(N+1):
		#fill with Sigma_a and subtract values for adjacent mesh cells as they're calculated
		diagonals[CENTER].append(avgSigma(i, j))
		
		for k in [BOTTOM, LEFT, RIGHT, TOP]:
			d = dirs[k]
			if i + d[0] >= 0 and i + d[0] <= N and j + d[1] >= 0 and j + d[1] <= M:
				diagonals[k].append(matrixElement(i, j, k))
				diagonals[CENTER][-1] -= diagonals[k][-1]
		
	diagonals[RIGHT].append(0)
	



#spacing for when the diagonals are built into a diagonal matrix
diagonals[LEFT] = diagonals[LEFT][1:] + [0]
diagonals[RIGHT].pop()
diagonals[BOTTOM] += [0] * (N+1)

#places rows at these indices relative to the main diagonal
offsets = [-(N+1), -1, 0, 1, N+1]

#TODO add diagonals for periodic boundaries (or TODON't)

mat = dia_matrix((np.matrix(diagonals), offsets), shape=((M+1)*(N+1), (M+1)*(N+1)))

source = []

for j in range(M+1):
	for i in range(N+1):
		source.append(avgSource(i, j))


#TODO: replace this with an iterative solver
fluxVector = linalg.spsolve(mat,  source)


flux = [fluxVector[i:i+N+1] for i in range(0, (M+1)*(N+1), N+1)]

cellAvgFlux = [[sum(flux[l][k] for k in [i, i+1] for l in [j, j+1])/4 for i in range(N)] for j in range(M)]


with open('out.txt', 'a') as out:
	print('\n\n', file=out)
	
	print('new run of diffusion solver\n', file=out)
	
	print('rows:', M, file=out)
	print('columns:', N, file=out)
	
	print('\n', file=out)
	
	print('cell widths', file=out)
	print(', '.join(str(val) for val in width), file=out)
	print('\n', file=out)
	
	print('cell heights', file=out)
	print('\n'.join(str(val) for val in reversed(height)), file=out)
	print('\n', file=out)
	
	#print all given physics values
	for key in mesh[0][0].__dict__:
		print(key, file=out)
		for row in reversed(mesh): #reversed so that it prints such that y = 0 is on the bottom
			print(', '.join(str(val.__dict__[key]) for val in row), file=out)
			
		print('\n\n', file=out)
	
	print('matrix to invert', file=out)
	for row in mat.toarray():
		print(', '.join(str(val) for val in row), file=out)
	
	print('\n\n', file=out)
	
	print('source vector', file=out)
	print('\n'.join(str(val) for val in source), file=out)
	
	print('\n\n', file=out)
	
	print('flux', file=out)
	for row in reversed(flux):
		print(', '.join(str(val) for val in row), file=out)
	
	print('\n\n', file=out)
	
	print('cell averaged flux', file=out)
	for row in reversed(cellAvgFlux):
		print(', '.join(str(val) for val in row), file=out)
	
	print('\n\n', file=out)

import matplotlib.pyplot as plt

plt.pcolormesh(xVals, yVals, cellAvgFlux)
plt.show()

"""
from math import cosh, sqrt

def fixedExact(x):
	return 8/0.2*(1 - cosh(sqrt(0.2/1)*x) / cosh(sqrt(0.2/1)*4))


fluxAnalytical = [[fixedExact(x-4) for x in xVals] for y in yVals]

cellAvgFluxAnalytical = [[sum(fluxAnalytical[l][k] for k in [i, i+1] for l in [j, j+1])/4 for i in range(N)] for j in range(M)]


plt.pcolormesh(xVals, yVals, cellAvgFluxAnalytical)
plt.show()

diff = [[cellAvgFlux[j][i]-cellAvgFluxAnalytical[j][i] for i in range(N)] for j in range(M)]

with open('out.txt', 'a') as out:
	print('analytical cell averaged flux', file=out)
	for row in reversed(cellAvgFluxAnalytical):
		print(', '.join(str(val) for val in row), file=out)
	
	print('\n\n', file=out)
	
	print('diff', file=out)
	for row in reversed(diff):
		print(', '.join(str(val) for val in row), file=out)
"""