"""
diffusion_solver.py
"""

class PhysicsVals:
	"""
	D: diffusion coefficient [cm]
	Sigma_a: macroscopic absorption cross section [cm^-1]
	S: fixed neutron source [neutrons / (cm^3 s)]
	"""
	def __init__(self, D=0, Sigma_a=0, S=0):
		self.D = D
		self.Sigma_a = Sigma_a
		self.S = S

class diffusion_solver:
	"""
	solves the neutron diffusion equation on a two-dimensional mesh
	"""
	
	BOTTOM, LEFT, CENTER, RIGHT, TOP = range(5)
	dirs = [(0, -1), (-1, 0), (0, 0), (1, 0), (0, 1)]
	REFLECTING = 'reflecting'
	PERIODIC = 'periodic'
	
	@classmethod
	def standardizeName(cls, x):
		if x in ['ref', 'reflecting', diffusion_solver.REFLECTING]:
			return diffusion_solver.REFLECTING
		if x in ['per', 'periodic', diffusion_solver.PERIODIC]:
			return diffusion_solver.PERIODIC
		
		return x
	
	def __init__(self, physicsValueMesh, widths, heights, \
				 leftBoundary=0, rightBoundary=0, topBoundary=0, bottomBoundary=0, \
				outputStream=None):
		"""
		physicsValueMesh is an M by N row-major array of PhysicsVals, ordered bottom to top and left to right
		widths is an array of length N whose entries correspond to mesh cell widths, ordered left to right
		heights is an array of length M whose entries correspond to mesh cell heights, ordered bottom to top
		
		for a fixed source boundary, assign a number; assign 0 for vacuum
		for reflecting, assign 'reflecting'
		periodic not implemented
		other boundary conditions not implemented
		
		outputStream is where the solver will print to, leave as None for standard output
		"""
		self.mesh = physicsValueMesh
		self.width = widths
		self.height = heights
		
		self.leftBoundary = self.standardizeName(leftBoundary)
		self.rightBoundary = self.standardizeName(rightBoundary)
		self.topBoundary = self.standardizeName(topBoundary)
		self.bottomBoundary = self.standardizeName(bottomBoundary)
		
		self.N = len(self.width)
		self.M = len(self.height)
		
		self.outStream = outputStream
		
	
	@classmethod
	def isNumeric(cls, x):
		return isinstance(x, (int, float))
	
	def avgSigma(self, i, j):
		"""
		avergaes values of Sigma_a from cells adjacent to vertex (i, j)
		"""
		
		result = 0
		
		for k in [i-1, i]:
			for l in [j-1, j]:
				if self.leftBoundary == self.PERIODIC:
					k %= self.N
				if self.topBoundary == self.PERIODIC:
					l %= self.M
				if k >= 0 and k < self.N and l >= 0 and l < self.M:
					result += self.mesh[l][k].Sigma_a * self.width[k]*self.height[l]/4
		
		#if we're at a fixed source boundary, the flux should equal the source at the boundary
		if i == 0 and self.isNumeric(self.leftBoundary):
			result = 1
		if i == self.N and self.isNumeric(self.rightBoundary):
			result = 1
		if j == 0 and self.isNumeric(self.bottomBoundary):
			result = 1
		if j == self.M and self.isNumeric(self.topBoundary):
			result = 1
		
		return result
	
	def avgSource(self, i, j):
		"""
		averages values of source from cells adjacent to vertex (i, j)
		"""
		
		#at ambiguous corners, use the average of the source term from each surface
		
		edgeConditions = []
		
		if i == 0 and self.isNumeric(self.leftBoundary):
			edgeConditions.append(self.leftBoundary)
		if i == self.N and self.isNumeric(self.rightBoundary):
			edgeConditions.append(self.rightBoundary)
		if j == 0 and self.isNumeric(self.bottomBoundary):
			edgeConditions.append(self.bottomBoundary)
		if j == self.M and self.isNumeric(self.topBoundary):
			edgeConditions.append(self.topBoundary)
		
		if len(edgeConditions) > 0:
			return sum(edgeConditions) / len(edgeConditions)
		
		#if we're not at a fixed source edge, sum over the adjacent cells
		
		result = 0
		
		for k in [i-1, i]:
			for l in [j-1, j]:
				if self.leftBoundary == self.PERIODIC:
					k %= self.N
				if self.topBoundary == self.PERIODIC:
					l %= self.M
				if k >= 0 and k < self.N and l >= 0 and l < self.M:
					result += self.mesh[l][k].S * self.width[k]*self.height[l]/4
		
		return result
	
	def matrixElement(self, i, j, direction):
		"""
		returns matrix element corresponding to the effect on the flux at (i, j)
		due to the flux in the adjacent cell in given direction
		"""
		
		
		#correct for off-by-one between edge-centered and cell-centered coordinates
		i_next = i
		j_next = j
		i = i-1
		j = j-1
		
		#if we're at a fixed source boundary, the flux should equal the source at the boundary
		if i == -1 and self.isNumeric(self.leftBoundary):
			return 0
		if i_next == self.N and self.isNumeric(self.rightBoundary):
			return 0
		if j == -1 and self.isNumeric(self.bottomBoundary):
			return 0
		if j_next == self.M and self.isNumeric(self.topBoundary):
			return 0
		
		if i == -1:
			if self.leftBoundary == self.PERIODIC:
				i = self.N-1 # wrap around to the right edge (negative list indexing in python would have had the same effect)
			elif self.leftBoundary == self.REFLECTING:
				if (direction == self.BOTTOM):
					return -(self.mesh[j][i_next].D*self.width[i_next]) / (2*self.height[j])
				if (direction == self.TOP):
					return -(self.mesh[j_next][i_next].D*self.width[i_next]) / (2*self.height[j_next])
		
		if j == -1:
			if self.bottomBoundary == self.PERIODIC:
				i = self.M-1
			elif self.bottomBoundary == self.REFLECTING:
				if (direction == self.LEFT):
					return -(self.mesh[j_next][i].D*self.height[j_next]) / (2*self.width[i])
				if (direction == self.RIGHT):
					return -(self.mesh[j_next][i_next].D*self.height[j_next]) / (2*self.width[i_next])
		
		if i_next == self.N:
			if self.rightBoundary == self.PERIODIC:
				i_next = 0 # wrap back to the left edge
			elif self.rightBoundary == self.REFLECTING:
				if (direction == self.BOTTOM):
					return -(self.mesh[j][i].D*self.width[i]) / (2*self.height[j])
				if (direction == self.TOP):
					return -(self.mesh[j_next][i].D*self.width[i]) / (2*self.height[j_next])
		
		if j_next == self.M:
			if self.topBoundary == self.PERIODIC:
				j_next = 0 # wrap back to the bottom edge
			elif self.topBoundary == self.REFLECTING:
				if (direction == self.LEFT):
					return -(self.mesh[j][i].D*self.height[j]) / (2*self.width[i])
				if (direction == self.RIGHT):
					return -(self.mesh[j][i_next].D*self.height[j]) / (2*self.width[i_next])
		
		
		#copied from lecture notes
		if (direction == self.LEFT):
			return -(self.mesh[j][i].D*self.height[j] + self.mesh[j_next][i].D*self.height[j_next]) / (2*self.width[i])
		if (direction == self.RIGHT):
			return -(self.mesh[j][i_next].D*self.height[j] + self.mesh[j_next][i_next].D*self.height[j_next]) / (2*self.width[i_next])
		if (direction == self.BOTTOM):
			return -(self.mesh[j][i].D*self.width[i] + self.mesh[j][i_next].D*self.width[i_next]) / (2*self.height[j])
		if (direction == self.TOP):
			return -(self.mesh[j_next][i].D*self.width[i] + self.mesh[j_next][i_next].D*self.width[i_next]) / (2*self.height[j_next])
		
		print('a bad direction was passed to matrixElement(i, j, direction)', file=self.outStream)
		print('i:', i, 'j:', j, 'direction:', direction, file=self.outStream)
	
	def buildMatrix(self):
		#5 diagonals, one for the flux at a point and 4 for the neighbors
		self.diagonals = [[] for d in range(5)]
		
		#spacing for when the diagonals are built into a diagonal matrix
		self.diagonals[self.RIGHT] += [0]
		self.diagonals[self.TOP] += [0] * (self.N+1)
		
		for j in range(self.M+1):
			self.diagonals[self.LEFT].append(0)
			for i in range(self.N+1):
				#fill with Sigma_a and subtract values for adjacent mesh cells as they're calculated
				self.diagonals[self.CENTER].append(self.avgSigma(i, j))
				
				for k in [self.BOTTOM, self.LEFT, self.RIGHT, self.TOP]:
					d = self.dirs[k]
					if i + d[0] >= 0 and i + d[0] <= self.N and j + d[1] >= 0 and j + d[1] <= self.M:
						self.diagonals[k].append(self.matrixElement(i, j, k))
						self.diagonals[self.CENTER][-1] -= self.diagonals[k][-1]
			
			self.diagonals[self.RIGHT].append(0)
		
		
		#spacing for when the diagonals are built into a diagonal matrix
		self.diagonals[self.LEFT] = self.diagonals[self.LEFT][1:] + [0]
		self.diagonals[self.RIGHT].pop()
		self.diagonals[self.BOTTOM] += [0] * (self.N+1)
		
		#places rows at these indices relative to the main diagonal
		offsets = [-(self.N+1), -1, 0, 1, self.N+1]
	
	def buildSource(self):
		self.source = []
		
		for j in range(self.M+1):
			for i in range(self.N+1):
				self.source.append(self.avgSource(i, j))
	
	def getCoarseSolver(self):
		'''
		returns a new solver by merging adjacent mesh cells
		'''
		coarseN = self.N//2 + self.N%2
		coarseM = self.M//2 + self.M%2
		
		coarseWidth = [self.width[2*i] + self.width[2*i+1] for i in range(coarseN-self.N%2)]
		if self.N%2 == 1:
			coarseWidth.append(self.width[-1])
		
		coarseHeight = [self.height[2*i] + self.height[2*i+1] for i in range(coarseM-self.M%2)]
		if self.M%2 == 1:
			coarseHeight.append(self.height[-1])
		
		coarseMesh = [[PhysicsVals() for i in range(coarseN)] for j in range(coarseM)]
		
		for j in range(coarseM-self.M%2):
			for i in range(coarseN-self.N%2):
				
				c = coarseMesh[j][i]
				
				c.D, c.Sigma_a, c.S = 0, 0, 0
				
				for j2, i2 in [(j*2, i*2), (j*2+1, i*2), (j*2, i*2+1), (j*2+1, i*2+1)]:
					V = self.height[j2] * self.width[i2]
					c.D += self.mesh[j2][i2].D * V
					c.Sigma_a += self.mesh[j2][i2].Sigma_a * V
					c.S += self.mesh[j2][i2].S * V
				
				c.D /= coarseHeight[j]*coarseWidth[i]
				c.Sigma_a /= coarseHeight[j]*coarseWidth[i]
				c.S /= coarseHeight[j]*coarseWidth[i]
			
		
		if self.M%2 == 1:
			for i in range(coarseN-self.N%2):
				
				c = coarseMesh[-1][i]
				
				c.D, c.Sigma_a, c.S = 0, 0, 0
				
				for j2, i2 in [(-1, i*2), (-1, i*2+1)]:
					V = self.height[j2] * self.width[i2]
					c.D += self.mesh[j2][i2].D * V
					c.Sigma_a += self.mesh[j2][i2].Sigma_a * V
					c.S += self.mesh[j2][i2].S * V
				
				c.D /= coarseHeight[-1]*coarseWidth[i]
				c.Sigma_a /= coarseHeight[-1]*coarseWidth[i]
				c.S /= coarseHeight[-1]*coarseWidth[i]
			
		
		if self.N%2 == 1:
			for j in range(coarseM-self.M%2):
				
				c = coarseMesh[j][-1]
				
				c.D, c.Sigma_a, c.S = 0, 0, 0
				
				for j2, i2 in [(j*2, -1), (j*2+1, -1)]:
					V = self.height[j2] * self.width[i2]
					c.D += self.mesh[j2][i2].D * V
					c.Sigma_a += self.mesh[j2][i2].Sigma_a * V
					c.S += self.mesh[j2][i2].S * V
				
				c.D /= coarseHeight[j]*coarseWidth[-1]
				c.Sigma_a /= coarseHeight[j]*coarseWidth[-1]
				c.S /= coarseHeight[j]*coarseWidth[-1]
			
		
		if self.N%2 == 1 and self.M%2 == 1:
			c = coarseMesh[-1][-1]
			c.D = self.mesh[-1][-1].D
			c.Sigma_a = self.mesh[-1][-1].Sigma_a 
			c.S = self.mesh[-1][-1].S
		
		coarseSolver = diffusion_solver(physicsValueMesh=coarseMesh, widths=coarseWidth, heights=coarseHeight,\
			 leftBoundary=self.leftBoundary, rightBoundary=self.rightBoundary, topBoundary=self.topBoundary, bottomBoundary=self.bottomBoundary,\
			outputStream=self.outStream)
		
		return coarseSolver
		
	
	def getInterpolatedFlux(self, coarseFlux):
		'''
		returns interpolated flux from flux for coarseSolver
		'''
		
		flux = [0 for i in range((self.M+1)*(self.N+1))]
		
		for j in range(self.M+1-self.M%2):
			for i in range(self.N+1-self.N%2):
				if i%2 == 0 and j%2 == 0:
					f = coarseFlux[j//2][i//2]
				elif i%2 == 1 and j%2 == 0:
					widthRatio = self.width[i-1] / (self.width[i-1] + self.width[i])
					f = coarseFlux[j//2][i//2]*(1-widthRatio) + coarseFlux[j//2][i//2+1]*widthRatio
				elif i%2 == 0 and j%2 == 1:
					heightRatio = self.height[j-1] / (self.height[j-1] + self.height[j])
					f = coarseFlux[j//2][i//2]*(1-heightRatio) + coarseFlux[(j//2+1)][i//2]*heightRatio
				else: #bilinear interpolation
					widthRatio = self.width[i-1] / (self.width[i-1] + self.width[i])
					heightRatio = self.height[j-1] / (self.height[j-1] + self.height[j])
					xDiff = coarseFlux[j//2][i//2+1] - coarseFlux[j//2][i//2]
					yDiff = coarseFlux[j//2+1][i//2] - coarseFlux[j//2][i//2]
					diagTerm = coarseFlux[j//2+1][i//2+1] + coarseFlux[j//2][i//2] - coarseFlux[j//2][i//2+1] - coarseFlux[j//2+1][i//2]
					f = coarseFlux[j//2][i//2] + xDiff*widthRatio + yDiff*heightRatio + diagTerm*widthRatio*heightRatio
				
				flux[j*(self.N+1) + i] = f
		
		if self.M%2 == 1:
			for i in range(self.N+1-self.N%2):
				if i%2 == 0:
					f = coarseFlux[(self.M+1)//2][i//2]
				else:
					widthRatio = self.width[i-1] / (self.width[i-1] + self.width[i])
					f = coarseFlux[(self.M+1)//2][i//2]*(1-widthRatio) + coarseFlux[(self.M+1)//2][i//2+1]*widthRatio
				
				flux[j*(self.N+1) + i] = f
		if self.N%2 == 1:
			for j in range(self.M+1-self.M%2):
				if j%2 == 0:
					f = coarseFlux[j//2][(self.N+1)//2]
				else:
					heightRatio = self.height[j-1] / (self.height[j-1] + self.height[j])
					f = coarseFlux[j//2][(self.N+1)//2]*(1-heightRatio) + coarseFlux[j//2+1][(self.N+1)//2]*heightRatio
				
				flux[j*(self.N+1) + i] = f
		
		if self.N%2 == 1 and self.M%2 == 1:
			flux[-1] = coarseFlux[-1]
		
		return flux
		
	
	def fiveBandedSOR(self, omega=1.95, tolerance=1e-7, maxIterations=100000):
		"""
		solves for the flux using red-black SOR
		solves on a coarser grid and interpolates before solving this grid if the matrix is large enough
		slow compared to scipy sparse solver
		"""
		
		if len(self.source) > 2**12:
			newSolver = self.getCoarseSolver()
			x_coarse = newSolver.solve(tolerance=tolerance**0.5)
			x = self.getInterpolatedFlux(x_coarse)
		else:
			x = [0 for i in range(len(self.source))]
		
		toleranceSq = tolerance**2
		iterations = 0
		errSq = toleranceSq+1
		
		#localizing variables to improve performance
		bottom, left, center, right, top = self.diagonals
		width = self.N+1
		source = self.source
		
		while (errSq > toleranceSq and iterations < maxIterations):
			errSq = 0
			
			sumTerm = (source[0]
				-right[1]*x[1]
				-top[width]*x[width])
			
			increment = omega*(sumTerm / center[0] - x[0])
			errSq += increment**2
			x[0] += increment
			
			for i in range(1, width):
				sumTerm = (source[i]
					-left[i-1]*x[i-1]
					-right[i+1]*x[i+1]
					-top[i+width]*x[i+width])
				
				increment = omega*(sumTerm / center[i] - x[i])
				errSq += increment**2
				x[i] += increment
			
			
			#red-black SOR: first solve for even points, then solve for odd points using new data for even points
			for parity in [0, 1]:
				for j in range(1, self.M):
					offset = 0 if j%2 == parity else 1 #offset for selecting even and then odd points
					
					b = width*(j-1) + offset #bottom index
					l = width*j - 1 + offset #left index
					r = width*j + 1 + offset #right index
					t = width*(j+1) + offset #top index
					
					for i in range(width*j + offset, width*(j+1), 2): #start of performance-critical section
						increment = omega*((source[i]
							-bottom[b]*x[b]
							-left[l]*x[l]
							-right[r]*x[r]
							-top[t]*x[t]) / center[i] - x[i])
						errSq += increment**2
						x[i] += increment
						
						b += 2
						l += 2
						r += 2
						t += 2
					#end of performance-critical section
				
				
			
			for i in range(len(x)-width, len(x)-1):
				sumTerm = (source[i]
					-bottom[i-width]*x[i-width]
					-left[i-1]*x[i-1]
					-right[i+1]*x[i+1])
				
				increment = omega*(sumTerm / center[i] - x[i])
				errSq += increment**2
				x[i] += increment
			
			i = len(x) - 1
			sumTerm = (source[i]
				-bottom[i-width]*x[i-width]
				-left[i-1]*x[i-1])
			
			increment = omega*(sumTerm / center[i] - x[i])
			errSq += increment**2
			x[i] += increment
			
			iterations += 1
		
		print('iterations:', iterations, file=self.outStream)
		if iterations >= maxIterations:
			print("failed to converge within maxIterations", file=self.outStream)
		
		return x
		
		
	
	def solve(self, tolerance=1e-7):
		self.buildMatrix()
		self.buildSource()
		
		import time
		start_time = time.time()
		
		from math import sin, cos, pi
		widthHeightRatio = sum(self.width)/len(self.width)*len(self.height)/sum(self.height)
		#Numerical Recipes: The Art of Scientific Computing, chapter 20.5
		rho_Jacobi = (cos(pi/(self.M+1)) + widthHeightRatio**2*cos(pi/(self.N+1))) / (1 + widthHeightRatio**2)
		omega_optimal = 2/(1 + (1-rho_Jacobi**2)**0.5)
		print('M, N:', self.M, self.N, file=self.outStream)
		print('rho_Jacobi:', rho_Jacobi, file=self.outStream)
		print('omega_optimal:', omega_optimal, file=self.outStream)
		self.fluxVector = self.fiveBandedSOR(omega=omega_optimal,tolerance=tolerance)
		print("%s seconds" % (time.time() - start_time), file=self.outStream)
		
		#turn 1d solution vector into 2d array
		self.flux = [self.fluxVector[i:i+self.N+1] for i in range(0, (self.M+1)*(self.N+1), self.N+1)]
		
		return self.flux
	
	def cellAvgFlux(self):
		return [[sum(self.flux[l][k] for k in [i, i+1] for l in [j, j+1])/4 for i in range(self.N)] for j in range(self.M)]
	
	def printInput(self, outputStream='default'):
		if outputStream == 'default':
			outputStream = self.outStream
		
		print('\ndiffusion solver inputs:\n', file=outputStream)
		print('rows:', self.M, file=outputStream)
		print('columns:', self.N, file=outputStream)
		print('\n', file=outputStream)
		
		print('cell widths', file=outputStream)
		print(', '.join(str(val) for val in self.width), file=outputStream)
		print('\n', file=outputStream)
		
		print('cell heights', file=outputStream)
		print('\n'.join(str(val) for val in reversed(self.height)), file=outputStream)
		print('\n', file=outputStream)
		
		print('boundary conditions', file=outputStream)
		print('left:', self.leftBoundary, file=outputStream)
		print('right:', self.rightBoundary, file=outputStream)
		print('top:', self.topBoundary, file=outputStream)
		print('bottom:', self.bottomBoundary, file=outputStream)
		print('\n', file=outputStream)
		
		#print all given physics values
		for key in self.mesh[0][0].__dict__:
			print(key, file=outputStream)
			for row in reversed(self.mesh): #reversed so that it prints such that y = 0 is on the bottom
				print(', '.join(str(val.__dict__[key]) for val in row), file=outputStream)
			
			print('\n\n', file=outputStream)
		
	
	def printFlux(self, outputStream='default'):
		if outputStream == 'default':
			outputStream = self.outStream
		
		print('flux', file=outputStream)
		for row in reversed(self.flux):
			print(', '.join(str(val) for val in row), file=outputStream)
		
		print('\n\n', file=outputStream)
	
	def printCellAvgFlux(self, outputStream='default'):
		if outputStream == 'default':
			outputStream = self.outStream
		
		print('cell averaged flux', file=outputStream)
		for row in reversed(self.cellAvgFlux()):
			print(', '.join(str(val) for val in row), file=outputStream)
		
		print('\n\n', file=outputStream)
	
	def printGeneralInfo(self, outputStream='default'):
		if outputStream == 'default':
			outputStream = self.outStream
		
		from datetime import datetime
		print("\n\ncurrent date and time:\n", datetime.now(), file=outputStream)
		
		print("""
NE 155 diffusion solver final project
version 1
Alexander Blank
Spring 2017

""", file=outputStream)
	