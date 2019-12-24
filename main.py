import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.colors import hsv_to_rgb
import numpy as np
from scipy import special
from collections import Counter
import math
import re
import time

class ANSIColor:
	"""Borrowed from Blender"""
	HEADER = '\033[95m'
	OKBLUE = '\033[94m'
	OKGREEN = '\033[92m'
	WARNING = '\033[93m'
	FAIL = '\033[91m'
	ENDC = '\033[0m'

class MyPolynomial:
	"""Custom polynomial container"""
	def __init__(self, coeff):
		self.degree = len(coeff) - 1
		self.coeff = coeff
		self.poly = np.poly1d(coeff)
	def __str__(self):
		out = ""
		i = self.degree
		for a_i in self.coeff:
			if (a_i == 0):
				i = i - 1
				continue
			if (i < self.degree):
				out = out + " + "
			out = out + str(a_i)
			if (i > 0):
				out = out + "x"
			if (i > 1):
				out = out + "^" + str(i)
			i = i - 1
		return out
	def eval(self, z, degree=0):
		if (degree == 0):
			return np.polyval(self.poly, z)
		return np.polyval(np.polyder(self.poly, degree), z)

quietMode = True		# Suppress most print statements
rectWidth = 2			# Create polygraph of region (-rectWidth - rectWidth*j) to (rectWidth + rectWidth*j)
maxIterations = 2**12	# Number of iterations to compute until declaring non-convergence to a root
deltaError = 1E-10		# Declare a root found when two iterations differ by this amount or less
criticalError = 1E-8	# Consider a point a critical point when it is this close to a critical point
resolution = 200			# Number of colored pixels in the output polynomiograph minus one and divided by two (i.e. an easily divisable even number like 4, 80, or 200)

def main(method):
	size = resolution*2 + 1
	x_min, x_max = -rectWidth, rectWidth
	y_min, y_max = -rectWidth, rectWidth
	dx = (x_max - x_min) / (size - 1)
	dy = (y_max - y_min) / (size - 1)
	x = np.linspace(x_min, x_max, size).round(6)
	y = np.linspace(y_min, y_max, size).round(6)
	A = []
	roots = []
	total_rounds = 0
	print("Beginning root computation...0%")
	progress = 0
	execTimeStart = time.time()
	for i in x:
		for u in range(5, 100, 5):
			if (i >= (x_min + (x_max - x_min)*(u / 100)) and progress < u):
				progress = u
				print(str(u) + "%")
		for j in y:
			root, rounds = findRoot(method, P, np.complex(i, j))
			total_rounds = total_rounds + rounds
			if np.isnan(root):
				root = np.nan
			if root not in roots:
				roots.append(root)
			A.append((i, j, roots.index(root)))		# change to index of root in root list to give unique value on colormap
	execTimeDelta = time.time() - execTimeStart
	A = np.array(A).T
	x, y, z = [item.flatten() for item in A]
	idx = np.round((x - x.min()) / dx).astype(np.int)
	idy = np.round((y - y.min()) / dy).astype(np.int)
	grid = np.zeros((size, size))
	grid[idx, idy] = z
	plt.imshow(np.flipud(grid.T), extent=(x.min(), x.max(), x.min(), x.max()), cmap='viridis', vmin=0, vmax=len(roots)-1, interpolation=None)
	plt.axvline(x=0,color='k', lw=0.2)
	plt.axhline(y=0,color='k', lw=0.2)
	rootI  = []
	for i in roots:
		rootI.append(roots.index(i))
	plt.scatter(np.real(roots), np.imag(roots), marker='.', s=80, linewidths=0.5, edgecolors='k', c=rootI, cmap='viridis', vmin=0, vmax=len(roots)-1)
	print("Roots: " + str(roots))
	print("Total iterations: " + str(total_rounds))
	print("Total time: " + str(execTimeDelta) + "s")
	plt.show()
	return None

def findRoot(method, poly, z0):
	iter = 1
	z1 = z0
	z2 = method(z0, poly)
	while (not np.isnan(z2) and abs(z2-z1) > deltaError):
		if (iter > maxIterations):
			if (method == hybridMethod):
				return findRoot(robustMethod, poly, z0)
			if (not quietMode): print(ANSIColor.WARNING + "Point " + str(z0) + " could not converge to a root." + ANSIColor.ENDC)
			return np.nan, iter
		iter = iter + 1
		z1 = z2
		z2 = method(z1, poly)
	#print(ANSIColor.OKGREEN + "Converged at root: " + str(z2) + " in " + str(iter) + " iteration(s)." + ANSIColor.ENDC)
	z2 = complex(round(np.real(z2), 6), round(np.imag(z2), 6))
	#print("Root from " + str(z0) + " is approximately: " + ANSIColor.OKBLUE + "{:.00000001f}".format(z2) + ANSIColor.ENDC)
	return z2, iter

def getInput():
	a = input(ANSIColor.HEADER + "Enter a polynomial's integer coefficients, from greatest degree to least degree, including zeroes:" + ANSIColor.ENDC + "\n")
	coeff = re.split(r" |,|;|\|", a.strip())
	coeff = [int(i) for i in coeff]
	b = input(ANSIColor.HEADER + "Enter the root-finding method you would like to use.\n(1) Newton, (2) Halley, (3) Robust Newton, (4) Hybrid Newton" + ANSIColor.ENDC + "\n")
	b = int(b)
	if (b == 1):
		b = newtonMethod
	elif (b == 2):
		b = halleyMethod
	elif (b == 3):
		b = robustMethod
	elif (b == 4):
		b = hybridMethod
	else:
		raise Exception(ANSIColor.FAIL + "Invalid method." + ANSIColor.ENDC + "\n")
	return MyPolynomial(coeff),b

def newtonMethod(z, p):
	p_1 = p.eval(z, 1)
	if (abs(p_1) <= criticalError or np.isnan(z)):
		return np.nan
	return z - p.eval(z)/p_1

def halleyMethod(z, p):
	p_1 = p.eval(z, 1)
	denominator = p_1*p_1 - p.eval(z)*(p.eval(z, 2)/2)
	if (abs(p_1) <= criticalError or np.isnan(z)):
		return np.nan
	return z - p.eval(z)*p_1/denominator

def robustMethod(z, p):
	p_0 = p.eval(z)
	if (p_0 == 0):			# Assert P(z0) != 0
		return z
	k = 1
	p_k = p.eval(z, k)
	while (abs(p_k) < criticalError):			# Find smallest degree of derivation 'k' for which d^k/dx^k[P(z0)] is not 0
		p_k = p.eval(z, k)
		k = k + 1
	u_k = 1.0/(special.factorial(k)) * p_0 * np.conjugate(p_k)
	y = 2.0 * np.real(u_k**(k-1))
	d = -2.0 * np.imag(u_k**(k-1))
	c_k = max(abs(y),abs(d))
	theta = None
	if ((c_k == abs(y)) and (y < 0)):
		theta = 0.0
	if ((c_k == abs(y)) and (y > 0)):
		 theta = np.pi/k
	if ((c_k == abs(d)) and (d < 0)):
		theta = np.pi/(2.0*k)
	if ((c_k == abs(d)) and (d > 0)):
		theta = 3.0*np.pi/(2.0*k)
	A = abs(p_0)
	for j in range(1, p.degree + 1):
		A_j = abs(p.eval(z, j))/(special.factorial(j))
		A = max(A, A_j)
	C_k = c_k * (abs(u_k))**(2-k)/(6.0*A*A)
	return z + C_k/3.0 * u_k/abs(u_k) * complex(math.cos(theta), math.sin(theta))

def hybridMethod(z, p):
	if (np.isnan(z)):
		return np.nan
	p_1 = p.eval(z, 1)
	if (abs(p_1) <= criticalError):
		return robustMethod(z, p)
	result = z - p.eval(z)/p_1
	if (np.isnan(result)):
		return robustMethod(z, p)
	return result

if __name__ == "__main__":
	P, method = getInput()
	print("Polynomial input:\n"+ str(P))
	main(method)