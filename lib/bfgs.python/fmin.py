from scipy.optimize import fmin_bfgs

def func1(x):
	y = x ** 2 + x
	return y

print fmin_bfgs(func1, 0)
