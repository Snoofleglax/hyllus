import numpy as np
import scipy as sp
from scipy import optimize

args = (2, 3, 7, 8, 9, 10)
def f(x, *args):
  u, v = x
  a, b, c, d, e, f = args
  return a * u ** 2 + b * u * v + c * v ** 2 + d * u + e * v + f
def gradf(x, *args):
  u, v = x
  a, b, c, d, e, f = args
  gu = 2 * a * u + b * v + d
  gv = b * u + 2 * c * v + e
  return np.array((gu, gv))
x0 = np.asarray((0, 0))

res1 = optimize.fmin_cg(f, x0, fprime=gradf, args=args)

print(res1)

