import numpy as np
import scipy.linalg

def matrixPrint(name, matrix):
  print(name, ":", sep = '')
  for line in matrix:
    print('  '.join(map(str, line)))
  print("\n", end = '')

pArray = [1.0, 0.1, 0.01, 0.0001, 0.000001]

for p in pArray:
  print("-------------------------------------------------------------------------")
  # Initializing A and b
  A = np.array([[p+27, -6, -1, -6, -3, -4, -3, -4], [-6, 35, -1, -6, -5, -6, -3, -8],
  [-1, -1, 19, -6, -8, -2, 0, -1], [-6, -6, -6, 36, -4, -3, -4, -7], [-3, -5, -8, -4, 25, 0, -1, -4],
  [-4, -6, -2, -3, 0, 28, -8, -5], [-3, -3, 0, -4, -1, -8, 21, -2], [-4, -8, -1, -7, -4, -5, -2, 31]])
  
  b = np.array([[8*p+140], [-91], [-7], [142], [7], [-99], [25], [-117]])
  
  LU, p_LU = scipy.linalg.lu_factor(A) # DECOMP to A
  x2 = scipy.linalg.lu_solve((LU, p_LU), b) # solving the system to get x2 value

  E = np.array([[1, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0],
  [0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0],
  [0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 1]])
  Ainv = np.array([[0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0],
  [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0],
  [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]])
  Ainv = scipy.linalg.lu_solve((LU, p_LU), E)
  x1 = Ainv.dot(b) # A^(-1) * b = x1
  
  condN = np.linalg.cond(A)
  delta = np.linalg.norm(x1 - x2) / np.linalg.norm(x1)

  matrixPrint("x1", x1)
  matrixPrint("x2", x2)
  print("Parameter P: %5.6f" % (p))
  print("Condition number: ", condN)
  print("δ = ||x1 - x2|| / ||x1||: ", delta)
