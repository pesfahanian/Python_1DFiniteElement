import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la


# RHS of BVP in form -u" + u = ?
def fn ( x ):
  value = 1 + x
  return value


def FEM ( ):

  # Problem domain interval and the number of nodes
  a = 0.0
  b = 1.0
  n = 4
  x = np.linspace ( a, b, n + 1 )

  print ( '' )
  print ( '  Nodes:' )
  for i in range ( 0, n + 1 ):
    print ( '  %d  %f' %( i, x[i] ) )
  
  # Seting a 3 point quadrature rule on the interval [0,1]
  ng = 3
  xg = np.array ( ( 0.11, 0.5, 0.88 ) )
  wg = np.array ( ( 0.25, 0.44, 0.25 ) )
  
  # Initializing the system matrix and the rhs
  A = np.zeros ( ( n + 1, n + 1 ) )
  rhs = np.zeros ( n + 1 )

  # Looking at elements 0 to n-1
  for e in range ( 0, n ):
    l = e
    r = e + 1
    x_left = x[l]
    x_right = x[r]

    # Considering quadrature point 0, 1, 2 in element E
    for i in range ( 0, ng ):

      # Transition mapping for simplicity
      x_i = x_left + xg[i] * ( x_right - x_left )
      w_i = wg[i] * ( x_right - x_left )

      # The basis functions and their derivatives
      psi_left = ( x_right - x_i  ) / ( x_right - x_left )
      d_psi_left = - 1.0 / ( x_right - x_left )=
      psi_right = ( x_i - x_left ) / ( x_right - x_left )
      d_psi_right = 1.0 / ( x_right - x_left )

      # Compute the system matrix and the rhs
      A[l][l] = A[l][l] + w_i * ( d_psi_left * d_psi_left + psi_left * psi_left ) 
      A[l][r] = A[l][r] + w_i * ( d_psi_left * d_psi_right + psi_left * psi_right )
      rhs[l]  = rhs[l]  + w_i *   psi_left * fn ( x_i )
      A[r][l] = A[r][l] + w_i * ( d_psi_right * d_psi_left + psi_right * psi_left )
      A[r][r] = A[r][r] + w_i * ( d_psi_right * d_psi_right + psi_right * psi_right )
      rhs[r]  = rhs[r]  + w_i *   psi_right * fn ( x_i )

  # Fitting the boundry conditions into the linear system
  A[0,0] = 1.0
  A[0,1:n+1] = 0.0
  rhs[0] = 1.0
  A[n,n] = 1.0
  A[n,0:n] = 0.0
  rhs[n] = 1.0

  # Solving the linear system
  u = la.solve ( A, rhs )

  # The exact solution at the nodes
  u_ex = [1, 1.0350, 1.0565, 1.0502, 1]

  # For each node, display Un and Uexact and calculate error
  print ( '' )
  print ( '  Node              Un          Uexact           Error' )
  for i in range ( 0, n + 1 ):
    if(i==0):
       err = 0
    else:
      err = abs ( u_ex[i] - u[i] )
    print ( '  %4d  %14.6g  %14.6g  %14.6g' % ( i, u[i], u_ex[i], err ) )

  # Plotting Un for each node
  plt.plot ( x, u, 'bo-', x, u_ex ,'r.')
  plt.show ( )

  return


FEM()