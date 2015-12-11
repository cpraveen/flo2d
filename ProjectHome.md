flo2d is a 2-dimensional flow and adjoint solver for inviscid and viscous fluids. It solves the Euler and Navier-Stokes equations on unstructured triangular grids using a vertex-centroid finite volume scheme. The adjoint solver is developed using the automatic differentiation tool TAPENADE. At present, only the flow solver is available, which is furthermore only for laminar flows.

flo2d is a further development of euler2d, which is a 2-D inviscid flow and adjoint solver. See the euler2d [website](http://euler2d.sourceforge.net)

### Some features ###

  * Inviscid and viscous (laminar only at present)
  * Finite volume method
  * Cell-centered scheme
  * Triangular grids (only)
  * Second order scheme (MUSCL-type)
  * Roe, KFVS
  * Explicit TVD-RK, LUSGS, GMRES (full newton using AD)