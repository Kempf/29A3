# ENGN2229 Assignment 3
ENGN2229 Assignment 3 by Paul Apelt and Stephen Lonergan.

`soln( a, n, l, dt, p0, w0 )` - Phase map and trajectory plotting:

*   a - angle of attack (rad)
*   n - resolution of the phase map
*   l - number of solutions to compute
*   dt - time step of the solutions
*   p0 - starting roll angle (rad)
*   w0 - starting angular velocity (rad/tick)

`pmap( a, tol, dt, dw )` - Poincare map:
*   a - angle of attack (rad)
*   tol - omega tolerance
*   dt - time step of the solutions
*   dw - resolution of the p-map

See `traj.txt` for notable trajectories and points of interest.
