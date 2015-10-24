# ENGN2229 Assignment 3
ENGN2229 Assignment 3 by Paul Apelt and Stephen Lonergan. Wing rock ODE analysis with lots of fancy graphs. Pls no plagiarise or we are kill (it's cool if you're not from ANU and not going to submit it as an assignment).

###Usage

See `traj.txt` for notable trajectories, usage examples and points of interest. Seriously, just do it. It won't run properly with random parameters.

###Functions:

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

`lcyc( a0, af, da, tol, dt, dw )` - Limit cycle finder:
*   a0, af - angle of attack range
*   da - angle of attack step
*   tol - omega tolerance
*   dt - time step of the solutions
*   dw - resolution of the p-map
