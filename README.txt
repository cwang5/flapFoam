# flapFoam
Initially written for OpenFOAM 1.6.x to simulate transient passive-flap motion. It can be compiled in 1.7.x but there appears to be some problem with the setup that would crash the simulation. It will not compile in OpenFOAM 2.0.0 and up due to changes in code structure. 
The equation of motion is based on beta_move = 0.5 * Moment_Inertia * realDeltaT * realDeltaT + beta_dot_0 * realDeltaT. 
I'm leaving the code on github so I can come back later and see how I dealt with OpenFOAM coding; you will be better served to use FSI code from the OF-dev branch, which is written by people more familiar with the FSI part of things than I was.
