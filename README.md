# rotorNoise
Simple matlab code to estimate the loading noise of a generic propeller. The geometrical properties of the rotor are sorted by the Propeller class, the aerodynamical properties (Reynolds number, velocity vector ecc.)  for each section are calculated through simple analytical function in the TTcalculator class, wich returns also the total Thrust of the propeller. 
The local thrust is then passed to the Schlegel function, wich implements the integral formula defined by [Schlegel](https://apps.dtic.mil/sti/citations/AD0645884).
returns the noise of the first 10 harmonics.
The local cl and cd are calculated with [Xfoil](https://it.mathworks.com/matlabcentral/fileexchange/49706-xfoil-interface-updated), thanks to the matlab interface.


# References
Louis Edelman (2022). Xfoil Interface Updated (https://www.mathworks.com/matlabcentral/fileexchange/49706-xfoil-interface-updated), MATLAB Central File Exchange. Retrieved September 13, 2022.

Helicopter Rotor Noise Generation And Propagation, 1966-10-01, R. Schlegel, R. King and H. Mull
