## ParkerSolarWind

This repository complete and extend the one with the exact same name that can be found at [https://github.com/STBadman/ParkerSolarWind](https://github.com/STBadman/ParkerSolarWind). It contains Python code which solves two-fluids hydrodynamic solar wind equations, for a 1D radial trans-sonic flow in spherical expansion, including super-radial expansion, based on Eugene Parker's theory of the solar wind ([Parker 1958](https://ui.adsabs.harvard.edu/abs/1958ApJ...128..664P/abstract), [Parker 1960](https://ui.adsabs.harvard.edu/abs/1960ApJ...132..821P/abstract)) and expansion factor modeling ([Kopp & Holzer 1976](https://ui.adsabs.harvard.edu/abs/1976SoPh...49...43K/abstract)).

This code follows the recent "isopoly" resolution by [Dakeyo et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...940..130D/abstract) and [Dakeyo et al. (2024b)](https://ui.adsabs.harvard.edu/abs/2022ApJ...940..130D/abstract), which model an isothermal evolution close to the Sun, followed by a polytropic evolution above the critical radius, including expansion factor influence in the near Sun region.

There are three main files here `function_iso_poly_dakeyo2024b`, `exe_function_iso_poly_dakeyo2024b` and `inputs_function_iso_poly_dakeyo2024b` . The former contains the functions to solve the equations themselves, a file that runs the solving and plots the solution, and a later that only takes the inputs of the model then displays the solutions and provides the outputs solar parameters (u, Tp, Te, n, f). The equations are solved by a finite difference scheme. 
There are three types of solutions that can be obtained :

* A fully isothermal solar wind ($\gamma_p =  \gamma_e = 1$) <details><p> - This follows [Parker 1958](https://ui.adsabs.harvard.edu/abs/1958ApJ...128..664P/abstract), in which the solar wind fluid is held at a fixed temperature. Mass flux conservation results in a negative density gradient and in turn an outwards directed pressure gradient force. For sufficiently hot $T_0$, this outwards force outcompetes gravitation, resulting in a trans-sonic solar wind flow out to infinity. While such a constant temperature is non-physical in the heliosphere, it is a reasonable first approximation to behavior in the solar corona where coronal heating operates.</p></details>
* A two fluid isopoly solar wind with single transition ($\gamma_p$ and $\gamma_e \neq 1$  with  $r_{iso|p} = r_{iso|e}$) <details><p> - Here, the solar wind temperature is allowed to cool with heliocentric distance, as is observed to actually occur in the solar wind (e.g. [Dakeyo et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022ApJ...940..130D/abstract). 
*This consists of an initial isothermal evolution (isothermal layer) out to some boundary distance called the "isothermal radius" $r_{iso}$, which can nominally be interpreted as defining a corona as the region in which coronal heating (as an abstract physical process) operates. In this solution case, both protons and electrons share the transition, i.e. $r_{iso} = r_{iso|p} = r_{iso|e}$. 
For $r \gt r_{iso}$, the solar wind is constrained to follow a polytropic evolution which is initialized by the outer boundary conditions of the isothermal region. Protons and electrons can follow differentiate polytropic evolution ($gamma_p \neq \gamma_e$ is possible). 
 For most combinations of physical conditions, the trans-sonic critical point is located within the isothermal region. As long as the isothermal boundary is sufficiently high that the solar wind stays super-sonic at the (unphysical discontinuity) transition to polytropic behavior, the solution remains on the asymptotically accelerating solution branch and a reasonable solar wind solution is obtained. The transition can be smoothed (no discontinuity anymore) considering slowy varying polytropic indixes at the transition between the two regions, but this feature is not adressed here and may required deepest work. </details></p> 
* A two fluid isopoly solar wind with double transition ($\gamma_p$ and $\gamma_e \neq 1$  with  $r_{iso|p} \neq r_{iso|e}$) <details><p> - This case is closely similar to the single transition, at the difference that protons and electrons do not share the same isothermal radius. </details></p>

Each of these solution functions returns an array of heliocentric distances, mass densities, fluid velocities and fluid temperatures, as well as any parameters that went into the solution.

Units are tracked with `astropy.units` and the outputs of the above functions are `astropy.units.Quantity` objects.

`plot_parkersolarwind` subsequently contains plotting functions which are expecting these same output arrays and parameters. 

In the following example, we solve and plot an isothermal layer solution :

```python
import parkersolarwind as psw
import numpy as np
import astropy.units as u
import astropy.constants as const
import matplotlib.pyplot as plt
sol = psw.solve_isothermal_layer(
    np.logspace(0,np.log10(200),400)*u.R_sun, # Radial grid points
    11.5*u.R_sun, # Isothermal layer outer radius
    1.64*u.MK, # Isothermal layer temperature
    1.38, # Polytropic layer polytropic index
    n0=11e6/u.cm**3 # Density normalization
) 
fig, axes = psw.plot_isothermal_layer(sol)
plt.show()
```
![image](IsoLayerExample.png)

Examples for the other two types of solution as a function of varying input parameters can be seen in `ExampleNotebook.ipynb`

As of (2/6/2023) : 

* The average corpuscular mass defaults to $\mu=0.5$ but may be tweaked as a parameter to the solve functions
* (Future enhancement) Radial flux tube expansion is assumed throughout all space 
* (Known Issue) For `solve_parker_polytropic` No check is currently done that the input parameters lie in the allowed region of $(T_\odot-\gamma)$ space ([Shi et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022PhPl...29l2901S/abstract)), so it is possible to generate an all-NaN solution for this reason.
* (Known Issue) For `solve_parker_polytropic`, there are some numerical instability issues which can result in solutions in the wrong branch or all-NaN solutions. This is likely due to the implementation of `scipy.optimize.root` with which we currently integrate from $1R_\odot$ outwards and use the previous solution as the initial guess, except at the critical point where we enforce the guess to jump to above the critical speed. A future implementation should integrate explicitly from the critical point down and up separately with adaptive radial grid spacing and the local solution gradient should be used to improve the guess.

Update (7/3/2023) :

Following [Shi et al. (2022)](https://ui.adsabs.harvard.edu/abs/2022PhPl...29l2901S/abstract) Section D we implement the ability to include an external force in the solutions. We do this for each of the isothermal, polytropic and isothermal layer solutions. Plotting options are included to display the force parameters where the figure of merit is the force per unit mass relative to $GM_\odot/R^2$.

A new set of example plots is appended to [ExampleNotebook.ipynb](https://github.com/STBadman/ParkerSolarWind/blob/main/ExampleNotebook.ipynb) to illustrate these updates. 
