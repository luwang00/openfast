
.. _AD_theory:

AeroDyn Theory
==============

This theory manual is work in progress, please refer to the AeroDyn 14 manual for more details :cite:`ad-AeroDyn:manual`. Many changes have occurred since AeroDyn 14 (e.g. BEM formulation, coordinate system used in the BEM equations, dynamic stall, dynamic BEM), but these changes are not yet documented here.



Steady BEM
~~~~~~~~~~

The steady blade element momentum (BEM) equations are solved as a constrained equation, and the formulation follows the description from Ning :cite:`ad-Ning:2014`.



.. _AD_DBEMT:

Dynamic BEM Theory (DBEMT)
~~~~~~~~~~~~~~~~~~~~~~~~~~



Two equivalent versions of Oye's dynamic inflow model are implemented in AeroDyn.
The first one uses discrete time, it can be used with the constant-tau1 model 
(``DBEMT_Mod=1``) or the varying-tau1 model (``DBEMT_Mod=2``), but it cannot be used for linearization.
The second version uses a continuous-time state-space formulation  (``DBEMT_Mod=1``), it assumes a constant-tau1, and can be used for linearization.  
For a same value of :math:`\tau_1`, the discrete-time and continuous-time formulations returns exactly the same results.





Oye's dynamic inflow model consists of two first-order differential equations (see :cite:`ad-Branlard:book`):

.. math::
   \begin{aligned}
       \boldsymbol{W}_\text{int}+\tau_1    \boldsymbol{\dot{W}}_\text{int} &= \boldsymbol{W}_\text{qs} + k \tau_1 \boldsymbol{\dot{W}}_\text{qs} \\
       \boldsymbol{W}+\tau_2 \boldsymbol{\dot{W}} &= \boldsymbol{W}_\text{int}
   \end{aligned}

where 
:math:`\boldsymbol{W}` is the dynamic induction vector at the rotor (at a given blade position and radial position), 
:math:`\boldsymbol{W}_\text{qs}` is the quasi-steady induction, 
:math:`\boldsymbol{W}_\text{int}` is an intermediate value coupling the quasi-steady and the actual inductions (may be discontinuous if the quasi-steady indution is discontinuous).
and
:math:`(\dot{\ })` represents the time derivative.
The coupling constant :math:`k`, with values between 0 and 1, is usually chosen as :math:`k=0.6`.
Oye's dynamic inflow model relies on two time constants, :math:`\tau_1` and :math:`\tau_2` :

.. math::
        \tau_1=\frac{1.1}{1-1.3 \min(\overline{a},0.5)} \frac{R}{\overline{U}_0}
        , \qquad
        \tau_2 =\left[ 0.39-0.26\left(\frac{r}{R}\right)^2\right] \tau_1

where :math:`R` is the rotor radius, :math:`\overline{U}_0` is the average wind speed over the rotor, :math:`\overline{a}` is the average axial induction over the rotor, and :math:`r` is the radial position along the blade.
For ``DBEMT_Mod=1`` or ``DBEMT_Mod=3``, the user needs to provide the value of :math:`\tau_1`.




The continuous-time state-space formulation of the dynamic inflow model (``DBEMT_Mod=3``) was derived in :cite:`ad-Branlard:2022`.

.. math::
   \begin{align}
      \begin{bmatrix}
      \boldsymbol{\dot{W}}_\text{red}\\
      \boldsymbol{\dot{W}}\\
      \end{bmatrix}
      =
      \begin{bmatrix}
      -\frac{1}{\tau_1}\boldsymbol{I}_2 & \boldsymbol{0} \\
       \frac{1}{\tau_2}\boldsymbol{I}_2 &
      -\frac{1}{\tau_2}\boldsymbol{I}_2 \\
      \end{bmatrix}
      \begin{bmatrix}
      \boldsymbol{W}_\text{red}\\
      \boldsymbol{W}\\
      \end{bmatrix}
      +
      \begin{bmatrix}
       \frac{1-k}{\tau_1} \\
       \frac{k}{\tau_2}\\
      \end{bmatrix}
     \boldsymbol{W}_\text{qs}
   \end{align}

where 
:math:`\boldsymbol{I}_2` is the 2x2 identity matrix,
:math:`\boldsymbol{W}_\text{red}` is the reduced induction which is a continuous, scaled, and lagged version of the quasi-steady induction, defined as:

.. math::
    \boldsymbol{W}_\text{int} = \boldsymbol{W}_\text{red} + k \boldsymbol{W}_\text{qs} 


The discrete-time version of the model is documented in the unpublished manual of DBEMT.
The current discrete-time formulation is complex and in the future it can be simplified by using :math:`\boldsymbol{W}_\text{red}`.






.. _AD_twr_influence:

Tower influence models
~~~~~~~~~~~~~~~~~~~~~~~~~~~

AeroDyn can model the influence of the tower on the flow reaching the blades through two
superimposable effects: a potential-flow disturbance of the flow around the tower
(``TwrPotent``) and a downstream shadow (wake) velocity deficit (``TwrShadow``). Both are
evaluated at each blade node from the position of the node relative to the nearest point on
the tower, expressed in a local tower reference frame centered on the tower axis:
:math:`\overline{x}` is the downstream (streamwise) coordinate, :math:`\overline{y}` the
lateral coordinate, and :math:`\overline{z}` the coordinate along the tower axis, each
normalized by the local tower radius. The resulting disturbance-velocity fractions are scaled by
:math:`W_\text{tower}`, the local incoming wind speed component normal to the tower axis, and
applied in the plane normal to the tower. :math:`C_d` is the local tower drag coefficient and
:math:`TI` (Eames model only) is the local turbulence intensity at the tower node. It is
convenient to define :math:`\overline{r} = \sqrt{ \overline{x}^2 + \overline{y}^2 }`.

.. _AD_twr_potent:

Tower potential-flow model
-----------------------------

The baseline potential-flow model (**TwrPotent=1**) represents the tower cross section as a
two-dimensional cylinder in potential flow (a doublet), giving the streamwise
(:math:`u_{TwrPotent}`) and lateral (:math:`v_{TwrPotent}`) disturbance-velocity fractions:

.. math::
   u_{TwrPotent} = \frac{-\overline{x}^2 + \overline{y}^2}{\overline{r}^4},
   \qquad
   v_{TwrPotent} = \frac{-2\,\overline{x}\,\overline{y}}{\overline{r}^4}

The Bak correction (**TwrPotent=2**), following Bak, Madsen, and Johansen (2001), adds a
drag-induced source term and offsets the streamwise coordinate by :math:`+0.1` (a fixed
empirical constant) to better represent the near wake. Writing
:math:`\overline{x}' = \overline{x} + 0.1` and
:math:`\overline{r}' = \sqrt{ \overline{x}'^2 + \overline{y}^2 }`,

.. math::
   u_{TwrPotent} = \frac{-\overline{x}'^2 + \overline{y}^2}{\overline{r}'^4}
                 + \frac{C_d\,\overline{x}'}{2\pi\,\overline{r}'^2},
   \qquad
   v_{TwrPotent} = \frac{-2\,\overline{x}'\,\overline{y}}{\overline{r}'^4}
                 + \frac{C_d\,\overline{y}}{2\pi\,\overline{r}'^2}

Alongside the tower, the potential-flow disturbance velocity is applied in full.
Beyond the tower ends the influence is tapered and cut off; this end treatment,
which is shared with the shadow model, is described in :numref:`AD_twr_ends`.

.. _AD_twr_shadow:

Tower shadow models
----------------------

Powles tower shadow model (**TwrShadow=1**) is given by:

.. math::
   u_{TwrShadow} = - \frac{C_d}{  \sqrt{\overline{r}}  }
               \cos\left( \frac{\pi}{2}\,\frac{\overline{y}}{\sqrt{\overline{r}}}\right)^2

where :math:`\overline{r} = \sqrt{ \overline{x}^2 + \overline{y}^2 }`.


Eames tower shadow model (**TwrShadow=2**) is given by:

.. math::
   u_{TwrShadow} = -\frac{C_d}{ TI \: \overline{x} \, \sqrt{2 \pi }  }
               \exp{\left(  -\frac{1}{2}  \left(\frac{ \overline{y}}{ TI \: \overline{x} } \right)^2 \right) }

where :math:`TI` is the turbulence intensity at the tower node. 

To avoid excessive flow reversal behind the tower, the shadow
deficit fraction is limited to :math:`u_{TwrShadow} \ge -0.5`.

The potential-flow and shadow contributions are superimposed and scaled by
:math:`W_\text{tower}` to give the velocity perturbation at the blade node in the local tower
frame:

.. math::
   v_x = \left( u_{TwrPotent} + u_{TwrShadow} \right) W_\text{tower},
   \qquad
   v_y = v_{TwrPotent}\, W_\text{tower}

This perturbation :math:`(v_x, v_y, 0)` is rotated into the earth-fixed frame and added to the
free-stream (undisturbed) inflow velocity to obtain the disturbed inflow velocity used by the
blade-element calculations.

.. _AD_twr_ends:

Tower ends and exclusion zones
------------------------------

Both tower-influence effects use the position of the blade node relative to the nearest point on
the tower. Alongside the tower this nearest point is the orthogonal projection of the blade node
onto the tower axis, for which the axial coordinate is :math:`\overline{z} = 0`. When a blade
node lies axially beyond a tower end, the nearest point becomes the tower end itself, and
:math:`\overline{z}` is then the axial distance from that end, normalized by the local tower
radius.

For each blade node, the tower clearance is computed as :math:`c = \lVert \mathbf{r} \rVert - R`,
where :math:`\mathbf{r}` is the vector from the nearest tower point to the blade node and
:math:`R = \tfrac{1}{2}` ``TwrDiam`` is the local tower radius. The disturbance is suppressed both
very close to the tower (:math:`c \le 0.01\,` ``TwrDiam``) and far from it (:math:`c > 20\,` ``TwrDiam``,
a far-field cutoff). Because the clearance beyond a tower end is measured to the end point, the
surface of minimum clearance for the tower influence models there is a hemispherical cap closing off
the cylinder. The near-tower exclusion zone is therefore a capsule --- a cylinder capped by hemispheres
at both ends --- rather than a bare cylinder.

Beyond a tower end the disturbance is faded out over one tower radius using a cosine-squared axial
taper. For :math:`|\overline{z}| < 1` the in-plane coordinates are divided by
:math:`\cos\!\left(\tfrac{\pi}{2}\overline{z}\right)`,

.. math::
   \overline{x} \;\rightarrow\; \frac{\overline{x}}{\cos\!\left(\tfrac{\pi}{2}\overline{z}\right)},
   \qquad
   \overline{y} \;\rightarrow\; \frac{\overline{y}}{\cos\!\left(\tfrac{\pi}{2}\overline{z}\right)}

The taper is applied purely through this coordinate scaling; how it attenuates each effect
depends on how that effect depends on :math:`\overline{x}` and :math:`\overline{y}`. Because the
potential-flow disturbance velocity scales as
:math:`1/\overline{r}^{\,2} = 1/(\overline{x}^2 + \overline{y}^2)`, dividing :math:`\overline{x}`
and :math:`\overline{y}` in this way multiplies the potential-flow velocity by exactly
:math:`\cos^2\!\left(\tfrac{\pi}{2}\overline{z}\right)`: it decreases smoothly from its full value
at the end plane (:math:`\overline{z} = 0`) to zero at :math:`\overline{z} = 1` with zero slope
there (:math:`C^1`-continuous). The shadow deficit is tapered by the same coordinate scaling, but
its dependence on :math:`\overline{x}` and :math:`\overline{y}` is different, so the attenuation is
not a clean :math:`\cos^2` factor (for the Powles model the deficit scales roughly as
:math:`\sqrt{\cos}` with an additional shift in its lateral argument, and for the Eames model it
scales linearly in :math:`\cos` with the lateral profile unchanged). Both shadow models
nonetheless fade to zero as :math:`\overline{z} \rightarrow 1`. No tower influence is applied for
:math:`|\overline{z}| \ge 1`.

.. note::
   The tower shadow deficit is currently tapered near the tower ends only as a byproduct of the
   coordinate scaling above, which does not give a principled axial wake taper. The end
   treatment of the tower shadow is expected to be improved in a future release.


.. _AD_buoyancy:

Buoyancy
~~~~~~~~

When a solid object is submerged in a fluid, it experiences a net force, buoyancy, 
from the hydrostatic fluid pressure acting on its surface. This force can often 
be neglected in less dense fluids, such as air, but can be significant in denser
fluids, such as water. To capture the effects of this force on MHK turbines, 
buoyant loads are calculated for the turbine blades, tower, hub, and nacelle. 
Marine growth is neglected for all components. :numref:`AD_buoy_coords` - 
:numref:`AD_buoy_hubnacelle` detail the coordinate systems and blade, tower, hub,
and nacelle buoyancy calculations.

.. _AD_buoy_coords:

Coordinate Systems
------------------
The buoyant force acting on an element depends on its instantaneous orientation 
and depth. The orientation is defined by heading and inclination angles, which
are calculated for each element at every time step. Total water depth is defined
by the user, relative to the still water level (or relative to the mean sea 
level when running AeroDyn in standalone mode with the AeroDyn driver). The
instantaneous depth of each element is based on its position in global coordinates
at each time step.

.. _AD_buoy_bladestower:

Blades and Tower
----------------
To allow for an efficient analytical solution, the blades and tower are modeled
as tapered cylinders. The cross-sectional area of the tapered cylinders is set
equal to the blade or tower cross-sectional area. Loads are estimated by breaking
the blade or tower into elements of a given length and integrating the hydrostatic
pressure over the wetted area of each element. For the blades, loads are applied
at a user-specified center of buoyancy. For the tower, loads are applied at the
centerline. When applicable, end effects are accounted for by calculating the
fluid pressure on the exposed axial face of the element. The tower is assumed to
be either embedded into the seabed or attached to another support structure member,
such that no end effects at the tower base are needed. For MHK turbines with a support
structure (i.e., any structure other than a simple tower embedded in the seabed), it
is currently recommended to model the entire support structure, including the tower, in HydroDyn.
Future releases will include the ability to neglect fluid loads at the interface between a tower
modeled in AeroDyn and a platform modeled in HydroDyn.


The buoyancy calculation for the blades and tower is completed according to the following steps:

1.	Calculate parameters related to element geometry that do not change with time
2.	Check that no elements cross the free surface or go beneath the seabed
3.	Calculate the instantaneous orientation and depth of each element
4.	Integrate hydrostatic fluid pressure over the wetted surface of each element and express as a force acting at the center of buoyancy
5.	For blades, calculate the buoyant force on the axial face of the blade root and tip; add the tip force to the adjacent element and store the root force
6.	For the tower, calculate and store the buoyant force on the axial face of the tower top
7.	Move buoyant loads from the center of buoyancy to the aerodynamic center
8.	Express buoyant loads in the form expected by OpenFAST
9.	Add buoyant loads to aerodynamic loads

Although the blade and tower buoyant loads are not based on volume, the volumes of these components are
written to the AeroDyn summary file for reference. The blade and tower volumes are calculated by summing
the volume of each element, assumed to be a tapered cylinder. The volume of a single element :math:`(V_{elem})`
is given by:

.. math::
   V_{elem} = \frac{\pi}{3} (r_i^2 + r_i r_{i+1} + r_{i+1}^2) dl

where :math:`r_i` is the element radius at node :math:`i`, :math:`r_{i+1}` is the element radius at node 
:math:`i+1`, and :math:`dl` is the element length.


.. _AD_buoy_hubnacelle:

Hub and Nacelle
---------------
The hub and nacelle are treated as separate components. The buoyant force is 
determined by the volume of either the hub or nacelle and applied at its 
user-specified center of buoyancy. Corrections are made to account for the joints
between the hub and blades and the nacelle and tower, as the joint locations are
not exposed to fluid pressure. No correction is made for the joint between the 
hub and nacelle.

The buoyancy calculation for the hub and nacelle is completed according to the following steps:

1.	Check that the component does not cross the free surface or go beneath the seabed
2.	Calculate the instantaneous depth of the component
3.	Calculate the buoyant force from the volume of the component
4.	Move buoyant loads from the center of buoyancy to the aerodynamic center
5.	For the hub, correct loads to account for the joints with each blade
6.	For the nacelle, correct loads to account for the joint with the tower

.. _AD_addedmass_inertia:

Added Mass and Fluid Inertia
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Added mass loads are caused by body and fluid accelerations.
These forces can often be neglected in less dense fluids, such as air, but can be significant in denser
fluids, such as water. To capture the effects of these forces on MHK turbines, 
added mass and fluid inertia loads are calculated for the turbine blades and tower.
Per-unit-length loads are estimated at each blade or tower node by calculating the added mass and fluid inertia
forces according to the appropriate terms from Morison's equation. The resulting loads are summed with the
previously calculated hydrodynamic and/or buoyant per-unit-length loads. 
Loads for the blades are applied at the aerodynamic center. Loads for the tower are applied at the centerline.
Marine growth and end effects are neglected, and members are not allowed to cross the free surface
(i.e., members are always fully submerged). Ballast is not considered. Nodes do not need to be uniformly spaced,
and axial loads are neglected. The tower is assumed to be axisymmetric (with the same coefficients used in both transverse directions),
but the blade is not (with different coefficients normal and tangential to the chord, as well as an added mass coefficient for pitch).

.. _AD_addedmass_inertia_Morison:

Morison's Equation
------------------
Added mass and fluid inertia loads are calculated according to the appropriate terms from Morison's equation. The added mass force is given as

.. math::
   F_{a} = \rho C_a V (\dot{u} - \dot{v})

where :math:`\rho` is the fluid density, :math:`C_a` is the added mass coefficient, :math:`V` is the element volume, :math:`\dot{u}` is the 
fluid acceleration, and :math:`\dot{v}` is the body acceleration.

The fluid inertia force is given as

.. math::
   F_{i} = \rho C_p V \dot{u}

where :math:`C_p` is the dynamic pressure coefficient.

The fluid density and added mass and dynamic pressure coefficients are user-specified. Added mass and fluid
inertia loads can be turned off by setting the relevant coefficients to zero. Additional information about calculating added mass coefficients can be
found in :numref:`AD_user_guide` ("Determination of Added Mass Coefficients for Floating Hydrokinetic Turbine Blades using Computational Fluid Dynamics").
The body and fluid accelerations are calculated internally and passed to AeroDyn. Body accelerations are available from the structural solver (or driver),
and fluid accelerations are calculated based on the inflow velocity time series. Added mass and fluid inertia loads are calculated as per-unit-length within
AeroDyn. Therefore, :math:`V` is taken as the cross-sectional area at the node of interest. For the blades, the reference cross-sectional area for the normal
and tangential terms is chord*thickness (:math:`ct`). This is expressed as :math:`(c^2)(t/c)`, where :math:`t/c` (i.e., ``t_c``) is specified
in the AeroDyn blade input file and cannot be less than 0. For the tower, the reference cross-sectional area is :math:`\pi r^2` where :math:`r` 
is calculated as (0.5 ``TwrDiam``). The normalization for the ``BlCpn``, ``BlCpt``, ``BlCan``, and ``BlCat`` coefficients should be :math:`\rho ct`;
the normalization for the ``BlCam`` coefficient should be :math:`(1/12)\rho ct(c^2+t^2)`; and the normalization for the ``TwrCp`` and ``TwrCa`` coefficients should
be :math:`\rho\pi(0.5` ``TwrDiam``) :math:`^2`.

Blade Added Mass and Fluid Inertia
----------------------------------
Added mass and fluid inertia loads are calculated for the normal-to-chord, tangential-to-chord, and pitch directions in the blade coordinate system.
The following coefficients are defined by the user in the AeroDyn blade input file:

-  ``BlCpn`` specifies the blade normal-to-chord dynamic pressure coefficient;
   to neglect normal-to-chord fluid inertia loads on the blade, set ``BlCpn`` to 0

-  ``BlCpt`` specifies the blade tangential-to-chord dynamic pressure coefficient;
   to neglect tangential-to-chord fluid inertia loads on the blade, set ``BlCpt`` to 0

-  ``BlCan`` specifies the blade normal-to-chord added mass coefficient, cannot be less than 0;
   to neglect normal-to-chord added mass loads on the blade, set ``BlCan`` to 0

-  ``BlCat`` specifies the blade tangential-to-chord added mass coefficient, cannot be less than 0;
   to neglect tangential-to-chord added mass loads on the blade, set ``BlCat`` to 0

-  ``BlCam`` specifies the blade pitch added mass coefficient, cannot be less than 0;
   to neglect pitch added mass loads on the blade, set ``BlCam`` to 0

Tower Added Mass and Fluid Inertia
----------------------------------
Added mass and fluid inertia loads are calculated for the transverse direction in the tower coordinate system.
The following coefficients are defined by the user in the AeroDyn primary input file:

-  ``TwrCp`` specifies the tower transverse dynamic pressure coefficient;
   to neglect fluid inertia loads on the tower, set ``TwrCp`` to 0

-  ``TwrCa`` specifies the tower transverse added mass coefficient, cannot be less than 0;
   to neglect added mass loads on the tower, set ``TwrCa`` to 0

