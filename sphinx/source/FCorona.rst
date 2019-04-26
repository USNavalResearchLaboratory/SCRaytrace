===============================
Physics for F Corona raytracing
===============================

PhysicsVSF (physics 4)
======================

Use the Lamy-Perrin VSF (1986), without any correction of the dust density based on the observer distance to the sun.

PhysicsVSFVaryDist (physics 5)
==============================

Use the Lamy-Perrin VSF (1986), with correction of the dust density based on the observer distance to the sun.

Validity is within 0.3AU to 1AU.


Line of sight integration for PhysicsVSFVaryDist, in equation
=============================================================
Constant factor from getConstFactors

.. math::

   btf &= ds \times RSUN_{CM} \times RSUN_{RSUN}^2 \times 1361 \\
   ds  &= (LosRangeStart - LosRangeEnd) / nbPoints \\
   Density &= C \frac{1}{r_p} \exp{A} \\
   B(r_{Obs}, r_p) &= Density \biggl( \frac{r_{Obs(in Rsun)}}{1AU(in Rsun)} \biggr) ^{-0.31} Vsf(\theta) \frac{1}{r_p^2}


Basic testing
=============
Using constant density of 1 (model 57), no VSF, no variation with radial distance, integration constant = ds
Integration result is equal to the depth of the integration; example: losrange=[0. , 215.] -> result = 215. This makes sense, no problem.


Calibration with VSF
--------------------
VSF units: [cm^-1.sr^-1]


Various F Corona dust distribution models
=========================================
* Various Fan model from Leinert, seen on page 275 of Lamy, Perrin 1986
* Lamy, Perrin 1986 model
* Dirbe/Siipola model


To Do
=====
.. todo::
   * Adjust the VSF on the KL model
