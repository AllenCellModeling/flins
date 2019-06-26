===============
Actin diffusion
===============

On first glance at a model visualization, it appears that actin is trucking
along at a prodigious rate when unbound. I was initially concerned that I had
overestimated the diffusive constants for actin filaments. Let's check to see
that they are within reason.

We find the diffusive coefficient for actin by treating it as a cylinder
diffusing about its long axis. From page 107 of `Howard, 2001`_ we get that
this will produce a drag coefficient of

.. math::  \gamma = \frac{2 \pi \eta L}{\ln(L/2r) - 0.20}

We start with :math:`\eta` as the viscosity of water, but multiply it by a
factor of 3.2 to account for the `increased crowding of cytoplasm`_.  :math:`L`
will vary from filament to filament, but let's take a short example length of
14 g-actin pairs or approximately :math:`40 nm`.

From this we can calculate the diffusion constant of our sample actin as

.. math::  D = \frac{k T}{ \gamma } 

Now, I'm always a bit paranoid about dimension checking, so let's explicitly
calculate the diffusion constant taking units into account. The `misu`_ package
will help here. Here is a summary of our units

+----------------+--------------------+-----------------+-------------+
| Variable       | Value              | Units           | Source      |
+================+====================+=================+=============+
| :math:`k`      | :math:`1.38E-23`   | :math:`J/K`     | `1`_        |
+----------------+--------------------+-----------------+-------------+
| :math:`T`      | :math:`277`        | :math:`K`       | choice      |
+----------------+--------------------+-----------------+-------------+
| :math:`\eta`   | :math:`0.00365`    | :math:`Pa\,s`   | `2`_, `3`_  |
+----------------+--------------------+-----------------+-------------+
| :math:`L`      | :math:`40`         | :math:`nm`      | choice      |
+----------------+--------------------+-----------------+-------------+
| :math:`r`      | :math:`3`          | :math:`nm`      | `4`_        |
+----------------+--------------------+-----------------+-------------+

.. _Howard, 2001: http://books.google.com/books?vid=ISBN9780878933334
.. _increased crowding of cytoplasm: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=102561
.. _misu: https://github.com/cjrh/misu
.. _1: https://www.wolframalpha.com/input/?i=boltzmann+constant
.. _2: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=102561
.. _3: https://www.wolframalpha.com/input/?i=3.2+*+water+viscosity+at+288K
.. _4: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=106827

.. code:: ipython3

    import misu as m
    import numpy as np
    
    k = 1.38E-23 * m.J * m.K**-1
    T = 288 * m.K
    eta = 0.00356 * m.Pa_s  
    L = 40 * m.nanometre
    r = 3 * m.nanometre
    
    D = (k * T)/((2 * np.pi * eta * L)/(np.log(L/(2*r)) - 0.20))
    print("D is %0.1f microns**2/sec"%(D >> (m.micron**2/m.s)))


.. parsed-literal::

    D is 7.5 microns**2/sec


That compares very reasonably with an `observed diffusion constant`_ of
:math:`5.2 \mu m ^2/s` for :math:`40 nm` long actin fragments in yeast cytosol.
We can be reassured that, while our actin diffusion appears to be quite quick,
it matches well to observed values.

.. _observed diffusion constant: https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111112

