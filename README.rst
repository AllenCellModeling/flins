============
flins
============


.. .. image:: https://img.shields.io/pypi/v/flins.svg
..         :target: https://pypi.python.org/pypi/flins

.. image:: https://github.com/AllenCellModeling/flins/workflows/Documentation/badge.svg
        :target: https://AllenCellModeling.github.io/flins
        :alt: Documentation Status

.. image:: https://codecov.io/gh/AllenCellModeling/flins/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/AllenCellModeling/flins
        :alt: Code Coverage Status

.. image:: https://github.com/AllenCellModeling/flins/workflows/Build%20Master/badge.svg
        :target: https://github.com/AllenCellModeling/flins/actions
        :alt: Continuous Integration Status

Overview
--------

``flins`` is a spatial, agent-based simulation of the movement of actin, cross-linkers, and adhesions in a pre-myofibril. It explores the conditions necessary to generate the emergent organization we see in differentiating and developing muscle cells. 

As a spatially explicit simulation ``flins`` recreates the movement and force responses of its proteins by treating them as springs, subject to deformation and able to generate force with variable rest lengths. Connectivity and binding withing this system of proteins is controlled by binding sites distributed along the proteins that link springs together. The kinetics of these binding sites are dependent upon the forces their parent proteins are subjected to. This produces a network of springs that transmits and generates forces with connectivity that changes depending on stochastic kinetics and the current deformation within the system. 
  
Minimal example::

    import flins
    world = flins.construct.create_test_world(
        radius=1,     # how wide the world is
        span=10000,   # how long
        n_actin=5,    # actins per slice
        n_actinin=20, # crosslinkers per slice
        n_motors=10,  # motors per slice
        )    
    world.step()

* Documentation: https://AllenCellModeling.github.io/flins
