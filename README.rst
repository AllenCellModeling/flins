============
stress_fiber
============


.. .. image:: https://img.shields.io/pypi/v/stress_fiber.svg
..         :target: https://pypi.python.org/pypi/stress_fiber

.. image:: https://readthedocs.org/projects/stress-fiber/badge/?version=latest
        :target: https://stress-fiber.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://codecov.io/gh/AllenCellModeling/stress_fiber/branch/master/graph/badge.svg
        :target: https://codecov.io/gh/AllenCellModeling/stress_fiber
        :alt: Code Coverage Status

.. image:: https://travis-ci.com/AllenCellModeling/stress_fiber.svg?branch=master
        :target: https://travis-ci.com/AllenCellModeling/stress_fiber
        :alt: Continuous Integration Status

Overview
--------

``stress_fiber`` is a spatial, agent-based simulation of the movement of actin, cross-linkers, and adhesions in a pre-myofibril. It explores the conditions necessary to generate the emergent organization we see in differentiating and developing muscle cells. 

As a spatially explicit simulation ``stress_fiber`` recreates the movement and force responses of its proteins by treating them as springs, subject to deformation and able to generate force with variable rest lengths. Connectivity and binding withing this system of proteins is controlled by binding sites distributed along the proteins that link springs together. The kinetics of these binding sites are dependent upon the forces their parent proteins are subjected to. This produces a network of springs that transmits and generates forces with connectivity that changes depending on stochastic kinetics and the current deformation within the system. 
  
* Documentation: https://stress-fiber.readthedocs.io.
