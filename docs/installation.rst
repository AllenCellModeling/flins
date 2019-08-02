.. highlight:: shell

============
Installation
============


Updateable install
------------------

To install ``stress_fiber`` in an updatable and editable form, run these commands in your terminal:

.. code-block:: console

    $ git clone git://github.com/AllenCellModeling/stress_fiber
    $ cd stress_fiber
    $ pip install -e .[all]

This is the preferred method to install ``stress_fiber`` in its current alpha form as it will allow updates via:

.. code-block:: console

    $ git pull

Direct pip install
------------------

It is also possible to install ``stress_fiber`` directly from the `Github repo`_ via:

.. code-block:: console

    $ pip install git+https://github.com/AllenCellModeling/stress_fiber.git

.. _Github repo: https://github.com/AllenCellModeling/stress_fiber
.. _tarball: https://github.com/AllenCellModeling/stress_fiber/tarball/master

Visualization support
---------------------

Display of SVG visualizations within notebooks is supported by Ipython's core display functions, but to render these to PNGs and then to animations we rely on `librsvg`_ and `ffmpeg`_. Installation with Anaconda is accomplished via:

.. code-block:: console

    $ conda install -c conda-forge librsvg ffmpeg

.. _librsvg: https://en.wikipedia.org/wiki/Librsvg
.. _ffmpeg: https://ffmpeg.org
