Box Least Squares
=================

These are Python bindings for the original Fortran implementation of Box Least
Squares (BLS) algorithm from `Kovács et al. (2002)
<http://arxiv.org/abs/astro-ph/0206099>`_.

Installation
------------

Clone the source code from `the GitHub repository
<https://github.com/dfm/python-bls>`_ and install using the standard python
tools:

.. code-block:: bash

    git clone https://github.com/dfm/python-bls.git
    cd python-bls
    python setup.py install

For testing in development, you can use

.. code-block:: bash

    python setup.py build_ext --inplace

to build the bindings directly in the ``bls`` directory.

Authors
-------

These Python bindings were developed—building directly on the code released by
`Kovács <http://www.konkoly.hu/staff/kovacs/eebls.f>`_—at the `SAMSI
<http://samsi.info>`_ workshop `Modern Statistical and Computational Methods
for Analysis of Kepler
<http://www.samsi.info/working-groups/kepler-working-group>`_ by

* `Ruth Angus (Oxford) <https://github.com/RuthAngus>`_
* `Dan Foreman-Mackey (NYU) <https://github.com/dfm>`_

License
-------

The Python bindings are licensed under the `MIT License
<https://github.com/dfm/python-bls/blob/master/LICENSE>`_.

Basic Usage
-----------

TODO: describe Pythonic binding.

Advanced Usage
--------------

You can get direct access to `the Fortran subroutine provided by Kovács
<http://www.konkoly.hu/staff/kovacs/eebls.f>`_ through the ``eebls()``
function:

.. code-block:: python

    import bls
    results = bls.eebls(time, flux, u, v, nf, fmin, df, nb, qmi, qma)

where

* ``time`` is an ``N``-dimensional array of timestamps for the light curve,
* ``flux`` is the ``N``-dimensional light curve array,
* ``u`` and ``v`` are ``N``-dimensional empty work arrays,
* ``nf`` is the number of frequency bins to test,
* ``fmin`` is the minimum frequency to test,
* ``df`` is the frequency grid spacing,
* ``nb`` is the number of bins to use in the folded light curve,
* ``qmi`` is the minimum transit duration to test, and
* ``qma`` is the maximum transit duration to test.

The returned values are

.. code-block:: python

    power, best_period, best_power, depth, q, in1, in2 = results

where

* ``power`` is the ``nf``-dimensional power spectrum array at frequencies ``f
  = fmin + arange(nf) * df``,
* ``best_period`` is the best-fit period in the same units as ``time``,
* ``best_power`` is the power at ``best_period``,
* ``depth`` is the depth of the transit at ``best_period``,
* ``q`` is the fractional transit duration,
* ``in1`` is the bin index at the start of transit, and
* ``in2`` is the bin index at the end of transit.
