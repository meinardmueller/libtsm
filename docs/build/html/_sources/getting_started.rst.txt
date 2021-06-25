Getting Started
===============

You can install libtsm using Anaconda:

.. code-block:: bash

    git clone https://github.com/meinardmueller/libtsm.git
    conda create -n libtsm python=3.8
    conda activate libtsm
    cd ./libtsm
    pip install .

For development purposes, install:

.. code-block:: bash

    pip install .[dev]
    pip install .[tests]
    pip install .[docs]
