**rfsed** is a software for receiver functions analysis and dealing with sediment effects 


License
==================

Copyright (c) 2023, Stephen Akinremi

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Installation
==================

Dependencies
-----------------
The current version works well with **Python >=3.6, <3.11**. It also require 
the following packages and their default dependencies:

* `rf <https://github.com/trichter/rf>`_
* `obspy <https://github.com/obspy/obspy/wiki>`_

Aside, these packages, other packages required for proper installation, 
testing, and running of the software are defined in the setup files. 
They are installed automatically alongside the software. 


Conda environment
------------------
A custom `conda environment <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html>`_ 
where **rfsed** can be installed is recommended:

.. code-block:: bash

    conda create -n rfsed python=3.9 -c conda-forge

Activate the newly created environment

.. code-block:: bash

    conda activate rfsed


Installing latest version from PyPi
-------------------------------------
Install the latest version from PyPi using pip:

.. code-block:: bash

    pip install rfsed


Installing development version from source
--------------------------------------------
* Clone the repository from github:

.. code-block:: bash

    git clone https://github.com/akinremisa/rfsed.git

* Change into the software directory:

.. code-block:: bash

    cd rfsed


* Install the package using pip:

.. code-block:: bash

    pip install .



Testing
-----------------
`pytest` is used for testing the software. It is installed as a dependency of the 
software.

* When `rfsed` is installed by cloning the repository from GitHub, the software can be tested thus:


In the `rfsed` directory, run the following command:

.. code-block:: bash

    cd tests

.. code-block:: bash

    pytest

This command will look for all available tests in the `tests` directory

Alternatively, you can run individual tests in the `tests` directory of the software.


* When `rfsed` is installed using `pip install rfsed`, the software can be tested thus:

The automated test files can be installed locally from the package using 
the following command in a ``python`` window:

.. code-block:: python

    from rfsed.util import install_rfsed_tests

    install_rfsed_tests(path='./rfsed_tests')

Change into the `rfsed_tests` directory and run the following command:

.. code-block:: bash

    cd rfsed_tests

.. code-block:: bash

    pytest
    

Usage
==================
Jupyter Notebooks examples are provided with the software.
These notebooks provide a step-by-step guide on how to use the software.
If **rfsed** is installed by cloning the repository from GitHub, the notebooks
can be found in the ``examples`` directory of the software. If **rfsed** is 
installed from PyPi, the notebooks can be installed locally from the package using 
the following command in a ``python`` window:

.. code-block:: python

    from rfsed.util import install_examples

    install_examples(path='./rfsed_examples')

To run the notebooks, you need to install ``jupyter``. Install from the terminal using:


.. code-block:: bash

    pip install jupyter

Then, run the following command in the terminal:    

.. code-block:: bash
    
    cd rfsed_examples

    jupyter notebook	

You can then set up your own usage of the software by following the examples in the notebooks.