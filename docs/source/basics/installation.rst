.. _installation:

============
Installation
============

Download
--------

The hydrangea library code is hosted on github. To get a copy, clone the repository with the command::
	
	$ git clone https://github.com/ymbahe/hydrangea.git

This will download all the library files to a subdirectory named ``hydrangea`` inside the current working directory.
If you already have an (older) version of the library, you can update it with the command::

    $ git pull

from *within* the ``hydrangea`` sub-directory to update the files to the latest version.

.. _install:

Installation
------------

Once you have downloaded (or updated) the repository, it can be installed with the following steps

#. Go to the base directory of the library package::

   $ cd hydrangea

#. [Optional] Edit file ``./hydrangea/local.py`` to specify the base directory of the Hydrangea/C-EAGLE simulation
   data on your system
   (i.e. change the default ``BASE_DIR = None`` to ``BASE_DIR = '/path/to/base/dir'``).
   This makes it slightly easier to specify simulation directories when using the library, and can be omitted when
   updating the library.

#. [Optional] Compile the external C libraries (only required for speed-up in a few specialised functions)::

    $ ./c_install.sh

#. Install the hydrangea package in your python installation::
   
   $ pip install --user .

You can now import the library in any python script by adding the line::

  import hydrangea [as hy]  # Part in square brackets is optional


