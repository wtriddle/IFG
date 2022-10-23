.. _collection-script-ref:

SMILES-Based Functional Group Data Collection Script
====================================================

A script called ``main.py`` (or ``main.ipynb``) generates functional group data from a set of target SMILES code structures
and exports the data to an output excel file. The ``main.py`` is a concise single runtime script, while the ``main.ipynb`` is a 
verbose cell-by-cell execution environment designed to help individuals better understand python code and the steps 
involved in the ``main.py`` script. 

.. note:: 
    If you would like to add more functional groups to the identification process,
    please add them to the `ifg/chem/data/functional_group_smiles_codes.csv` file in the repository in the format shown in the file.

Set Up The Script
-----------------

First, follow the :ref:`installation-ref` steps to set up the scripting environment.

SMILES codes
++++++++++++

To choose a valid set of SMILES-defined molecular structures to draw functional group data from, create a SMILES csv file in `/scripts/smiles`
that reflects the model of the other SMILES csv files already present in the folder. You can target this new file by adjusting the 
``STRUCTURES_PATH`` variable to target your SMILES codes, located in the ``main.py`` (or ``main.ipynb``) script. Please see the
:ref:`supported smiles codes <smiles-code-support-ref>` to see what types of SMILES codes can be used.

.. note:: 
    If the csv file you wish to use is in a different format, either 
    convert the csv file format to fit the format shown in the `/scripts/smiles`, or feel free to adjust 
    the ``main.py`` python script itself to load data from an entirley different format into the script.
    For example, data loading from excel sheets, a text file, with other variables, 
    and so on.

Excel Sheet
+++++++++++

The excel sheet that is created as an output from the ``main.py`` (or ``main.ipynb``) scripts can be configured by name or file destination in either 
script via the ``MAIN_OUTPUT_PATH`` variable (shown in both scripts).

Run The Script
--------------

To run the ``main.py`` script, run the following `CLI <https://en.wikipedia.org/wiki/Command-line_interface>`_ command in your `/ifg/scripts` directory:

.. code-block:: 
   :caption: Script Run Command

    (active-ifg-poetry-environment) C:\IFG\ifg\scripts> poetry run python main.py

The script targets your desired csv set of SMILES codes and processes each one through the functional group identification process, then produces an
excel sheet of the collected data afterwards. 

To run the ``main.ipynb``, open a jupyter notebooks-like environment such as jupyter-labs (which is available with this repository through the poetry environment)
and execute each cell. Note that the set-up is the same for this script and will help show users how the python script works step by step.

View The Output
---------------

The output of the ``main.py`` (or ``main.ipynb``) script is an excel file at the configured output ``MAIN_OUTPUT_PATH`` file path with two excel sheets called ``all_data`` and ``exact_data``.
Each sheet lists the functional groups and rings present per structure according to the `all` and `exact` functional group type definitions for molecules.
Read about the two types in the `Script Outputs`_ section below.


.. _fg-definitions-ref:

Script Outputs
--------------

The script outputs two data sheets in its output excel file according to the two types of functional group definitions in this project: ``all`` and ``exact``

``All Definition``
    Includes overlapped functional groups in the overall number of times a particular functional group appears in a structure. 
    A common example is including a ketone count while it overlaps with an ester.


``Exact Definition``
    Excludes overlapped functional groups in the overall number of times a particular functional group appears in a structure. 
    A common example is excluding a ketone count while it overlaps with an ester. 



