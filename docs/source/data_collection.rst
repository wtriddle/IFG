.. _collection-script-ref:

Collect Functional Groups Data From SMILES Structures Script
============================================================

A script called ``main.py`` generates functional group data from a set of target SMILES code structures
and exports the data to an output excel file.

Set Up The Script
-----------------

Follow the steps in :ref:`script-setup-ref` to set up the scripting environment


Configure The Script
--------------------

The ``main.py`` script is configurable `via` its input structures list and its output excel file location. 
View the :ref:`configurations-ref` doc page for the ``main.py`` file configuration information.

Run The Script
--------------

To run the script, ensure the scripting environment has been set up,
then run the following command line interface command in your IFG root directory:

.. code-block:: 
   :caption: Script Run Command

    C:\IFG> poetry run python main.py

This will execute the data collection process for all target SMILES structures and generate an excel sheet at the configured output file path with the identified data.

.. warning::
    Ensure that all SMILES codes are hydrogen-suppressed.
    Otherwise, the script will fail.

View The Output
---------------

The output of the ``main.py`` script is an excel file at the configured output path with two sheets: ``all_data`` and ``exact_data``.
Each sheet lists the functional groups and rings present per structure according to the definitions of `all` and `exact`.
View the `all` and `exact` definitions in the `Script Outputs`_ section below.


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



