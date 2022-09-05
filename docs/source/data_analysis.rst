.. _analysis-script-ref:

Analyze Functional Groups For Statistical Trends Script
=======================================================

The ``analysis.py`` script generates an excel data sheet with two or more statistical matricies 
indicating the relationship between functional groups and the optical molecular bandgap (also called the bandgap). 
View the `Script Outputs`_ section below to see what statistics are computed by this script.

.. note:: 
    This script uses the ``exact`` definition of functional groups during statistical analysis.
    See the :ref:`fg-definitions-ref` of ``main.py`` for detailed definition information.


Set Up The Script
-----------------

Follow the steps in :ref:`script-setup-ref` to set up the scripting environment

Configure The Script
--------------------

To adjust the relevant file paths for this script, view the :ref:`configurations-ref` doc page.
To create a new relational matrix, add a new set of functional groups to the ``SETS`` data structure in the ``analysis.py``. 
View the relational matrix section in the `Script Outputs`_ section below for details about the sets.


Run The Script
--------------

To run the script, ensure that all python packages have been installed, and that the bandgap data is loaded properly into the script.
Then, run the following command in a command line interface

.. code-block:: 
    :caption: CLI Script Run Command

    C:\IFG> poetry run python analysis.py

View the Outputs
----------------

Once the script has completed running, you can view the output excel file based on the target excel output path that
you configured before running the script, ``ANALYSIS_OUTPUT_PATH``


Script Outputs
--------------

There are two types of matricies that this script outputs: relational matricies and an instance matrix:

`Instance Matrix`
    A matrix that counts the number of structures which exhibit a particular functional group
    and exhibit a bandgap in a particular range.

`Relational Matrix`
    A matrix that counts the number of structures which exhibit a particular set of functional groups
    and exhibit a bandgap in a particular range.

.. note:: 
    A relational matrix is generated for each functional group set listed in the ``SETS`` data structure by creating 
    its `powerset <https://en.wikipedia.org/wiki/Power_set>`_, then for each set in the powerset, querying for all structures 
    which explicitly have those functional groups and explicitly do not have the remaining complemented structures as well. 
    That is, for the set of all structures which contain the functional groups listed in the ``SETS`` entry,
    a functional groups categorical column counts all structures which have those functional groups
    present `and` which do `not` have any of the other remaining functional groups from the original set listed in the ``SETS`` entry.

The total number of structures, bandgap mean and bandgap standard deviations per column are computed and shown in both types of matricies.
Additionally, each relational matrix includes a stacked percentage chart for the distribution of functional group cateogires in each bandgap range.