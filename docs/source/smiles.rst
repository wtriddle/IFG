.. _smiles-code-support-ref:

SMILES Code Support
===================

This software has been designed for canonicalized non-isomeric non-isotopic hydrogen-suppressed SMILES codes to execute its internal
processing for functional groups collection. To see examples of these types of SMILES codes, please view the smiles.csv file in the repository
under `ifg/scripts/smiles/smiles.csv`.

See examples of unsupported SMILES codes with type and property references below:

.. list-table:: Unsupported SMILES type references 
   :widths: 25 25 50
   :header-rows: 1

   * - SMILES Type
     - Description
     - Examples
   * - Isomeric
     - Chirality and Stereochemitry properties
     - F/C=C/F, N[C@@](F)(C)C(=O)O
   * - Isotopic
     - Isotopic mass property
     - [14c]1ccccc1
   * - Hydrogen-inclusive
     - Explicit hydrogen property
     - C1=CN=C[NH]C(=O)1, [CH3]C(=O)[OH]
   * - Disconnected Structures
     - Multiple molecule graphs
     - [I-].[Na+].C=CCBr, C1.C1
   * - Percentage Defined Rings
     - 10 Ring exceeding structures
     - C2%13%24
  
If you would like to contribute to the project to extend the range of these supported SMILES codes, please see the :ref:`contributing-ref`.
