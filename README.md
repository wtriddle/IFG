# IFG (Identify Functional Groups)

Python algorithms which extract functional group data from hydrogen-suppressed 
SMILES (Simplified molecular-input line-entry system) codes <br>
Project and Research Overview: https://youtu.be/yOdvyQ0seAc <br>
<p align="center">
  <!-- <img src="https://github.com/wtriddle/IFG/blob/master/MoleculeGifSmall.gif" />
  <img src="https://github.com/wtriddle/IFG/blob/master/MoleculeGifSmall.gif" /> -->
</p>

# Installation

To set up IFG for functional group data collection and analysis from SMILES codes, follow the steps below:

1. Install [git](<https://git-scm.com/downloads>)
2. Install [poetry](<https://python-poetry.org/docs/>)
3. Clone this github repository into a desired folder with: 

    ```cmd
      C:\> git clone https://github.com/wtriddle/IFG.git
    ```

4. Install the python packages required for IFG with poetry:

    ```cmd
      C:\IFG> poetry install
    ```
    
5. Add IFG/src to your system PYTHONPATH environemnt variable. See this [tutorial](https://www.simplilearn.com/tutorials/python-tutorial/python-path)

You are now ready to use IFG! Read usage to find out how to use the program.

# Usage

There are two scripts that are available in the root directory of IFG to retireve and analyze functional group data.

## Main.py

The main script creates an excel data sheet with functional groups data by using the algorithm source code.
Three file path variables in ``config.py`` are relevant for the script to run: ``STRUCTURES_PATH``, ``FUNCTIONAL_GROUPS_PATH``, and ``MAIN_OUTPUT_PATH``. Please view the ``config.py`` file for descriptions of each variable. 

To execute the ``main.py`` file with the desires structure list, functional groups list, and excel file output path, 
run the following command after installation:

```cmd
  poetry run python main.py
```

Allow the script to run and view the output in the desired excel sheet location. 
For more information about the output, please read the paper.

## Analysis.py

The analysis script creates an excel data sheet with functional groups property relationship data by using data from ``main.py`` 
in a staistical manner. This is a detailed procedure which has extensive discussion in the paper, so please read the paper for this.
Three file path variables in ``config.py`` are relevant for the script to run: ``BANDGAPS_PATH``, ``ANALYSIS_PATH``, and ``MAIN_OUTPUT_PATH``. Please view the ``config.py`` file for descriptions of each variable. 

To execute the ``analysis.py`` file with the desires structure list, functional groups list, and excel file output path, 
run the following command after installation:

```cmd
  poetry run python analysis.py
```

Allow the script to run and view the output in the desired excel sheet location. 
For more information about the output, please read the paper.

## Software Interface

If the library is to be ported into other code, then utilze the following example to do so:

Molecule Class:
```python
    >>> from Molecule import Molecule
    >>> mol = Molecule('O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O', 'ABEGOH')
    >>> print(mol)
    ABEGOH : O=C1NC2C(N(CN2N(=O)=O)N(=O)=O)N1N(=O)=O
```

# Contributing

When contributing, open an issue for suggestions and follow the commenting style observed across the repository <br>
New formats of SMILES codes can be supported by updating and contributing to the repository <br>

# License
[MIT](https://choosealicense.com/licenses/mit/)
