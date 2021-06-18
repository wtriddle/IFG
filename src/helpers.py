""" Module where helping functiona, used across the program, can be defined """

from collections import defaultdict


def createFgDataDict(functionalGroups):
    """ Return a dictionary count of every functional group in a valid functional groups list

        functionalGroups (list) : List of Molecule obejects that represent the functioanl groups in a smiles code

        Notes:
            defaultdict allows for an unrecognized key to be initalized as an integer, then instantly incremented.
    """
    fgDataDict = defaultdict(int)
    for group in functionalGroups:
        fgDataDict[group.NAME] += 1
    return fgDataDict
