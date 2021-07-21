""" Module where helping functions are defined"""

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


def formatSmiles(smiles):
    """ Remove [H+] symbols entirley from a smiles code and returns a DLA-SAR converted smiles code 

        Return a double letterd atom (DLA) tranformed SMILES code using single atom representations (SAR)
        Required to perform sybmol by symbol analysis on the SMILES
    """

    DLA_TO_SAR = {                                  
        "Br": "X",                                       
        "Cl": "Z"
    }

    while 'H' in smiles:                               # Filter H if present, otherwise skip
        for pos, symbol in enumerate(smiles):

            if symbol == '[':
                startBracketPos = pos

            if symbol == 'H':
                cutPos = pos
                while smiles[cutPos] != ']':
                    cutPos += 1
                smiles = smiles[0:startBracketPos] + \
                    smiles[startBracketPos+1] + smiles[cutPos+1:len(smiles)]
                break
    
    if "Br" in smiles or "Cl" in smiles:
        pos = -1                               
        newSmiles = ""
        while pos != len(smiles) - 1:       

            pos += 1
            DLA = smiles[pos:pos+2]         # 2 char string is a potential DLA

            if DLA in DLA_TO_SAR:                # Convert DLA'S via the legend if matched
                newSmiles += DLA_TO_SAR[DLA]     # Add SAR in place of DLA
                pos += 1                     

            else:                                # No DLA match keeps smiles the same
                newSmiles += smiles[pos]
        return newSmiles

    return smiles
