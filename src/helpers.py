from collections import defaultdict


def createFgDataDict(functionalGroups):
    fgDataDict = defaultdict(int)
    for group in functionalGroups:
        fgDataDict[group.NAME] += 1
    return fgDataDict
