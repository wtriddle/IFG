def createFgDataDict(functionalGroups):
    fgDataDict = {}
    for group in functionalGroups:
        if group.NAME in fgDataDict.keys():
            fgDataDict[group.NAME] += 1
        else:
            fgDataDict.update({group.NAME: 1})
    return fgDataDict
