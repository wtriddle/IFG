import re

class dataManipulators:

	# Create dictionary which counts the number of list elements into dictionary
	def createFGDataDict(self, dataList):
		dict = {}
		for group in dataList:
			keys = []
			for key in dict.keys():
				keys.append(key)
			if group[0] not in keys:
				dict.update({group[0] : 1})
			else:
				dict[group[0]] += 1
			del(keys) # Delete the data after usage to prevent potential data leaks
		return dict

	# Add the contents of additiveDict to dataDict, either as new key:value pairs or to key:value+additiveValue
	def addDictData(self, additiveDict, dataDict):
		additiveKeys = additiveDict.keys()
		dataKeys = dataDict.keys()
		for addKey in additiveKeys:
			if addKey not in dataKeys:
				dataDict.update({addKey : additiveDict[addKey]})
				continue
			for key in dataKeys:
				if key == addKey:
					dataDict[key] += additiveDict[addKey]
