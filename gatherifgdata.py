from ifg import ifg
from dataManipulators import dataManipulators

class gatherData(dataManipulators):

	def __init__(self, smilescode):
		self.smiles = smilescode
		self.functionalGroupData = self.functionalGroups(smilescode)

	def functionalGroups(self, SMILES):
	
		# Individual SMILEScode portions are split by periods
		# Overaching functional groups, no including all smaller component groups
		containmentFGDict = {}
		allFGDict = {}  # Every functional group that appears, regardless of size or containment in other larger groups
		structures = SMILES.split('.')
		for componentSmiles in structures:
			ifgObject = ifg(componentSmiles)
			data = ifgObject.FGDATAFINAL # Retrieve functional groups
			# Add slice of total structure to total functional groups present
			self.addDictData(data[1], containmentFGDict)
			self.addDictData(data[3], allFGDict)
		return (containmentFGDict, allFGDict)