import ifg
import re

# for line in open('SMILES.txt', 'r'):
#     lineInfo = re.compile(r'\S+').findall(line)
#     smiles = lineInfo[0]
#     print("EVALUATING ", lineInfo[1])
#     functionalGroupData = ifg.ifg(smiles)
#     print(functionalGroupData)
functionalGroupData = ifg.ifg("[N-]=[N+]=NC1=NN=NN1N=C1CCC(CC1)=NN1N=NN=C1N=[N+]=[N-]")
print(functionalGroupData)
print("[N-]=[N+]=NC1=NN=NN1N=C1CCC(CC1)=NN1N=NN=C1N=[N+]=[N-]")
