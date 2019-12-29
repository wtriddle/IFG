# IFG
(Identify Functional Groups)
A python script which extracts functional group data from SMILES (Simplified molecular-input line-entry system) codes

The data parsing functions written in IFG.py slice and group the characters in
the SMILEScode strings in various ways to obtain a string of a potential 
functional group, represented as "templates" in the file FGlist.txt. 
The template is just a string of characters, which appears in the
SMILEScode strings in different ways. The goal of this program is to
search the SMILEScode, with respect to how the ring structure and
bond types appear, for these functional group templates.


