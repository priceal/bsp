

# analyze a dataframe on 'RE' for putative/ORF content

df = siteDf

putative=df[ df['RE'].str.endswith('p') ]
orfs=df[ df['RE'].str.count('orf')>0 ]
noOrf=df[ df['RE'].str.count('orf')==0 ]

putNoOrf = putative[ putative['RE'].str.count('orf')==0 ]
putOrf = putative[ putative['RE'].str.count('orf')>0 ]

noPutOrf = orfs[ ~orfs['RE'].str.endswith('p') ] 
noPutNoOrf = noOrf[ ~noOrf['RE'].str.endswith('p')]
