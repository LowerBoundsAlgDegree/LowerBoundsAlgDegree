from sage.crypto.boolean_function import BooleanFunction
import sys
load("aux_function.sage")
load("MILP_function.sage")


#Modelization of the Super Sbox of an AES-like block cipher
#The code is a bit messy with duplicated code but it was easier than setting control ifs everywhere, and it works
#It is probably better/easier to use the functions provided in the c++ code

def genBaseModelSARKMCS(name,S,M,sboxSize,blockSize,linAsSbox):
	M_ncols = M.ncols()
	nbSbox = blockSize//sboxSize
	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)
	sboxSize = anf[0].parent().n_variables() #size of the Sbox in bits

	if linAsSbox:
		#Precompute inequalities for MC
		V = vectorF2n(4)
		L = [M*vector(x) for x in V]
		MCbox = [None for _ in range(16)]
		for i in range(16):
			x = 0
			for j in range(4):
				if L[i][j] == 1:
					x += 2^j
			MCbox[i] = x

		(BPR,anf) = anfSbox(MCbox)
		TMC = divPropANFBinTable(anf)
		ineqMC= sboxReducedInequalities(TMC)

	m = Model(name+"SSB")

	u  = [m.addVar(vtype=GRB.BINARY, name="u"+str(i)) for i in range(blockSize)]
	x  = [m.addVar(vtype=GRB.BINARY, name="x"+str(i)) for i in range(blockSize)]
	k  = [m.addVar(vtype=GRB.BINARY, name="k"+str(i)) for i in range(blockSize)]
	y  = [m.addVar(vtype=GRB.BINARY, name="y"+str(i)) for i in range(blockSize)]
	z  = [m.addVar(vtype=GRB.BINARY, name="z"+str(i)) for i in range(blockSize)]
	v  = [m.addVar(vtype=GRB.BINARY, name="v"+str(i)) for i in range(blockSize)]

	#Sbox constraints
	for j in range(nbSbox):
		inputvar  = [u[sboxSize*j+i] for i in range(sboxSize)]
		outputvar = [x[sboxSize*j+i] for i in range(sboxSize)]
		addSboxConstr(m,ineqSbox,inputvar,outputvar)

	#Constraints for ARK
	for i in range(blockSize):
		m.addConstr(x[i] + k[i] == y[i])

	#Constraints for MC
	if linAsSbox:
		for offset in range(sboxSize):
			inputvar  = [y[j*sboxSize + offset] for j in range(M_ncols)]
			outputvar = [z[j*sboxSize + offset] for j in range(M_ncols)]
			addSboxConstr(m,ineqMC,inputvar,outputvar)
	else:
		#First, create the temporary variables for the Copy + XOR modeling
		t = [[None for __ in range(M_ncols)] for _ in range(M.nrows())]
		for i in range(M.nrows()):
			for j in range(M_ncols):
				if M[i][j] == 1:
					t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(i)+"_"+str(j))

		#Copy constraints
		for j in range(M_ncols):
			m.addGenConstrOr(y[j], [t[i][j] for i in range(M.nrows()) if M[i][j] == 1])

		#XOR constraints
		for i in range(M.nrows()):
			m.addConstr(z[i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

	#Sbox constraints
	for j in range(nbSbox):
		inputvar  = [z[sboxSize*j+i] for i in range(sboxSize)]
		outputvar = [v[sboxSize*j+i] for i in range(sboxSize)]
		addSboxConstr(m,ineqSbox,inputvar,outputvar)

	return (m,u,x,k,y,z,v)

def genBaseModelSMCARKS(name,S,M,sboxSize,blockSize,linAsSbox,dedicatedPRESENTlastLayer):
	M_ncols = M.ncols()
	nbSbox = blockSize//sboxSize
	#Precompute inequalities for the sbox
	(BPR,anfS) = anfSbox(S)
	T = divPropANFBinTable(anfS)
	ineqSbox = sboxReducedInequalities(T)
	sboxSize = anfS[0].parent().n_variables() #size of the Sbox in bits

	if linAsSbox:
		#Precompute inequalities for MC
		V = vectorF2n(4)
		L = [M*vector(x) for x in V]
		MCbox = [None for _ in range(16)]
		for i in range(16):
			x = 0
			for j in range(4):
				if L[i][j] == 1:
					x += 2^j
			MCbox[i] = x

		(BPR,anf) = anfSbox(MCbox)
		TMC = divPropANFBinTable(anf)
		ineqMC= sboxReducedInequalities(TMC)


	m = Model(name+"SSB")

	u  = [m.addVar(vtype=GRB.BINARY, name="u"+str(i)) for i in range(blockSize)]
	x  = [m.addVar(vtype=GRB.BINARY, name="x"+str(i)) for i in range(blockSize)]
	k  = [m.addVar(vtype=GRB.BINARY, name="k"+str(i)) for i in range(blockSize)]
	y  = [m.addVar(vtype=GRB.BINARY, name="y"+str(i)) for i in range(blockSize)]
	z  = [m.addVar(vtype=GRB.BINARY, name="z"+str(i)) for i in range(blockSize)]
	v  = [m.addVar(vtype=GRB.BINARY, name="v"+str(i)) for i in range(blockSize)]

	#Sbox constraints
	for j in range(nbSbox):
		inputvar  = [u[sboxSize*j+i] for i in range(sboxSize)]
		outputvar = [x[sboxSize*j+i] for i in range(sboxSize)]
		addSboxConstr(m,ineqSbox,inputvar,outputvar)

	#Constraints for MC
	if linAsSbox:
		for offset in range(sboxSize):
			inputvar  = [x[j*sboxSize + offset] for j in range(M_ncols)]
			outputvar = [y[j*sboxSize + offset] for j in range(M_ncols)]
			addSboxConstr(m,ineqMC,inputvar,outputvar)
	else:
		#First, create the temporary variables for the Copy + XOR modeling
		t = [[None for __ in range(M_ncols)] for _ in range(M.nrows())]
		for i in range(M.nrows()):
			for j in range(M_ncols):
				if M[i][j] == 1:
					t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(i)+"_"+str(j))

		#Copy constraints
		for j in range(M_ncols):
			m.addGenConstrOr(x[j], [t[i][j] for i in range(M.nrows()) if M[i][j] == 1])

		#XOR constraints
		for i in range(M.nrows()):
			m.addConstr(y[i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

	#Constraints for ARK
	for i in range(blockSize):
		m.addConstr(y[i] + k[i] == z[i])

	if not dedicatedPRESENTlastLayer:
		#Sbox constraints
		for j in range(nbSbox):
			inputvar  = [z[sboxSize*j+i] for i in range(sboxSize)]
			outputvar = [v[sboxSize*j+i] for i in range(sboxSize)]
			addSboxConstr(m,ineqSbox,inputvar,outputvar)
	else:
		#dedicated last layer (y0,y0+y1+y3,y2,y0+y3+y2)
		anfS[1] = anfS[0] + anfS[1] + anfS[3]
		anfS[3] = anfS[0] + anfS[2] + anfS[3]
		T = divPropANFBinTable(anfS)
		ineqSboxDecicated = sboxReducedInequalities(T)
		for j in range(nbSbox):
			inputvar  = [z[sboxSize*j+i] for i in range(sboxSize)]
			outputvar = [v[sboxSize*j+i] for i in range(sboxSize)]
			addSboxConstr(m,ineqSboxDecicated,inputvar,outputvar)

	return (m,u,x,k,y,z,v)

load("SSBModel_parameters.sage")
if keyAfterMC:
	(m,u,x,k,y,z,v) = genBaseModelSMCARKS(name,S,M,sboxSize,blockSize,linAsSbox,dedicatedPRESENTlastLayer)
else:
	(m,u,x,k,y,z,v) = genBaseModelSARKMCS(name,S,M,sboxSize,blockSize,linAsSbox)
m.write(name+".mps")