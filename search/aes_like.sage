from sage.crypto.boolean_function import BooleanFunction
import sys
load("aux_function.sage")
load("MILP_function.sage")

#Modelization of an AES-like block cipher
#The code is a bit messy with duplicated code but it was easier than setting control ifs everywhere, and it works
#It is probably better/easier to use the functions provided in the c++ code

def computeInequalitiesMC(M):
	#Precompute inequalities for MC
	M_nrows = M.nrows()
	M_ncols = M.ncols()

	V = vectorF2n(M_nrows)
	L = [M*vector(x) for x in V]
	MCbox = [None for _ in range(1 << M_nrows)]
	for i in range(1 << M_nrows):
		x = 0
		for j in range(M_nrows):
			if L[i][j] == 1:
				x += 2^j
		MCbox[i] = x

	(BPR,anf) = anfSbox(MCbox)
	TMC = divPropANFBinTable(anf)
	return sboxReducedInequalities(TMC)

def applySbox(m, xr, yr, ineq, sboxSize, nbSbox):
	#Add the constraints for the Sbox layer from xr to yr, using the inequalities in ineq
	#Constraints for Sbox
	for j in range(nbSbox):
		inputvar = [xr[sboxSize*j+i] for i in range(sboxSize)]
		outputvar = [yr[sboxSize*j+i] for i in range(sboxSize)]
		addSboxConstr(m, ineq, inputvar, outputvar)

def applyMCAsLbox(m, wr, xrp1, ineqMC, M_ncols, nbCol, sboxSize):
	colSize = M_ncols*sboxSize
	for col in range(nbCol):
		for offset in range(sboxSize):
			inputvar = [wr[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
			outputvar = [xrp1[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
			addSboxConstr(m,ineqMC,inputvar,outputvar)

def applyMCAsCopyXor(m, wr, xrp1, M, nbCol, r):
	M_nrows = M.nrows()
	M_ncols = M.ncols()
	for col in range(nbCol):
		#First, create the temporary variables for the Copy + XOR modeling
		t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
		for i in range(M_nrows):
			for j in range(M_ncols):
				if M[i][j] == 1:
					t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r)+"_"+str(col)+"_"+str(i)+"_"+str(j))

		#Copy constraints
		for j in range(M_ncols):
			m.addGenConstrOr(wr[col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

		#XOR constraints
		for i in range(M_nrows):
			m.addConstr(xrp1[col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))


def genBaseModelSPARKMC(rMax, S, P, M, name, linAsSbox,noLastMC):

	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)

	M_nrows = M.nrows()
	M_ncols = M.ncols()
	if linAsSbox:
		#Precompute inequalities for MC
		ineqMC = computeInequalitiesMC(M)

	#Some useful constants
	sboxSize = anf[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	nbCol = blockSize//M_nrows #Number of column in the state if linAsSbox = false
	if linAsSbox:
		nbCol = nbSbox//M_nrows
	colSize = M_ncols*sboxSize

	#Compute the permutation extended at the bit level
	expandP = [0 for i in range(blockSize)]
	for isbox in range(nbSbox):
		for j in range(sboxSize):
			expandP[isbox*sboxSize+j] = P[isbox]*sboxSize+j

	#Create the models and the variables
	m = Model(name)

	if noLastMC:
		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
		#No MC on last round, stops at y[rMax-1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		
	else:
		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]

	m.update()
	#Round Function constraints (except last round)
	for r in range(rMax-1):
		#Constraints for Sbox
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(z[r][i] + k[r][i] == w[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, w[r], x[r+1], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, w[r], x[r+1], M, nbCol, r)

	#last round
	r = rMax-1
	#Constraints for Sbox
	applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

	if not noLastMC: #we need the last linear layer
		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(z[r][i] + k[r][i] == w[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, w[r], x[r+1], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, w[r], x[r+1], M, nbCol, r)


	m.update()

	return (m,x,y,z,w,k)

def genBaseModelSPARKMC_startAfterSSB(rMax, S, P, M, name, linAsSbox,noLastMC):

	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)

	M_nrows = M.nrows()
	M_ncols = M.ncols()
	if linAsSbox:
		#Precompute inequalities for MC
		ineqMC = computeInequalitiesMC(M)

	#Some useful constants
	sboxSize = anf[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	nbCol = blockSize//M_nrows #Number of column in the state if linAsSbox = false
	if linAsSbox:
		nbCol = nbSbox//M_nrows
	colSize = M_ncols*sboxSize

	#Compute the permutation extended at the bit level
	expandP = [0 for i in range(blockSize)]
	for isbox in range(nbSbox):
		for j in range(sboxSize):
			expandP[isbox*sboxSize+j] = P[isbox]*sboxSize+j

	#Create the models and the variables
	m = Model(name)

	
	#first variables for the first round
	y0 = [m.addVar(vtype=GRB.BINARY, name="y_"+str(0)+"_"+str(j)) for j in range(blockSize)]
	z0 = [m.addVar(vtype=GRB.BINARY, name="z_"+str(0)+"_"+str(j)) for j in range(blockSize)]
	w0 = [m.addVar(vtype=GRB.BINARY, name="w_"+str(0)+"_"+str(j)) for j in range(blockSize)]
	k0 = [m.addVar(vtype=GRB.BINARY, name="k_"+str(0)+"_"+str(j)) for j in range(blockSize)]

	if noLastMC:
		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
		#No MC on last round, stops at y[rMax-1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		
	else:
		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]

	m.update()

	#Partial first round
	#Constraints for Permutation
	for i in range(blockSize):
		m.addConstr(z0[expandP[i]] == y0[i])

	#Constraints for ARK
	for i in range(blockSize):
		m.addConstr(z0[i] + k0[i] == w0[i])

	#Constraints for MixColumn
	if linAsSbox:
		applyMCAsLbox(m, w0, x[0], ineqMC, M_ncols, nbCol, sboxSize)

	else:
		applyMCAsCopyXor(m, w0, x[0], M, nbCol, r)

	#Round Function constraints (except last round)
	for r in range(rMax-1):
		#Constraints for Sbox
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(z[r][i] + k[r][i] == w[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, w[r], x[r+1], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, w[r], x[r+1], M, nbCol, r)

	#last round
	r = rMax-1
	#Constraints for Sbox
	applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

	if not noLastMC: #we need the last linear layer
		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(z[r][i] + k[r][i] == w[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, w[r], x[r+1], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, w[r], x[r+1], M, nbCol, r)

	

	m.update()

	return (m,x,y,z,w,k)

def genBaseModelSPMCARK(rMax, S, P, M, name, linAsSbox,noLastMC,dedicatedPRESENTlastLayer):

	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)

	M_nrows = M.nrows()
	M_ncols = M.ncols()
	if linAsSbox:
		#Precompute inequalities for MC
		ineqMC = computeInequalitiesMC(M)

	#Some useful constants
	sboxSize = anf[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	nbCol = blockSize//M_nrows #Number of column in the state if linAsSbox = false
	if linAsSbox:
		nbCol = nbSbox//M_nrows
	colSize = M_ncols*sboxSize

	#Compute the permutation extended at the bit level
	expandP = [0 for i in range(blockSize)]
	for isbox in range(nbSbox):
		for j in range(sboxSize):
			expandP[isbox*sboxSize+j] = P[isbox]*sboxSize+j

	#Create the models and the variables
	m = Model(name)

	if noLastMC:
		#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
		#No MC on last round, stops at y[rMax-1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
	else:
		#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]

	m.update()

	#Round Function constraints (except last round)
	for r in range(rMax-1):
		#Constraints for Sbox
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, z[r], w[r], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, z[r], w[r], M, nbCol, r)

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(w[r][i] + k[r][i] == x[r+1][i])

	#last round
	r = rMax-1
	#Constraints for Sbox
	if not dedicatedPRESENTlastLayer:
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)
	else:
		anf[1] = anf[0] + anf[1] + anf[3]
		anf[3] = anf[0] + anf[2] + anf[3]
		T = divPropANFBinTable(anf)
		ineqSboxDecicated = sboxReducedInequalities(T)
		applySbox(m, x[r], y[r], ineqSboxDecicated, sboxSize, nbSbox)

	if not noLastMC: #we need the last linear layer
		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, z[r], w[r], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, z[r], w[r], M, nbCol, r)

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(w[r][i] + k[r][i] == x[r+1][i])

	m.update()

	return (m,x,y,z,w,k)


def genBaseModelSPMCARK_startAfterSSB(rMax, S, P, M, name, linAsSbox,noLastMC,dedicatedPRESENTlastLayer):

	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)

	M_nrows = M.nrows()
	M_ncols = M.ncols()
	if linAsSbox:
		#Precompute inequalities for MC
		ineqMC = computeInequalitiesMC(M)

	#Some useful constants
	sboxSize = anf[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	nbCol = blockSize//M_nrows #Number of column in the state if linAsSbox = false
	if linAsSbox:
		nbCol = nbSbox//M_nrows
	colSize = M_ncols*sboxSize

	#Compute the permutation extended at the bit level
	expandP = [0 for i in range(blockSize)]
	for isbox in range(nbSbox):
		for j in range(sboxSize):
			expandP[isbox*sboxSize+j] = P[isbox]*sboxSize+j

	#Create the models and the variables
	m = Model(name)

	
	#No MC on last round, stops at y[rMax-1]
	#first variables for the first (partial) round
	y0 = [m.addVar(vtype=GRB.BINARY, name="y_"+str(0)+"_"+str(j)) for j in range(blockSize)]
	z0 = [m.addVar(vtype=GRB.BINARY, name="z_"+str(0)+"_"+str(j)) for j in range(blockSize)]
	w0 = [m.addVar(vtype=GRB.BINARY, name="w_"+str(0)+"_"+str(j)) for j in range(blockSize)]
	k0 = [m.addVar(vtype=GRB.BINARY, name="k_"+str(0)+"_"+str(j)) for j in range(blockSize)]

	if noLastMC:
		#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]

	else:
		#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]

	m.update()

	#Partial round
	#Constraints for Permutation
	for i in range(blockSize):
		m.addConstr(z0[expandP[i]] == y0[i])

	#Constraints for MixColumn
	if linAsSbox:
		applyMCAsLbox(m, z0, w0, ineqMC, M_ncols, nbCol, sboxSize)

	else:
		applyMCAsCopyXor(m, z0, w0, M, nbCol, r)

	#Constraints for ARK
	for i in range(blockSize):
		m.addConstr(w0[i] + k0[i] == x[0][i])


	#Round Function constraints (except last round)
	for r in range(rMax-1):
		#Constraints for Sbox
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)

		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, z[r], w[r], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, z[r], w[r], M, nbCol, r)

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(w[r][i] + k[r][i] == x[r+1][i])

	#last round
	r = rMax-1
	#Constraints for Sbox
	if not dedicatedPRESENTlastLayer:
		applySbox(m, x[r], y[r], ineqSbox, sboxSize, nbSbox)
	else:
		anf[1] = anf[0] + anf[1] + anf[3]
		anf[3] = anf[0] + anf[2] + anf[3]
		T = divPropANFBinTable(anf)
		ineqSboxDecicated = sboxReducedInequalities(T)
		applySbox(m, x[r], y[r], ineqSboxDecicated, sboxSize, nbSbox)

	if not noLastMC: #we need the last linear layer
		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z[r][expandP[i]] == y[r][i])

		#Constraints for MixColumn
		if linAsSbox:
			applyMCAsLbox(m, z[r], w[r], ineqMC, M_ncols, nbCol, sboxSize)

		else:
			applyMCAsCopyXor(m, z[r], w[r], M, nbCol, r)

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(w[r][i] + k[r][i] == x[r+1][i])
	

	m.update()

	return (m,x,y,z,w,k)



load("aes_like_parameters.sage");
M = matrix(GF(2), M)
if not keyAfterMC:
		if startAfterSSB:
			(m,x,y,z,w,k) = genBaseModelSPARKMC_startAfterSSB(rMax, S, P, M, name, linAsSbox,noLastMC)
		else:
			(m,x,y,z,w,k) = genBaseModelSPARKMC(rMax, S, P, M, name, linAsSbox,noLastMC)
else:
	if startAfterSSB:
		(m,x,y,z,w,k) = genBaseModelSPMCARK_startAfterSSB(rMax, S, P, M, name, linAsSbox,noLastMC,dedicatedPRESENTlastLayer)
	else:
		(m,x,y,z,w,k) = genBaseModelSPMCARK(rMax, S, P, M, name, linAsSbox,noLastMC,dedicatedPRESENTlastLayer)
m.write(name+".mps")