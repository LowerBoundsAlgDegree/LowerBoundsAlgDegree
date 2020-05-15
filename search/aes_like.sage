from sage.crypto.boolean_function import BooleanFunction
import sys
load("aux_function.sage")
load("MILP_function.sage")

#Modelization of an AES-like block cipher
#The code is a bit messy with duplicated code but it was easier than setting control ifs everywhere, and it works
#It is probably better/easier to use the functions provided in the c++ code

def genBaseModelSPARKMC(rMax, S, P, M, name, linAsSbox,noLastMC):

	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)

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

	#Some useful constants
	sboxSize = anf[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	M_nrows = M.nrows()
	M_ncols = M.ncols()
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
		m.update()

		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox
			for j in range(nbSbox):
				inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
				outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
				addSboxConstr(m, ineqSbox, inputvar, outputvar)


			if r < rMax-1:
				#Constraints for Permutation
				for i in range(blockSize):
					m.addConstr(z[r][expandP[i]] == y[r][i])

				#Constraints for ARK
				for i in range(blockSize):
					m.addConstr(z[r][i] + k[r][i] == w[r][i])

				#Constraints for MixColumn
				if linAsSbox:
					for col in range(nbCol):
						for offset in range(sboxSize):
							inputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							outputvar = [x[r+1][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							addSboxConstr(m,ineqMC,inputvar,outputvar)
				else:
					for col in range(nbCol):
						#First, create the temporary variables for the Copy + XOR modeling
						t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
						for i in range(M_nrows):
							for j in range(M_ncols):
								if M[i][j] == 1:
									t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r)+"_"+str(col)+"_"+str(i)+"_"+str(j))

						#Copy constraints
						for j in range(M_ncols):
							m.addGenConstrOr(w[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

						#XOR constraints
						for i in range(M_nrows):
							m.addConstr(x[r+1][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))
	else:
		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		m.update()

		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox
			for j in range(nbSbox):
				inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
				outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
				addSboxConstr(m, ineqSbox, inputvar, outputvar)

			#Constraints for Permutation
			for i in range(blockSize):
				m.addConstr(z[r][expandP[i]] == y[r][i])

			#Constraints for ARK
			for i in range(blockSize):
				m.addConstr(z[r][i] + k[r][i] == w[r][i])

			#Constraints for MixColumn
			if linAsSbox:
				for col in range(nbCol):
					for offset in range(sboxSize):
						inputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						outputvar = [x[r+1][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						addSboxConstr(m,ineqMC,inputvar,outputvar)
			else:
				for col in range(nbCol):
					#First, create the temporary variables for the Copy + XOR modeling
					t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
					for i in range(M_nrows):
						for j in range(M_ncols):
							if M[i][j] == 1:
								t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r)+"_"+str(col)+"_"+str(i)+"_"+str(j))

					#Copy constraints
					for j in range(M_ncols):
						m.addGenConstrOr(w[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

					#XOR constraints
					for i in range(M_nrows):
						m.addConstr(x[r+1][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

	m.update()

	return (m,x,y,z,w,k)

def genBaseModelSPARKMC_startAfterSSB(rMax, S, P, M, name, linAsSbox,noLastMC):

	#Precompute inequalities for the sbox
	(BPR,anf) = anfSbox(S)
	T = divPropANFBinTable(anf)
	ineqSbox = sboxReducedInequalities(T)

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

	#Some useful constants
	sboxSize = anf[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	M_nrows = M.nrows()
	M_ncols = M.ncols()
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
		#first variables for the first round
		y0 = [m.addVar(vtype=GRB.BINARY, name="y_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		z0 = [m.addVar(vtype=GRB.BINARY, name="z_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		w0 = [m.addVar(vtype=GRB.BINARY, name="w_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		k0 = [m.addVar(vtype=GRB.BINARY, name="k_"+str(0)+"_"+str(j)) for j in range(blockSize)]

		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
		#No MC on last round, stops at y[rMax-1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		m.update()

		#Partial round
		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z0[expandP[i]] == y0[i])

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(z0[i] + k0[i] == w0[i])

		#Constraints for MixColumn
		if linAsSbox:
			for col in range(nbCol):
				for offset in range(sboxSize):
					inputvar = [w0[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					outputvar = [x[0][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					addSboxConstr(m,ineqMC,inputvar,outputvar)
		else:
			for col in range(nbCol):
				#First, create the temporary variables for the Copy + XOR modeling
				t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
				for i in range(M_nrows):
					for j in range(M_ncols):
						if M[i][j] == 1:
							t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(0)+"_"+str(col)+"_"+str(i)+"_"+str(j))

				#Copy constraints
				for j in range(M_ncols):
					m.addGenConstrOr(w0[col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

				#XOR constraints
				for i in range(M_nrows):
					m.addConstr(x[0][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))



		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox
			for j in range(nbSbox):
				inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
				outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
				addSboxConstr(m, ineqSbox, inputvar, outputvar)

			if r < rMax-1:
				#Constraints for Permutation
				for i in range(blockSize):
					m.addConstr(z[r][expandP[i]] == y[r][i])

				#Constraints for ARK
				for i in range(blockSize):
					m.addConstr(z[r][i] + k[r][i] == w[r][i])

				#Constraints for MixColumn
				if linAsSbox:
					for col in range(nbCol):
						for offset in range(sboxSize):
							inputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							outputvar = [x[r+1][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							addSboxConstr(m,ineqMC,inputvar,outputvar)
				else:
					for col in range(nbCol):
						#First, create the temporary variables for the Copy + XOR modeling
						t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
						for i in range(M_nrows):
							for j in range(M_ncols):
								if M[i][j] == 1:
									t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r+1)+"_"+str(col)+"_"+str(i)+"_"+str(j))

						#Copy constraints
						for j in range(M_ncols):
							m.addGenConstrOr(w[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

						#XOR constraints
						for i in range(M_nrows):
							m.addConstr(x[r+1][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

	else:
		#first variables for the first round
		y0 = [m.addVar(vtype=GRB.BINARY, name="y_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		z0 = [m.addVar(vtype=GRB.BINARY, name="z_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		w0 = [m.addVar(vtype=GRB.BINARY, name="w_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		k0 = [m.addVar(vtype=GRB.BINARY, name="k_"+str(0)+"_"+str(j)) for j in range(blockSize)]

		#Variables are x[r] --S--> y[r] --P--> z[r] --ARK--> w[r] --MC--> x[r+1]
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

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(z0[i] + k0[i] == w0[i])

		#Constraints for MixColumn
		if linAsSbox:
			for col in range(nbCol):
				for offset in range(sboxSize):
					inputvar = [w0[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					outputvar = [x[0][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					addSboxConstr(m,ineqMC,inputvar,outputvar)
		else:
			for col in range(nbCol):
				#First, create the temporary variables for the Copy + XOR modeling
				t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
				for i in range(M_nrows):
					for j in range(M_ncols):
						if M[i][j] == 1:
							t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(0)+"_"+str(col)+"_"+str(i)+"_"+str(j))

				#Copy constraints
				for j in range(M_ncols):
					m.addGenConstrOr(w0[col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

				#XOR constraints
				for i in range(M_nrows):
					m.addConstr(x[0][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))



		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox
			for j in range(nbSbox):
				inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
				outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
				addSboxConstr(m, ineqSbox, inputvar, outputvar)

			#Constraints for Permutation
			for i in range(blockSize):
				m.addConstr(z[r][expandP[i]] == y[r][i])

			#Constraints for ARK
			for i in range(blockSize):
				m.addConstr(z[r][i] + k[r][i] == w[r][i])

			#Constraints for MixColumn
			if linAsSbox:
				for col in range(nbCol):
					for offset in range(sboxSize):
						inputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						outputvar = [x[r+1][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						addSboxConstr(m,ineqMC,inputvar,outputvar)
			else:
				for col in range(nbCol):
					#First, create the temporary variables for the Copy + XOR modeling
					t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
					for i in range(M_nrows):
						for j in range(M_ncols):
							if M[i][j] == 1:
								t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r+1)+"_"+str(col)+"_"+str(i)+"_"+str(j))

					#Copy constraints
					for j in range(M_ncols):
						m.addGenConstrOr(w[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

					#XOR constraints
					for i in range(M_nrows):
						m.addConstr(x[r+1][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

	m.update()

	return (m,x,y,z,w,k)

def genBaseModelSPMCARK(rMax, S, P, M, name, linAsSbox,noLastMC,dedicatedPRESENTlastLayer):

	#Precompute inequalities for the sbox
	(BPR,anfS) = anfSbox(S)
	T = divPropANFBinTable(anfS)
	ineqSbox = sboxReducedInequalities(T)

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

	#Some useful constants
	sboxSize = anfS[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	M_nrows = M.nrows()
	M_ncols = M.ncols()
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
		m.update()

		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox (ugly management for the dedicated present model but I need a short coding time for this)
			if r < rMax-1:
				for j in range(nbSbox):
					inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
					outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
					addSboxConstr(m, ineqSbox, inputvar, outputvar)
			elif not dedicatedPRESENTlastLayer:
				#no dedicated last layer
				for j in range(nbSbox):
					inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
					outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
					addSboxConstr(m, ineqSbox, inputvar, outputvar)
			else:
				#dedicated last layer (y0,y1+y3,y2,y3)
				anfS[1] = anfS[1] + anfS[3]
				T = divPropANFBinTable(anfS)
				ineqSboxDecicated = sboxReducedInequalities(T)
				for j in range(nbSbox):
					inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
					outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
					addSboxConstr(m, ineqSboxDecicated, inputvar, outputvar)

			if r < rMax-1:
				#Constraints for Permutation
				for i in range(blockSize):
					m.addConstr(z[r][expandP[i]] == y[r][i])

				#Constraints for MixColumn
				if linAsSbox:
					for col in range(nbCol):
						for offset in range(sboxSize):
							inputvar = [z[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							outputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							addSboxConstr(m,ineqMC,inputvar,outputvar)
				else:
					for col in range(nbCol):
						#First, create the temporary variables for the Copy + XOR modeling
						t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
						for i in range(M_nrows):
							for j in range(M_ncols):
								if M[i][j] == 1:
									t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r)+"_"+str(col)+"_"+str(i)+"_"+str(j))

						#Copy constraints
						for j in range(M_ncols):
							m.addGenConstrOr(z[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

						#XOR constraints
						for i in range(M_nrows):
							m.addConstr(w[r][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

				#Constraints for ARK
				for i in range(blockSize):
					m.addConstr(w[r][i] + k[r][i] == x[r+1][i])
	else:
		#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax+1)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		m.update()

		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox
			for j in range(nbSbox):
				inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
				outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
				addSboxConstr(m, ineqSbox, inputvar, outputvar)

			#Constraints for Permutation
			for i in range(blockSize):
				m.addConstr(z[r][expandP[i]] == y[r][i])

			#Constraints for MixColumn
			if linAsSbox:
				for col in range(nbCol):
					for offset in range(sboxSize):
						inputvar = [z[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						outputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						addSboxConstr(m,ineqMC,inputvar,outputvar)
			else:
				for col in range(nbCol):
					#First, create the temporary variables for the Copy + XOR modeling
					t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
					for i in range(M_nrows):
						for j in range(M_ncols):
							if M[i][j] == 1:
								t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r)+"_"+str(col)+"_"+str(i)+"_"+str(j))

					#Copy constraints
					for j in range(M_ncols):
						m.addGenConstrOr(z[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

					#XOR constraints
					for i in range(M_nrows):
						m.addConstr(w[r][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

			#Constraints for ARK
			for i in range(blockSize):
				m.addConstr(w[r][i] + k[r][i] == x[r+1][i])

	m.update()

	return (m,x,y,z,w,k)


def genBaseModelSPMCARK_startAfterSSB(rMax, S, P, M, name, linAsSbox,noLastMC,dedicatedPRESENTlastLayer):

	#Precompute inequalities for the sbox
	(BPR,anfS) = anfSbox(S)
	T = divPropANFBinTable(anfS)
	ineqSbox = sboxReducedInequalities(T)

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

	#Some useful constants
	sboxSize = anfS[0].parent().n_variables() #size of the Sbox
	blockSize = sboxSize*len(P) #block size
	nbSbox = len(P) #number of sboxes
	M_nrows = M.nrows()
	M_ncols = M.ncols()
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
		#No MC on last round, stops at y[rMax-1]
		#first variables for the first (partial) round
		y0 = [m.addVar(vtype=GRB.BINARY, name="y_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		z0 = [m.addVar(vtype=GRB.BINARY, name="z_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		w0 = [m.addVar(vtype=GRB.BINARY, name="w_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		k0 = [m.addVar(vtype=GRB.BINARY, name="k_"+str(0)+"_"+str(j)) for j in range(blockSize)]

		#Variables are x[r] --S--> y[r] --P--> z[r] --MC--> w[r] --ARK--> x[r+1]
		x = [[m.addVar(vtype=GRB.BINARY, name="x_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		y = [[m.addVar(vtype=GRB.BINARY, name="y_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax)]
		z = [[m.addVar(vtype=GRB.BINARY, name="z_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		w = [[m.addVar(vtype=GRB.BINARY, name="w_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		k = [[m.addVar(vtype=GRB.BINARY, name="k_"+str(i+1)+"_"+str(j)) for j in range(blockSize)] for i in range(rMax-1)]
		m.update()

		#Partial round
		#Constraints for Permutation
		for i in range(blockSize):
			m.addConstr(z0[expandP[i]] == y0[i])
		#Constraints for MixColumn
		if linAsSbox:
			for col in range(nbCol):
				for offset in range(sboxSize):
					inputvar = [z0[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					outputvar = [w0[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					addSboxConstr(m,ineqMC,inputvar,outputvar)
		else:
			for col in range(nbCol):
				#First, create the temporary variables for the Copy + XOR modeling
				t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
				for i in range(M_nrows):
					for j in range(M_ncols):
						if M[i][j] == 1:
							t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(0)+"_"+str(col)+"_"+str(i)+"_"+str(j))

				#Copy constraints
				for j in range(M_ncols):
					m.addGenConstrOr(z0[col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

				#XOR constraints
				for i in range(M_nrows):
					m.addConstr(w0[col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(w0[i] + k0[i] == x[0][i])


		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox (ugly management for the dedicated present model but I need a short coding time for this)
			if r < rMax-1:
				for j in range(nbSbox):
					inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
					outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
					addSboxConstr(m, ineqSbox, inputvar, outputvar)
			elif not dedicatedPRESENTlastLayer:
				#no dedicated last layer
				for j in range(nbSbox):
					inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
					outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
					addSboxConstr(m, ineqSbox, inputvar, outputvar)
			else:
				#dedicated last layer (y0,y1+y3,y2,y3)
				anfS[1] = anfS[1] + anfS[3]
				T = divPropANFBinTable(anfS)
				ineqSboxDecicated = sboxReducedInequalities(T)
				for j in range(nbSbox):
					inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
					outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
					addSboxConstr(m, ineqSboxDecicated, inputvar, outputvar)

			if r < rMax-1:
				#Constraints for Permutation
				for i in range(blockSize):
					m.addConstr(z[r][expandP[i]] == y[r][i])

				#Constraints for MixColumn
				if linAsSbox:
					for col in range(nbCol):
						for offset in range(sboxSize):
							inputvar = [z[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							outputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
							addSboxConstr(m,ineqMC,inputvar,outputvar)
				else:
					for col in range(nbCol):
						#First, create the temporary variables for the Copy + XOR modeling
						t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
						for i in range(M_nrows):
							for j in range(M_ncols):
								if M[i][j] == 1:
									t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r+1)+"_"+str(col)+"_"+str(i)+"_"+str(j))

						#Copy constraints
						for j in range(M_ncols):
							m.addGenConstrOr(z[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

						#XOR constraints
						for i in range(M_nrows):
							m.addConstr(w[r][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

				#Constraints for ARK
				for i in range(blockSize):
					m.addConstr(w[r][i] + k[r][i] == x[r+1][i])
	else:
		#first variables for the first (partial) round
		y0 = [m.addVar(vtype=GRB.BINARY, name="y_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		z0 = [m.addVar(vtype=GRB.BINARY, name="z_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		w0 = [m.addVar(vtype=GRB.BINARY, name="w_"+str(0)+"_"+str(j)) for j in range(blockSize)]
		k0 = [m.addVar(vtype=GRB.BINARY, name="k_"+str(0)+"_"+str(j)) for j in range(blockSize)]

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
			for col in range(nbCol):
				for offset in range(sboxSize):
					inputvar = [z0[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					outputvar = [w0[col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
					addSboxConstr(m,ineqMC,inputvar,outputvar)
		else:
			for col in range(nbCol):
				#First, create the temporary variables for the Copy + XOR modeling
				t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
				for i in range(M_nrows):
					for j in range(M_ncols):
						if M[i][j] == 1:
							t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(0)+"_"+str(col)+"_"+str(i)+"_"+str(j))

				#Copy constraints
				for j in range(M_ncols):
					m.addGenConstrOr(z0[col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

				#XOR constraints
				for i in range(M_nrows):
					m.addConstr(w0[col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

		#Constraints for ARK
		for i in range(blockSize):
			m.addConstr(w0[i] + k0[i] == x[0][i])

		#Round Function constraints
		for r in range(rMax):
			#Constraints for Sbox
			for j in range(nbSbox):
				inputvar = [x[r][sboxSize*j+i] for i in range(sboxSize)]
				outputvar = [y[r][sboxSize*j+i] for i in range(sboxSize)]
				addSboxConstr(m, ineqSbox, inputvar, outputvar)

			#Constraints for Permutation
			for i in range(blockSize):
				m.addConstr(z[r][expandP[i]] == y[r][i])

			#Constraints for MixColumn
			if linAsSbox:
				for col in range(nbCol):
					for offset in range(sboxSize):
						inputvar = [z[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						outputvar = [w[r][col*colSize + offset + j*sboxSize] for j in range(M_ncols)]
						addSboxConstr(m,ineqMC,inputvar,outputvar)
			else:
				for col in range(nbCol):
					#First, create the temporary variables for the Copy + XOR modeling
					t = [[0 for __ in range(M_ncols)] for _ in range(M_nrows)]
					for i in range(M_nrows):
						for j in range(M_ncols):
							if M[i][j] == 1:
								t[i][j] = m.addVar(vtype=GRB.BINARY, name="t_"+str(r+1)+"_"+str(col)+"_"+str(i)+"_"+str(j))

					#Copy constraints
					for j in range(M_ncols):
						m.addGenConstrOr(z[r][col*M_ncols + j], [t[i][j] for i in range(M_nrows) if M[i][j] == 1])

					#XOR constraints
					for i in range(M_nrows):
						m.addConstr(w[r][col*M_ncols + i] == quicksum(t[i][j] for j in range(M_ncols) if M[i][j] == 1))

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