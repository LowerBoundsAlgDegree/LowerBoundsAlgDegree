from sage.crypto.boolean_function import BooleanFunction
import sys

def divPropANFBinTable(ANF):
	"""
	Given the ANF of the sbox
	Return a dict T such that T[(x0,x1,x2,x3)] contains all possible transitions (x0,x1,x2,x3) -> y
	e.g. T[(1,0,0,1)] = [(0,0,0,0),(0,1,0,0),(1,1,0,0)]
	Essentially the same as function divPropANF, but put directly into a binary table format
	"""

	BPR = ANF[0].parent()
	x = BPR.gens()
	n = len(x)
	m = len(ANF)
	F2n = vectorF2n(n)

	#xu -> u map
	mapMonomial = dict()
	for u in F2n:
		xu = prod(x[i] if u[i] != 0 else BPR(1) for i in range(n))
		mapMonomial[xu] = u

	#Prepare the T map
	T = dict()
	for u in F2n:
		T[u] = []

	for v in F2n: #Loop on v first as the product y^v is costly

		yv = prod(ANF[i] if v[i] != 0 else BPR(1) for i in range(m))
		yvMonomials = yv.monomials()
		nbmon = len(yvMonomials)
		
		for xu in yvMonomials:
			T[mapMonomial[xu]].append(v)

	return T

def anfSbox(SBOX):
	"""
	Given a list SBOX of 2^n value
	return (P, y) where
	y is a list of n ANF for each output coordinate
	P is a common BooleanPolynomialRing for all element of y
	"""

	#size of the sbox
	n = max([x.nbits() for x in SBOX])

	#Get the bit representation of each value
	SBOX_bits = [x.bits() + [0 for i in range(n-len(x.bits()))] for x in SBOX]

	#Create the boolean function corresponding to each bit of the output from its truth table
	#i.e. all bits of index i in SBOX-table for the i-th output bit
	B = [ BooleanFunction([x[i] for x in SBOX_bits]) for i in range(n)]

	#Set a common BooleanPolynomialRing for each output function
	y0 = B[0].algebraic_normal_form()
	P = y0.ring()

	#Compute the ANF of each output bit
	y = [P(b.algebraic_normal_form()) for b in B]

	return (P, y)

def vectorF2n(n):
	"""
	Create the list of all vector in {0,1}^n
	"""
	return [tuple(Integer(c).bits() + [0 for i in range(n-Integer(c).nbits())]) for c in range(2^n)]

def vecToInt(x):
	vx = 0
	for i in range(len(x)):
		if x[i] == 1:
			vx += (1 << i)

	return vx

def greater(a,b):
	#return True if a[i] >= b[i] for all i
	#False otherwise
	for i in range(len(a)):
		if a[i] < b[i]:
			return False
	return True

#Pretty printing stuff
def cmpvec(x,y):
	"""
	Relation between two boolean vectors to get a nice print
	"""
	if x == y:
		return int(0)
	elif sum(x) < sum(y):
		return int(-1)
	elif sum(x) == sum(y) and x > y:
		return int(-1)
	else:
		return int(1)

def printDivTable(D):
	"""
	Print the table when in binary format
	"""
	Lx = D.keys()
	Lx.sort(cmpvec)

	for dx in Lx:
		s = ""
		for ex in dx:
			s+=str(ex)
		s+= " : "
		Ly = D[dx]
		Ly.sort(cmpvec)
		for dy in Ly:
			for ey in dy:
				s+= str(ey)
			s += " "
		print(s)

def printStringTable(M):
	"""
	Print the table when in string format
	"""
	Mkeys = M.keys()
	Mkeys = [[int(i) for i in list(m)] for m in Mkeys]
	Mkeys.sort(cmpvec)

	for dx in Mkeys:
		s = ""
		for ex in dx:
			s+=str(ex)
		s+= " : "
		Ly = M["".join(str(tmp) for tmp in dx)]
		Ly = [[int(i) for i in list(m)] for m in Ly]
		Ly.sort(cmpvec)
		for dy in Ly:
			for ey in dy:
				s+= str(ey)
			s += " "
		print(s)

def expandMatrix(M, modulus=[1,1,0,1,1,0,0,0,1]):
	"""
	Expand the matrix M given on GF(2)[X]/modulus over GF(2)
	M is given as a matrix of integer, no need to coerce it to the actual field
	The modulus is given as a list of coefficient, in ascending degree order
	(i.e. what would be obtained with p.coefficients(sparse=False) where p is a Sage Polynomial)
	e.g : p = x^8 + x^4 + x^3 + x + 1 <=> modulus = [1, 1, 0, 1, 1, 0, 0, 0, 1]
	"""

	deg_mod = len(modulus)-1

	#Precomptuation of the matrices for multiplication by x^i
	#Matrix for multiplication by x
	Mx = Matrix(GF(2),deg_mod, deg_mod)
	#First subdiagonal set to 1
	for j in range(1,deg_mod):
		Mx[j,j-1] = 1
	#Last column set to the coefficients of the modulus
	for j in range(deg_mod):
		Mx[j,deg_mod-1] = modulus[j]

	#Compute each Mxi = matrix of the multiplication by x^i
	L_Mxi = [0 for i in range(deg_mod)]
	L_Mxi[0] = identity_matrix(GF(2),deg_mod)
	for i in range(1,deg_mod):
		L_Mxi[i] = Mx*L_Mxi[i-1]



	#Precomputation for the matrix conversion
	# multMatrix[a] will contains the matrix representing the multiplication of an element by a € GF(2^m) at a bit level
	multMatrix = dict()
	#Get the set of all coefficient appearing in M
	set_coeff = set([int(M[i][j]) for j in range(M.ncols()) for i in range(M.nrows())])
	#Compute each multiplication matrix we will need
	for a in set_coeff:
		Ma = matrix(GF(2),deg_mod,deg_mod)
		for b in range(deg_mod):
			if (a >> b) & 1:
				Ma += L_Mxi[b]
		multMatrix[a] = Ma

	#Now we can convert the diffusion matrix M on GF(2^m) to a matrix Mbit on GF2
	Mbit = matrix(GF(2), 0, M.nrows()*deg_mod)
	for i in range(M.nrows()):
		Mrow = matrix(GF(2), deg_mod, 0)
		for j in range(M.ncols()):
			Mrow = Mrow.augment(multMatrix[M[i][j]])
		Mbit = Mbit.stack(Mrow)

	return Mbit