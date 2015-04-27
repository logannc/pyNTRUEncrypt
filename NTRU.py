from itertools import chain
import math, random

def invModQ(d,q):
	for i in xrange(1,q):
		if (d*i)%q == 1:
			return i
	return None

class ConvPoly(object):
	def __init__(self, coef=[0], N=None):
		if N is None:
			self.N = len(coef)
			self.coef = coef
		else:
			self.N = N
			self.coef = coef + [0]*(N-len(coef))
	def __repr__(self):
		return type(self).__name__ + str(self.coef)
	def __add__(self, other):
		if isinstance(other, type(self)) and (self.N == other.N):
			return ConvPoly(map(sum,zip(self.coef, other.coef)))
		elif isinstance(other, int):
			return ConvPoly([self.coef[0]+other] + self.coef[1:])
		else:
			return NotImplemented
	def __radd__(self, other):
		return self + other
	def __eq__(self, other):
		if self.coef == other.coef:
			return True
		return False
	def __ne__(self, other):
		return not (self == other)
	def __neg__(self):
		return type(self)(map(lambda x: -x, other.coef), other.q)
	def __sub__(self, other):
		return self + (-other)
	def __rsub__(self, other):
		return self - other
	def __mul__(self, other):
		if isinstance(other, type(self)) and self.N == other.N:
			coefs = []
			for k in xrange(self.N):
				s = 0
				for i in xrange(self.N):
					s += self.coef[i] * other.coef[(k-i)%self.N]
				coefs.append(s)
			return type(self)(coefs)
		elif isinstance(other, int):
			return type(self)([ other*c for c in self.coef ])
		else:
			return NotImplemented
	def __rmul__(self, other):
		return self*other

class PolyModQ(object):
	def __init__(self, coef=[0], q=3):
		coef = map(lambda x: x%q, coef)
		for index, val in enumerate(coef[::-1]):
			if val != 0:
				break
		self.coef = coef[:len(coef)-index]
		self.degree = len(self.coef)-1
		self.q = q
	def __repr__(self):
		return type(self).__name__ + "(" + str(self.coef) + ", " + str(self.q) + ")"
	def __eq__(self, other):
		if self.degree == other.degree:
			for pair in zip(self.coef, other.coef):
				if pair[0] != pair[1]:
					return False
			return True
		return False
	def __ne__(self, other):
		return not (self == other)
	def __neg__(self):
		return type(self)(map(lambda x: -x, self.coef), self.q)
	def __add__(self, other):
		if isinstance(other, type(self)) and self.q == other.q:
			if self.degree > other.degree:
				return type(self)(map(sum,zip(self.coef, other.coef + [0]*(self.degree-other.degree))),self.q)
			elif self.degree < other.degree:
				return other + self
			else:
				return type(self)(map(sum,zip(self.coef, other.coef)),self.q)
		elif isinstance(other, ConvPoly):
			return self + ConvModQ(other.coef, self.q, self.N)
		elif isinstance(other, int):
			return type(self)([self.coef[0]+other] + self.coef[1:],self.q)
		else:
			return NotImplemented
	def __radd__(self, other):
		return self + other
	def __sub__(self, other):
		if isinstance(other, type(self)) or isinstance(other, int):
			return self + (-other)
		else:
			return NotImplemented
	def __rsub__(self, other):
		return self - other
	def __mul__(self, other):
		if isinstance(other, type(self)) and self.q == other.q:
			coef = [0]*(self.degree+other.degree+1)
			for index1,c1 in enumerate(self.coef):
				for index2,c2 in enumerate(other.coef):
					coef[index1+index2] += c1*c2
			return type(self)(coef, self.q)
		elif isinstance(other, int):
			return type(self)([ other*c for c in self.coef], self.q)
		else:
			return NotImplemented
	def __rmul__(self, other):
		return self*other
	def centerLift(self):
		coefs = []
		for c in self.coef:
			if c>self.q/2.0:
				coefs.append(c-self.q)
			else:
				coefs.append(c)
		return ConvPoly(coefs)


class ConvModQ(PolyModQ):
	def __init__(self, coef, q=3, N=None):
		coef = map(lambda x: x%q, coef)
		if N is None:
			self.N = len(coef)
			self.coef = coef
		else:
			self.N = N
			self.coef = coef + [0]*(N-len(coef))
		self.degree = len(self.coef)-1
		self.q = q
	def __repr__(self):
		return type(self).__name__ + "(" + str(self.coef) + ", " + str(self.q) + ")"
	def __mul__(self, other):
		if isinstance(other, type(self)) and self.N == other.N:
			coefs = []
			for k in xrange(self.N):
				s = 0
				for i in xrange(self.N):
					s += self.coef[i] * other.coef[(k-i)%self.N]
				coefs.append(s)
			return type(self)(coefs, self.q, self.N)
		elif isinstance(other, ConvPoly):
			other = ConvModQ(other.coef, self.q, self.N)
			return self*other
		elif isinstance(other, int):
			return type(self)([ other*c for c in self.coef ], self.q, self.N)
		else:
			return NotImplemented
	def __div__(self, other):
		if isinstance(other, type(self)):
			return self*other.inverse()
		elif isinstance(other, int):
			otherinv = invModQ(other, self.q)
			if otherinv is None:
				raise Exception("{} not invertible mod {}".format(other, self.q))
			return self*otherinv
		else:
			return NotImplemented
	def modQ(self, q):
		return ConvModQ(self.coef, q, self.N)
	def inverse(self, N=None, debug=False):
		if self.q == 2048:
			self.q = 2
		FAIL = 100000
		i = 0
		if N is None:
			N = self.N
		quotients = []
		# Extended Euclidean Algorithm
		# q = b*k + r
		q = PolyModQ([-1] + [0]*(N-1) + [1], self.q)
		k = PolyModQ([0], self.q)
		b = PolyModQ(self.coef, self.q)
		r = q
		# repeat below while r!=0 for gcd/inverse
		bdinv = invModQ(b.coef[-1], self.q)
		# print r.coef[r.degree], b.coef[-1]
		if bdinv is None:
			return None
		while r.degree >= b.degree and i<FAIL:
			rcoef = r.coef
			kp = PolyModQ([0]*(r.degree-b.degree) + [rcoef[r.degree]*bdinv], self.q)
			k = k + kp
			r = r - kp*b
			i += 1
		quotients.append(k)
		if debug:
			print "{} = {}*{} + {}".format(q,b,k,r)
		while r != PolyModQ() and i<FAIL:
			q = b
			b = r
			k = PolyModQ([0]*(N+1), self.q)
			r = q
			# print r.coef[r.degree], b.coef[-1]
			bdinv = invModQ(b.coef[-1], self.q)
			if bdinv is None:
				# print r,b
				return None
			while r.degree >= b.degree and r != PolyModQ() and i<FAIL:
				rcoef = r.coef
				kp = PolyModQ([0]*(r.degree-b.degree) + [rcoef[r.degree]*bdinv], self.q)
				k = k + kp
				r = r - kp*b
				i += 1
			quotients.append(k)
			i += 1
			if debug:
				print "{} = {}*{} + {}".format(q,b,k,r)
		if i >= FAIL:
			print "Failed to generate inverse in {} steps, stopping.".format(FAIL)
			return None
		x = [ PolyModQ([0],self.q), PolyModQ([1], self.q)]
		y = [ PolyModQ([1], self.q), PolyModQ([0],self.q)]
		for index,quot in enumerate(quotients):
			x.append(quot*x[index+1]+x[index])
			y.append(quot*y[index+1]+y[index])
		if self.q == 2:
			n = 2
			self.q = 2048
			tinv = ConvModQ(x[-2].coef, self.q, N)
			while n <= 2048:
				tinv = 2*tinv - self * tinv * tinv
				n *= 2
			return tinv
		tinv = ConvModQ(x[-2].coef, self.q, N)
		tinv = self.q*tinv - self*tinv*tinv
		return 2*tinv

class NTRUParams(object):
	def __init__(self, k, choice="speed"): # NTRU Paramgen Paper page 16
		if k==112:
			if choice == "space":
				self.N, self.d_f = (401, 113)
			elif choice == "hybrid":
				self.N, self.d_f = (541, 49)
			elif choice == "speed":
				self.N, self.d_f = (659, 38)
		elif k==128:
			if choice == "space":
				self.N, self.d_f = (449, 134)
			elif choice == "hybrid":
				self.N, self.d_f = (613, 55)
			elif choice == "speed":
				self.N, self.d_f = (761, 42)
		elif k==192:
			if choice == "space":
				self.N, self.d_f = (677, 157)
			elif choice == "hybrid":
				self.N, self.d_f = (887, 81)
			elif choice == "speed":
				self.N, self.d_f = (1087, 63)
		elif k==256:
			if choice == "space":
				self.N, self.d_f = (1087, 120)
			elif choice == "hybrid":
				self.N, self.d_f = (1171, 106)
			elif choice == "speed":
				self.N, self.d_f = (1499, 79)
		else:
			raise Exception("Not implemented. :(")
		self.q = 2048
		self.p = 3
		self.d_g = self.N/3
		self.d_r = self.d_f

class NTRUKey(object):
	def __init__(self, ring=None, f=None, g=None):
		if ring is None:
			ring = NTRUParams(256, "speed")
		elif isinstance(ring, int):
			ring = NTRUParams(ring, "speed")
		self.ring = ring
		if f is None:
			self.f = self.randomTrinary(self.ring.d_f+1, self.ring.d_f)
		if g is None:
			self.g = self.randomTrinary(self.ring.d_g, self.ring.d_g)
		self.finvq = self.f.inverse()
		while self.finvq is None:
			print "finv was None. Retrying."
			self.f = self.randomTrinary(self.ring.d_f+1, self.ring.d_f)
			self.finv = self.f.inverse()
		self.finvp = ConvModQ(self.f.centerLift().coef, self.ring.p).inverse()
		while self.finvq is None or self.finvp is None:
			print "finv was None. Retrying."
			self.f = self.randomTrinary(self.ring.d_f+1, self.ring.d_f)
			self.finv = self.f.inverse()
			if self.finv is None:
				continue
			self.finvp = ConvModQ(self.f.coef, self.ring.p).inverse()
		self.h = self.finvq * self.g
	def randomTrinary(self, d1,d2):
		arr = [1]*d1 + [-1]*d2 + [0]*(self.ring.N - d1 - d2)
		random.shuffle(arr)
		return ConvModQ(arr, self.ring.q, self.ring.N)
	def publicKey(self):
		return (self.ring, self.h)

def strToBin(m):
	return "".join(format(ord(c), 'b').zfill(8) for c in m)

def binToStr(b):
	bs = []
	for i in xrange(len(b)/8):
		bs.append(chr(int("".join([ str(x) for x in b[i*8:(i+1)*8] ]),2)))
	return "".join(bs)
	
def chunk(N, iter):
	for i in xrange(int(math.ceil(len(iter)/float(N)))):
		yield iter[i*N:(i+1)*N]

def NTRUBlockEncrypt(ring, h, m):
	rvars = [1]*ring.d_r + [-1]*ring.d_r + [0]*(ring.N-2*ring.d_r)
	random.shuffle(rvars)
	r = ConvModQ(rvars, ring.q)
	c = ring.p*r*h+m
	return c

def NTRUEncrypt(ring, pub, m):
	m = map(int,list(strToBin(m)))
	if len(m) > ring.N:
		msplit = [ m for m in chunk(ring.N, m) ]
		n = len(msplit[-1])
		m = map(lambda m: ConvPoly(m, ring.N), msplit)
	else:
		n = len(m)
		m = [ConvPoly(m, ring.N)]
	menc = map(lambda m: NTRUBlockEncrypt(ring, pub, m), m)
	return (menc,n)

def NTRUBlockDecrypt(key, c):
	a = key.f*c
	aprime = a.centerLift()
	m = key.finvp*aprime
	return m.centerLift()

def NTRUDecrypt(key, cl, n):
	if len(cl)==1:
		m = NTRUBlockDecrypt(key, cl[0])
		return binToStr(m.coef[:n])
	cl = map(lambda c: NTRUBlockDecrypt(key, c), cl)
	mlist = [ c.coef for c in cl[:-1] ]
	mlist.append(cl[-1].coef[:n])
	m = list(chain.from_iterable(mlist))
	return binToStr(m)



if __name__=='__main__':
	print "Generating key"
	key = NTRUKey(NTRUParams(256,'speed'))
	print "encrypting message"
	morig = "Hello world!"
	import time
	for i in xrange(5): # warming up the JIT for PyPy
		enc = NTRUEncrypt(key.ring, key.h, morig)
		m = NTRUDecrypt(key, *enc)
	start = time.time()
	N = 100
	for i in xrange(N):
		enc = NTRUEncrypt(key.ring, key.h, morig)
		m = NTRUDecrypt(key, *enc)
	total = time.time() - start
	timeper = total/N
	print "Time per Encryption/Decryption cycle: {}".format(timeper)
	assert(m == morig)
