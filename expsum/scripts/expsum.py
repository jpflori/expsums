def hexprint(a):
    return '0x'+hex(a)

def listprint(L):
    r = ""
    for a in L[:-1]:
        r += myprint(a)
        r += ", "
    r += myprint(L[-1])
    return r

def myprint(a):
    if isinstance(a, list):
        return listprint(a)
    elif a in ZZ:
        return hexprint(ZZ(a))
    else:
        return str(a)

def dictprint(replacements):
    for (name, value) in replacements.items():
        print name, value

def myfilter(name, value, types):
    return not (name == 'allvars'
                or name.startswith('_')
                or inspect.isfunction(value)
                or inspect.isclass(value)
                or inspect.ismodule(value)
                or isinstance(value, types))

def write_header(replacements):
    infile = open('./src/expsum.h.in', 'r')
    outfile = open('./src/expsum.h', 'w')

    for line in infile:
        for src, target in replacements.iteritems():
            line = line.replace(src, target)
        outfile.write(line)

    infile.close()
    outfile.close()

import inspect
import sage.all
from sage.rings.integer_ring import ZZ
from sage.rings.finite_rings.constructor import GF
from sage.rings.polynomial.polynomial_ring import polygen
from sage.categories.homset import Hom
from sage.rings.arith import euler_phi

MAX_THREADS=ZZ(input("max_threads = "))

modulus = input("modulus (1->minimal_weight, 2->conway, other->str(other)) = ")
if modulus == 1:
    modulus = 'minimal_weight'
elif modulus == 2:
    modulus = 'conway'
else:
    modulus = str(modulus)

m2 = ZZ(input("m2 = "))
m1 = 2*m2
m0 = 2*m1
M = [m0, m1, m2]
M0BS, M0BSL = m0.quo_rem(8)
POWM1 = 2**m1
POWM1D2 = 2**(m1-1)
POWM2 = 2**m2
POWM2T2 = 2**(m2+1)
POWM2M2D3 = ZZ((2**m2-2)/3)
POWM2P1D3 = ZZ((2**m2+1)/3)
K1MASK = 2**m1 - 1
K1HIGHBIT = 2**(m1-1)
K0MASK = 2**m0 - 1
K0HIGHBIT = 2**(m0-1)
NECKLACES = 1/m1*sum(euler_phi(d)*2**(ZZ(m1/d)) for d in m1.divisors())
if m2 > 8:
    USE_DWORD = 1
else:
    USE_DWORD = 0

GF2 = GF(2)
x = polygen(GF2)
K2 = GF(2**m2, name='g2', modulus=modulus)
g2 = K2.gen()
z2 = K2.multiplicative_generator()
K1 = GF(2**m1, name='g1', modulus=modulus)
g1 = K1.gen()
z1 = K1.multiplicative_generator()
Z1 = z1.polynomial().change_ring(ZZ)(2)
K0 = GF(2**m0, name='g0', modulus=modulus)
g0 = K0.gen()
z0 = K0.multiplicative_generator()
L = [K0, K1, K2]
moduli = [k.modulus() for k in L]
F = [f.change_ring(ZZ)(2) for f in moduli]
weights = [f.hamming_weight()-2 for f in moduli]
degrees = [(f-x**f.degree()-1).exponents() for f in moduli]
FWEIGHTS = weights
F0D = degrees[0]
if F0D == [1]:
    F0X = 1
else:
    F0X = 0
F1D = degrees[1]
if F1D == [1]:
    F1X = 1
else:
    F1X = 0
if modulus == 'minimal_weight':
    SPARSE_MODULI = 1
else:
    SPARSE_MODULI = 0

k1hp = [g1**i for i in xrange(m1, 2*m1)]
K1HP = [z.polynomial().change_ring(ZZ)(2) for z in k1hp]
k0hp = [g0**i for i in xrange(m0, 2*m0)]
K0HP = [z.polynomial().change_ring(ZZ)(2) for z in k0hp]

t1 = z0**(2**m1-1)
t1inv = t1**-1
t1pow = [t1**(2*i) for i in xrange(MAX_THREADS + 1)]
t1invpow = [t1**(-2*i) for i in xrange(MAX_THREADS + 1)]
T1 = [z.polynomial().change_ring(ZZ)(2) for z in t1pow]
T1INV = [z.polynomial().change_ring(ZZ)(2) for z in t1invpow]

t2 = z1**(2**m2-1)
t2pow = [K1(1), t2, t2**2]
T2 = [z.polynomial().change_ring(ZZ)(2) for z in t2pow]

zeta3 = z1**((2**m1-1)/3)
ZETA3 = zeta3.polynomial().change_ring(ZZ)(2)

H = Hom(K1, K0)
g1im = [h(g1) for h in H]
weights = [z.polynomial().hamming_weight() for z in g1im]
index = weights.index(min(weights))
k0g1 = g1im[index]
h = H[index]
powers = [k0g1**i for i in xrange(m1)]
K0G1 = [z.polynomial().change_ring(ZZ)(2) for z in powers]
K0T2 = [h(z).polynomial().change_ring(ZZ)(2) for z in t2pow]
K0ZETA3 = h(zeta3).polynomial().change_ring(ZZ)(2)

TR = [(g0**i).trace() for i in xrange(m0)]
K0TRACES = ZZ(sum(2**i*tr.lift() for i,tr in enumerate(TR)))

GF4TR = [sum((g1**i)**(2**(2*j)) for j in xrange(m2)) for i in xrange(m1)]
K1GF4TRACES = [ZZ(sum(2**i for i,tr in enumerate(GF4TR) if tr == zeta3**j)) for j in xrange(3)]

allvars = vars().items()
replacements = dict(('%'+name+'%', myprint(value)) for (name, value) in allvars
                    if myfilter(name, value, type(GF)))

if __name__ == "__main__":
    write_header(replacements)
