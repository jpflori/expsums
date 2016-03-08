def traceL(s):
    return sum([s**(2**i) for i in xrange(m1)])

def halftraceL(s):
    return sum([s**(2**(2*i)) for i in xrange(m2)])

def KS(a):
    r"""
    Kloosterman sum
    """
    S = 1
    t = K(1)
    for i in xrange(2**m1-1):
        t *= y
        s = a*t + t**(2**m1-2)
        S += (-1)**(ZZ(traceL(s)))
    return S

def KSpartial(a):
    r"""
    Partial Kloosterman sums
    """
    S = {}
    S[0] = 1
    S[beta] = 0
    S[beta**2] = 0
    S[beta**3] = 0
    t = K(1)
    for i in xrange(2**m1-1):
        t *= y
        s = a*t + t**(2**m1-2)
        r = t**((2**m1-1)/3)
        S[r] += (-1)**(ZZ(traceL(s)))
    return S

def KSpartialw1(a, w1):
    r"""
    Twisted partial Kloosterman sums shifted by w1**-1
    """
    S = {}
    S[K(0)] = 0
    S[beta] = 0
    S[beta**2] = 0
    S[beta**3] = 0
    tw1 = w1**-1
    tt = K(1)
    for j in xrange(2**m1):
        tw1 *= t1
        tt *= t1
        r = (tt+tt**-1)**((2**n-1)/3)
        s = (tw1+tw1**-1)**2
        s = a*s + s**(2**m1-2)
        S[r] += (-1)**(ZZ(sum([s**(2**i) for i in xrange(m1)])))
    return S

def KSpartialw1_inverse(a, w1):
    r"""
    Twisted partial Kloosterman sums shifted by w1**-1
    """
    S = {}
    R = {}
    S[K(0)] = 0
    S[beta] = 0
    S[beta**2] = 0
    S[beta**3] = 0
    R[K(0)] = 0
    R[beta] = 0
    R[beta**2] = 0
    R[beta**3] = 0
    t = z**(2**m1+1)
    tw1 = w1**-1
    tt = K(1)
    for j in xrange(1,2**m1-1):
        tt *= t
        tw1 *= t
        r = (tw1**-1+tw1)**((2**n-1)/3)
        s = (tt+tt**-1)**2
        ss = a*s
        R[r] += (-1)**(ZZ(sum([ss**(2**i) for i in xrange(m1)])))
        s = ss + s**(2**m1-2)
        S[r] += (-1)**(ZZ(sum([s**(2**i) for i in xrange(m1)])))
    tw1 = w1**-1
    tt = K(1)
    for j in xrange(1,2**m1+1):
        tw1 *= t1
        tt *= t1
        r = (tt+tt**-1)**((2**n-1)/3)
        s = (tw1+tw1**-1)**2
        ss = a*s
        R[r] += (-1)**(ZZ(sum([ss**(2**i) for i in xrange(m1)])))
        s = ss + s**(2**m1-2)
        S[r] += (-1)**(ZZ(sum([s**(2**i) for i in xrange(m1)])))
    return S, R

def cubic(a):
    r"""
    Cubic sum
    """
    S = 1
    t = K(1)
    for i in xrange(2**m1-1):
        t *= y
        s = a*t**3
        S += (-1)**(ZZ(traceL(s)))
    return S

def char(a):
    r"""
    Non-principal cubic character
    """
    return a**((2**n-1)/3)

def chi(a, w):
    r"""
    (-1)^(Tr^n_1(a w^-1))
    """
    return (-1)**(ZZ((a*w**(1-2**m1)).trace()))

def chif(a, b, w):
    r"""
    Walsh--Hadamard transform at w of f_{a, b}
    """
    S = 0
    for t in K:
        s = b*t**((2**n-1)/3)
        S += (-1)**(ZZ((a*t**(2**m1-1)).trace() + s + s**2 + (w*t).trace()))
    return S

def S12(a, b):
    r"""
    Sum on U_1 and U_2 (for w = 0)
    """
    u1 = K(1)
    u2 = K(1)
    S1 = 0
    S2 = 0
    for i in xrange(2**m1+1):
        u1 *= t1
        s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
        S1 += s1
    for j in xrange(2**m2+1):
        u2 *= t2
        s2 = b*u2**((2**n-1)/3)
        S2 += (-1)**(ZZ(s2 + s2**2))
    return S1*S2

def S12w(a, b, w):
    r"""
    Sum on (u_1, u_2) \in U_1 x U_2 such that the trace is 0
    """
    u1 = K(1)
    u2 = K(1)
    S = 0
    for i in xrange(2**m1+1):
        u1 *= t1
        s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
        for j in xrange(2**m2+1):
            u2 *= t2
            s = w*u1*u2
            if s + s**(2**m2) + s**(2**m1) + s**(2**(m1+m2)) == 0:
                s2 = b*u2**((2**n-1)/3)
                S += (-1)**(ZZ(s2 + s2**2))*s1
    return S

def S12w1eq(a, b, w1):
    r"""
    Sum on (u_1, u_2) \in U_1 x U_2 such that the trace is 0 with u_1 = w_1^{-1}
    """
    u1 = w1**-1
    u2 = K(1)
    S = 0
    s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
    for j in xrange(2**m2+1):
        u2 *= t2
        s2 = b*u2**((2**n-1)/3)
        S += (-1)**(ZZ(s2 + s2**2))*s1
    return S

def S12w1diffstupid(a, b, w1, w2):
    r"""
    Sum on (u_1, u_2) \in U_1 x U_2 such that the trace is 0 with u_1 \neq w_1^{-1}
    Testing all possible u_2
    """
    u1 = w1**-1
    u2 = K(1)
    S = 0
    for i in xrange(2**m1):
        u1 *= t1
        s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
        for j in xrange(2**m2+1):
            u2 *= t2
            s = w1*w2*u1*u2
            if s + s**(2**m2) + s**(2**m1) + s**(2**(m1+m2)) == 0:
                s2 = b*u2**((2**n-1)/3)
                S += (-1)**(ZZ(s2 + s2**2))*s1
    return S

def S12w1diff(a, b, w1, w2):
    r"""
    Sum on (u_1, u_2) \in U_1 x U_2 such that the trace is 0 with u_1 \neq w_1^{-1}
    Only picking the correct value for u_2
    """
    u1 = w1**-1
    u2 = K(1)
    S = 0
    for i in xrange(2**m1):
        u1 *= t1
        s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
        u2 = u1*w1
        u2 = (u2 + u2**-1)**-1*w2**-1
        s2 = b*u2**((2**n-1)/3)
        S += (-1)**(ZZ(s2 + s2**2))*s1
    return S

def S12U1(a, b, w1):
    r"""
    Sum on U_1, not including w_1^{-1}, whatever the corresponding u_2
    """
    u1 = w1**-1
    S = 0
    for i in xrange(2**m1):
        u1 *= t1
        s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
        S += s1
    return S

def S12U1left(a, b, w1, w2):
    r"""
    Superfluous sum on U_1
    """
    u1 = w1**-1
    u2 = K(1)
    S = 0
    for i in xrange(2**m1):
        u1 *= t1
        u2 = u1*w1
        u2 = (u2 + u2**-1)*w2
        s2 = b*u2**(-(2**n-1)/3)
        if s2 == 1:
            s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))
            S += s1
    return S

def S12U1leftdistribution(a, b, w2):
    r"""
    Distribution of superfluous sums on U_1
    """
    D = {}
    w1 = K(1)
    for i in xrange(2^(m1-1)):
        w1 *= t1
        s = S12U1left(a, b, w1, w2)
        try:
            D[s] += 1
        except KeyError:
            D[s] = 1
    return D

def S12U1leftKS():
    r"""
    Values of superfluous sums for a given Kloosterman sum
    """
    DD = {}
    B = 2**(m2+1)
    k = 1 - B
    while ((k-4) % 12) != 0:
        k += 1
    while k < 1 + B:
        DD[k] = set()
        k += 12
    a = K(1)
    b = K(1)
    for i in xrange((2**m1-1)/3):
        b *= y
        a *= y**3
        if traceL(b) != 0:
            k = KS(a)
            l = S12U1left(a, K(1), K(1), K(1))
            D = S12U1leftdistribution(a, K(1), K(1))
            DD[k].add(D[l])
    return DD

def S12U1leftabs(a, b, w1, w2):
    r"""
    Superfluous sum of gathered absolute values on U_1
    """
    u1 = w1**-1
    uu1 = w1
    u2 = K(1)
    S = 0
    for i in xrange(2**(m1-1)):
        u1 *= t1
        uu1 *= t1
        u2 = u1*w1
        u2 = (u2 + u2**-1)*w2
        s2 = b*u2**(-(2**n-1)/3)
        if s2 == 1:
            s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))+(-1)**(ZZ((a*(uu1)**(2**m1-1)).trace()))
            S += abs(s1)
    return S

def S12U1leftsquare(a, b, w1, w2):
    r"""
    Superfluous sum of gathered squares on U_1
    """
    u1 = w1**-1
    uu1 = w1
    u2 = K(1)
    S = 0
    for i in xrange(2**(m1-1)):
        u1 *= t1
        uu1 *= t1
        u2 = u1*w1
        u2 = (u2 + u2**-1)*w2
        s2 = b*u2**(-(2**n-1)/3)
        if s2 == 1:
            s1 = (-1)**(ZZ((a*(u1)**(2**m1-1)).trace()))+(-1)**(ZZ((a*(uu1)**(2**m1-1)).trace()))
            S += s1**2
    return S

def my_init(t):
    global m0, m1, m2, n, K, z, y, beta, t1, t2, L, a, b, w1, w2
    m2 = t
    m1 = 2*m2
    n = m0 = 2*m1
    K.<z> = GF(2**n)
    y = z^(2**m1+1)
    beta = z**((2**n-1)/3)
    t1 = z**(2**m1-1)
    t2 = z**((2**m1+1)*(2**m2-1))

    b = K(1)
    w1 = K(1)
    w2 = K(1)
    w = w1*w2

    a = K.random_element()
    a += a**(2**m1)

def test():
    while True:
        w1 = K.random_element()**(2**m1-1)
        if ((w1 + w1**-1)**((2**n-1)/3) == K(1) and (a*w1**-2).trace() == K(1)):
            print S12U1left(a,b,w1,w2)
