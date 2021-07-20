def f_lensmaker(r1,r2,n, d):
    if r1 == 0:
        invr1 = 0
    else:
        invr1 = 1./r1

    if r2 == 0:
        invr2 = 0
    else:
        invr2 = 1./r2

    if r1 == 0 or r2 == 0:
        invr1r2 = 0
    else:
        invr1r2 = 1./(r1*r2)

    dummy = (n-1) * (invr1 - invr2 + (n-1)*d*invr1r2/n)

    if dummy == 0:
        invdummy = 0
    else:
        invdummy = 1./dummy
    
    return invdummy



def petzvalrad_surfaces(ns, rs):
    prods = []
    for i in range(len(ns)-1):
        nom = ns[i+1] - ns[i]
        den = ns[i]*ns[i+1]*rs[i]
        prods.append(nom/den)
    psum = ns[-1] * sum(prods)
    return -1./psum

def petzvalrad_thinlens(ns, fs):
    prods = [1/(x*y) for x,y in zip(ns, fs)]
    psum = sum(prods)
    return 1./psum