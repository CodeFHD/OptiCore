def f_lensmaker(r1,r2,n, d):
    dummy = (n-1) * (1./r1 - 1./r2 + (n-1)*d/(n*r1*r2))
    return 1/dummy



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