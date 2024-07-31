from math import log10, sqrt
from scipy import optimize

def calc_Reynolds(rho, U, D, mu):
    return rho*U*D/mu

def atrito_Colebrook_func(f, e, D, Re):
    termo1 = e/(3.7*D)
    termo2 = 2.51/(Re*sqrt(f))
    func = 1/sqrt(f) + 2*log10(termo1 + termo2)
    return func

def atrito_Colebrook(e, D, Re):
    return optimize.root_scalar(atrito_Colebrook_func,
                                args=(e, D, Re),
                                method='secant',
                                x0=0.01,
                                x1=0.00001)


if __name__ == "__main__":

    #REP*** 30m conditions using HCs
    L = 30
    U = 20.14
    D = 0.320675*2
    rho = 2.1314776275527
    mu = 9.71207824840234E-6

    Re = calc_Reynolds(rho, U, D, mu)
    f = atrito_Colebrook(0, D, Re)

    hl = f.root*L/D*(U**2)/2
    dP = rho*hl

    dP_h = rho*9.8067*L
