import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, x, wG, wB, k, N, low, kappa, v_GB, v_BG = sy.symbols('p pp q s a x_L x_S x_LL x_LS x_SS x wG wB k N low kappa v_GB v_BG')

B = ((1 - q) * wG * v_GB - q * wB * v_BG) * (wB - wG)
C = wG * wB * (v_GB + v_BG)
V = wG * v_GB + wB * v_BG

sigma = sy.Sum(((-V) ** k ) * sy.binomial(N, k + 1) * x ** k, (k, 0, N - 1))
fitness = - (a * s * x + 1 - a) * (C * N / V - (B / V) * sigma) + N



subsses = [
    (wB, sy.Rational(4, 5)),
    (q, sy.Rational(3, 10)),
    (wG, sy.Rational(1, 5)),
    (v_GB, sy.Rational(3, 10)),
    (v_BG, sy.Rational(7, 10))
]

sulm_round = fitness.subs(subsses)
print(sulm_round)

fig = plt.figure(figsize=(12, 10))
ax = []

xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(3, 5)
Nv = [2, 3, 10, 70, 1000]
ax.append(fig.add_subplot(2, 3, 1))   
ax[0].set_title("$N=1$")
ax[0].set_xlabel("$a$")
ax[0].set_ylabel("$s$")
ax[0].set_xlim(0, 0.8)
ax[0].set_ylim(0, 1)
ax[0].axvspan(0, 1, color="b", alpha=0.4)


for l in range(5):
    ax.append(fig.add_subplot(2, 3, l + 2))     
    sulm_roundN = sulm_round.subs(N, Nv[l]).doit()
    ax[l+1].set_title("$N=$" + str(Nv[l]))
    ax[l+1].set_xlabel("$a$")
    ax[l+1].set_ylabel("$s$")
    ax[l+1].set_xlim(0, 0.8)
    ax[l+1].set_ylim(0, 1)
    LL = sulm_roundN.subs(x, xLL)
    LS = sulm_roundN.subs(x, xLS)
    SS = sulm_roundN.subs(x, xSS)
    SSLS = sy.solve(SS - LS, s)[0]
    LSLL = sy.solve(LS - LL, s)[0]
    SSLL = sy.solve(SS - LL , s)[0]
    aa = np.array([])
    SSLSbox = np.array([])
    LSLLbox = np.array([])
    SSLLbox = np.array([])
    for a_jj in range(0, 81):
        a_j = sy.Rational(a_jj, 100)
        aa = np.append(aa, float(a_j))
        SSLS_s = SSLS.subs(a, a_j)
        LSLL_s = LSLL.subs(a, a_j)
        SSLL_s = SSLL.subs(a, a_j)
        SSLSbox = np.append(SSLSbox, float(sy.re(SSLS_s)))
        LSLLbox = np.append(LSLLbox, float(sy.re(LSLL_s)))
        SSLLbox = np.append(SSLLbox, float(sy.re(SSLL_s))) 
        print(SSLSbox)

    ax[l+1].fill_between(aa, LSLLbox, 1, color="b", where=LSLLbox > SSLSbox, alpha=0.4)

    ax[l+1].fill_between(aa, SSLSbox, 1, color="b" , where=LSLLbox < SSLSbox, alpha=0.4)

    ax[l+1].fill_between(aa, SSLSbox, LSLLbox, color="w", where = SSLSbox > LSLLbox, hatch="//",alpha=0.6, edgecolor="r")

    ax[l+1].fill_between(aa, SSLSbox, LSLLbox, color="w", where = SSLSbox > LSLLbox, hatch="\\",alpha=0.6, edgecolor="b")

    ax[l+1].fill_between(aa, SSLSbox, LSLLbox, color="y", edgecolor="0.2", where = SSLSbox < LSLLbox, hatch="..", alpha = 0.6)

    ax[l+1].fill_between(aa, 0, SSLSbox, color="r", where = LSLLbox > SSLSbox , alpha = 0.4)

    ax[l+1].fill_between(aa, 0, LSLLbox, color="r", where = LSLLbox < SSLSbox , alpha = 0.4)
 
plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__))+"/Figure3.pdf", format="pdf")
plt.show()
