import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, v_GB, v_BG, x, wG, wB, k, N, low, kappa, ip = sy.symbols('p pp q s a x_L x_S x_LL x_LS x_SS v_GB v_BG x wG wB k N low kappa ip')


B = ((1 - q) * wG * v_GB - q * wB * v_BG) * (wB - wG)
C = wG * wB * (v_GB + v_BG)
V = wG * v_GB + wB * v_BG
hx = low * (x / sy.Rational(4, 5)) ** kappa
sigma = sy.summation(((-V) ** k ) * sy.binomial(N, k + 1) * x ** k, (k, 0, N - 1))
fitness = 1 + (- (a * s * (x + hx) + 1 - a) * (C * N / V - (B / V) * sigma) + N)*ip

subsses = [
    (wB, sy.Rational(4, 5)),
    (q, sy.Rational(1, 2)),
    (wG, sy.Rational(1, 5)),
    (v_GB, sy.Rational(1, 10)),
    (v_BG, sy.Rational(1, 2)),
    (N, 30),
    (ip, sy.Rational(1, 100))
]

sulm_round = fitness.subs(subsses)
fig = plt.figure(figsize=(5, 10))
ax = []

xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(4, 5)
ax2 = fig.add_subplot(4, 2, 1)
ax2.set_xlabel("$x$")
ax2.set_ylabel("$h(x)$")
box = [np.array([]) for _ in range(7)] 
linebox=[":", "--", "-.", "-", "--", "-.", "-"]
colorbox=["k", "0.3", "0.3","0.3","0.6","0.6","0.6"]
#colorbox=["black", "b", "r","y","pink","green","purple"]
low_l = [0, sy.Rational(1, 10), sy.Rational(1, 5), sy.Rational(1, 10), sy.Rational(1, 5), sy.Rational(1, 10), sy.Rational(1, 5)]
kappa_l = [0, sy.Rational(1, 30), sy.Rational(1, 30), 1, 1, 20, 20]

for l in range(7):
    ss = hx.subs([(low, low_l[l]), (kappa, kappa_l[l])])
    xx=np.array([])
    for x_jj in range(81):
        x_j = sy.Rational(x_jj, 100)
        xx = np.append(xx, float(x_j))        
        box[l] = np.append(box[l], float(ss.subs(x, x_j)))
    ax2.plot(xx, box[l], linestyle = linebox[l], color = colorbox[l])
ll=[2, 3, 4, 5, 6, 7, 8]
for l in range(7):
    ax.append(fig.add_subplot(4, 2, ll[l]))     
    sulm_roundN = sulm_round.subs([(low, low_l[l]), (kappa, kappa_l[l])])
    ax[l].set_title("$p=$" + str(low_l[l])+"   "+"$K=$" + str(kappa_l[l]))
    ax[l].set_xlabel("$a$")
    ax[l].set_ylabel("$s$")
    ax[l].set_xlim(0, 0.8)
    ax[l].set_ylim(0, 1)
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
    for a_jj in range(501):
        a_j = sy.Rational(a_jj, 500)
        aa = np.append(aa, float(a_j))  
        SSLS_s = SSLS.subs(a, a_j)
        LSLL_s = LSLL.subs(a, a_j)
        SSLL_s = SSLL.subs(a, a_j)
        SSLSbox = np.append(SSLSbox, float(sy.re(SSLS_s)))
        LSLLbox = np.append(LSLLbox, float(sy.re(LSLL_s)))
        SSLLbox = np.append(SSLLbox, float(sy.re(SSLL_s))) 
        #print(SSLSbox)
    ax[l].set_ylim(0, 1)


    ax[l].fill_between(aa, LSLLbox, 1, color="b", where=LSLLbox > SSLSbox, alpha = 0.4)

    ax[l].fill_between(aa, SSLSbox, 1, color="b" , where=LSLLbox < SSLSbox, alpha = 0.4)

    ax[l].fill_between(aa, SSLSbox, LSLLbox, color="k", where = SSLSbox > LSLLbox, alpha = 0.6)

    ax[l].fill_between(aa, SSLSbox, LSLLbox, color="y", edgecolor="0.2", where = SSLSbox < LSLLbox, hatch="..", alpha = 0.6)

    ax[l].fill_between(aa, 0, SSLSbox, color="r", where = LSLLbox > SSLSbox , alpha = 0.4)

    ax[l].fill_between(aa, 0, LSLLbox, color="r", where = LSLLbox < SSLSbox , alpha = 0.4)
plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__))+"/Figure5.pdf", format="pdf")
plt.show()
