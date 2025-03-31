import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

sy.var("p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, v_GB, v_BG, x, wG, wB, N")
p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, x, wG, wB, k, N, u_GB, u_BG = sy.symbols('p pp q s a x_L x_S x_LL x_LS x_SS  x wG wB k N u_GB, u_BG')
Ad = (1 / (u_GB+u_BG)) * ((wG * u_BG + wB * u_GB) * N - ((1 - q) * u_GB- q * u_BG)* (wB - wG) *((1- (1-u_GB-u_BG) ** N)/(u_GB + u_BG)))
fitness = - a * Ad * s * x+N - (1 - a) * Ad
fig = plt.figure()
xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(3, 5)
substitutions = [(wB, sy.Rational(9, 10)),
                 (q, sy.Rational(4, 5)), 
                 (wG, sy.Rational(1, 20)), 
                 (a, sy.Rational(4, 10)), 
                 (s, sy.Rational(1, 2)),
                 (u_BG, sy.Rational(1, 10)),
                 (u_GB,  sy.Rational(1, 2))]
ax = fig.add_subplot(1, 1, 1)
ax.set_xlabel("$N$")
ax.set_ylabel("$Fitness$")
NN = []
SSbox = []
LSbox = []
LLbox = []
subs_round = fitness.subs(substitutions)
for N_jj in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
    N_j= N_jj
    NN.append(N_j)
    su_round = subs_round.subs(N, N_j)
    LL = su_round.subs(x, xLL)
    LS = su_round.subs(x, xLS)
    SS = su_round.subs(x, xSS)
    print(LL)
    LLbox.append(float(sy.re(LL)/N_j))
    LSbox.append(float(sy.re(LS)/N_j))
    SSbox.append(float(sy.re(SS)/N_j))
    print(N_j)
ax.plot(NN, LLbox, "o", color='none', markersize=4, markeredgewidth=1, markeredgecolor='k', alpha=0.8,label='LL')
ax.plot(NN, LSbox, "s", color='none', markersize=4, markeredgewidth=1, markeredgecolor='k', alpha=0.8,label='LS')
ax.plot(NN, SSbox, "v", color='none', markersize=4, markeredgewidth=1, markeredgecolor='k', alpha=0.8,label='SS')
plt.legend()
plt.tight_layout()
plt.xscale('log')
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure2d.pdf", format="pdf")
plt.show()
