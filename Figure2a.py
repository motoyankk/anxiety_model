import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

sy.var("p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, v_GB, v_BG, x, wG, wB, N")
def nextqd(N):
    return (q - v_GB/(v_BG + v_GB)) * (1 - v_BG - v_GB) ** (N - 1) + v_GB/(v_BG + v_GB)

fig = plt.figure()
substitutions = [
    (wB, sy.Rational(9, 10)),
    (q, sy.Rational(4, 5)), 
    (wG, sy.Rational(1, 20)), 
    (a, sy.Rational(4, 10)), 
    (s, sy.Rational(1, 2)),
    (v_BG, sy.Rational(1, 10)),
    (v_GB,  sy.Rational(1, 2))]

ax = fig.add_subplot(1, 1, 1)
v_BG_l =  sy.Rational(1, 2)
v_GB_m = sy.Rational(1, 10)
ax.set_xlabel("$N$")
ax.set_ylabel("$q(t)$")
NN = []
DD=[]
SSbox = []
LSbox = []
LLbox = []
for N_jj in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]:
    N_j = N_jj
    NN.append(N_j)

    sud_round = nextqd(N_j).subs(substitutions)
    sulmd_round = sud_round.subs([(v_BG, v_BG_l), (v_GB, v_GB_m)])
    DD.append(float(sulmd_round))
    print(N_jj)

ax.plot(NN, DD, "o", color='k', markersize=4, markeredgewidth=1, markeredgecolor='k', alpha=0.8,label='Dormancy')
plt.xscale('log')
plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure2a.pdf", format="pdf")
plt.show()
