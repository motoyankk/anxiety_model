import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

# Declare symbolic variables
#sy.var("p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, v_GB, v_BG, x, wG, wB, N")
sy.var("p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, x, wG, wB, N, v_GB, v_BG")
#v_GB = q
#v_BG= 1-q

# Define the functions 
def nextq(v_GB, x, N):
    #xLL=0の時
    #return (rd) ** n-1 *q +(1- (rd)** n-1) * Q
    #xLL>0の時
    return  (1 - (wB * v_BG + wG * v_BG) * x) ** (N-1) *q +(1- (1 - (wB * v_BG +wG * v_GB) * x ) ** (N - 1) ) * (v_GB/(v_GB + (wB / wG)* v_BG))

xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(4, 5)
substitutions = [(wB, sy.Rational(9, 10)),
                 (q, sy.Rational(4, 5)), 
                 (wG, sy.Rational(1, 2)), 
                 (a, sy.Rational(4, 10)), 
                 (s, sy.Rational(35, 100)),
                 (v_BG, sy.Rational(4, 5)),
                 (v_GB,  sy.Rational(1, 5))]

fig = plt.figure(facecolor='none') 
ax = fig.add_subplot(1, 1, 1)
ax.set_facecolor('none')
NN = []
SSbox = []
LSbox = []
LLbox = []
subs_round = nextq(v_GB, x, N).subs(substitutions)
LL_color = 'c'
LS_color = 'yellow'
SS_color = 'red'

for N_jj in range(1, 41):
    N_j = N_jj
    NN.append(N_j)
    su_round = subs_round.subs(N, N_j).doit()
    LL = su_round.subs(x, xLL)
    LS = su_round.subs(x, xLS)
    SS = su_round.subs(x, xSS)
    LLbox.append(float(sy.re(LL)/N_j))
    LSbox.append(float(sy.re(LS)/N_j))
    SSbox.append(float(sy.re(SS)/N_j))
    print(N_j)

ax.axvspan(1, 2, color=LL_color, alpha=0.4)  
ax.axvspan(3, 32, color=LS_color, alpha=0.1) 
ax.axvspan(33, 40, color=LL_color, alpha=0.4)  

# データポイントのプロット
ax.plot(NN, LLbox, "o", color=LL_color, markeredgecolor='black', markeredgewidth=1.5, 
        markersize=4, alpha=0.8, label='LL')
ax.plot(NN, LSbox, "s", color=LS_color, markeredgecolor='black', markeredgewidth=1.5, 
        markersize=4, alpha=0.8, label='LS')
ax.plot(NN, SSbox, "v", color=SS_color, markeredgecolor='black', markeredgewidth=1.5, 
        markersize=4, alpha=0.8, label='SS')
plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure2c.pdf", 
            format="pdf", 
            transparent=True,
            bbox_inches='tight')
plt.show()