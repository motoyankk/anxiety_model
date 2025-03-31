import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os


sy.var("p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, v_GB, v_BG, x, wG, wB, N")


p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, x, wG, wB, k, N, low, kappa, v_GB, v_BG = sy.symbols('p pp q s a x_L x_S x_LL x_LS x_SS  x wG wB k N low kappa v_GB v_BG')


B = ((1 - q) * wG * v_GB - q * wB * v_BG) * (wB - wG)
C = wG * wB * (v_GB + v_BG)
V = wG * v_GB + wB * v_BG

hx = low * (x / sy.Rational(4, 5)) ** kappa

sigma = sy.summation(((-V) ** k ) * sy.binomial(N, k + 1) * x ** k, (k, 0, N - 1))

h = low * (x /sy.Rational(3, 5))**kappa
fitness = - (a * s * (x + hx) + 1 - a) * (C * N / V - (B / V) * sigma) + N
xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(4, 5)

substitutions = [(wB, sy.Rational(9, 10)),
                 (q, sy.Rational(4, 5)), 
                 (wG, sy.Rational(1, 2)), 
                 (a, sy.Rational(4, 10)), 
                 (s, sy.Rational(35, 100)),
                 (v_BG, sy.Rational(4, 5)),
                 (v_GB,  sy.Rational(1, 5)),
                  (low, sy.Rational(1, 2)), 
                  (kappa, 20)]


fig = plt.figure(facecolor='none')
ax = fig.add_subplot(1, 1, 1)
ax.set_facecolor('none') 
NN = []
SSbox = []
LSbox = []
LLbox = []
subs_round = fitness.subs(substitutions)
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

current_winner = None
start_index = None

for i in range(len(NN)):
    values = [LLbox[i], LSbox[i], SSbox[i]]
    max_val = max(values)
    winner = None
    
    if max_val == LLbox[i]:
        winner = 'LL'
    elif max_val == SSbox[i]:
        winner = 'SS'
    elif max_val == LSbox[i]:
        winner = 'LS'
    

    if winner != current_winner:

        if current_winner is not None and start_index is not None:
            if current_winner == 'LL':
                ax.axvspan(NN[start_index], NN[i-1], color=LL_color, alpha=0.4)
            elif current_winner == 'SS':
                ax.axvspan(NN[start_index], NN[i-1], color=SS_color, alpha=0.4)
            elif current_winner == 'LS':
                ax.axvspan(NN[start_index], NN[i-1], color=LS_color, alpha=0.1)
        
        current_winner = winner
        start_index = i


if current_winner is not None and start_index is not None:
    if current_winner == 'LL':
        ax.axvspan(NN[start_index], NN[-1], color=LL_color, alpha=0.4)
    elif current_winner == 'SS':
        ax.axvspan(NN[start_index], NN[-1], color=SS_color, alpha=0.4)
    elif current_winner == 'LS':
        ax.axvspan(NN[start_index], NN[-1], color=LS_color, alpha=0.1)

ax.plot(NN, LLbox, "o", color=LL_color, markeredgecolor='black', markeredgewidth=1, 
        markersize=4, alpha=0.8, label='LL')
ax.plot(NN, LSbox, "s", color=LS_color, markeredgecolor='black', markeredgewidth=1, 
        markersize=4, alpha=0.8, label='LS')
ax.plot(NN, SSbox, "v", color=SS_color, markeredgecolor='black', markeredgewidth=1, 
        markersize=4, alpha=0.8, label='SS')


plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure2f.pdf", format="pdf", 
            transparent=True)
plt.show()