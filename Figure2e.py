import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os

# Declare symbolic variables
sy.var("p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, v_GB, v_BG, x, wG, wB, N")

# シンボルの定義
p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, x, wG, wB, k, N, low, kappa, v_GB, v_BG = sy.symbols('p pp q s a x_L x_S x_LL x_LS x_SS  x wG wB k N low kappa v_GB v_BG')

# 関数の定義
#v_GB = q
#v_BG = 1 - q

B = ((1 - q) * wG * v_GB - q * wB * v_BG) * (wB - wG)
C = wG * wB * (v_GB + v_BG)
V = wG * v_GB + wB * v_BG

# シグマの定義
sigma = sy.Sum(((-V) ** k ) * sy.binomial(N, k + 1) * x ** k, (k, 0, N - 1))

fitness = - (a * s * x  + 1 - a) * (C * N / V - (B / V) * sigma) + N




# Substitute constants
xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(4, 5)
substitutions = [(wB, sy.Rational(9, 10)),
                 (q, sy.Rational(4, 5)), 
                 (wG, sy.Rational(1, 2)), 
                 (a, sy.Rational(4, 10)), 
                 (s, sy.Rational(35, 100)),
                 (v_BG, sy.Rational(4, 5)),
                 (v_GB,  sy.Rational(1, 5))]



# 元のコードの前半部分は同じなので省略...

# Prepare for plotting
fig = plt.figure(facecolor='none')
ax = fig.add_subplot(1, 1, 1)
ax.set_facecolor('none') 
#ax.set_xlabel("$N$")
#ax.set_ylabel("$Fitness$")
NN = []
SSbox = []
LSbox = []
LLbox = []
subs_round = fitness.subs(substitutions)

# 色の定義
LL_color = 'c'
LS_color = 'yellow'
SS_color = 'red'

# データポイントの計算と格納
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

# 戦略の優位性を判定し、連続する区間を特定
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
    
    # 新しい優位戦略の開始
    if winner != current_winner:
        # 前の区間を終了
        if current_winner is not None and start_index is not None:
            if current_winner == 'LL':
                ax.axvspan(NN[start_index], NN[i-1], color=LL_color, alpha=0.4)
            elif current_winner == 'SS':
                ax.axvspan(NN[start_index], NN[i-1], color=SS_color, alpha=0.4)
            elif current_winner == 'LS':
                ax.axvspan(NN[start_index], NN[i-1], color=LS_color, alpha=0.1)
        
        # 新しい区間の開始
        current_winner = winner
        start_index = i

# 最後の区間の処理
if current_winner is not None and start_index is not None:
    if current_winner == 'LL':
        ax.axvspan(NN[start_index], NN[-1], color=LL_color, alpha=0.4)
    elif current_winner == 'SS':
        ax.axvspan(NN[start_index], NN[-1], color=SS_color, alpha=0.4)
    elif current_winner == 'LS':
        ax.axvspan(NN[start_index], NN[-1], color=LS_color, alpha=0.1)

# データポイントのプロット
ax.plot(NN, LLbox, "o", color=LL_color, markeredgecolor='black', markeredgewidth=1, 
        markersize=4, alpha=0.8, label='LL')
ax.plot(NN, LSbox, "s", color=LS_color, markeredgecolor='black', markeredgewidth=1, 
        markersize=4, alpha=0.8, label='LS')
ax.plot(NN, SSbox, "v", color=SS_color, markeredgecolor='black', markeredgewidth=1, 
        markersize=4, alpha=0.8, label='SS')

#plt.yscale('log')
#plt.legend()
plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure2e.pdf", format="pdf", 
            transparent=True)
plt.show()