# 依存性を確認
import matplotlib.pyplot as plt
import numpy as np
import sympy as sy
import os
import resource

p, pp, q, s, a, x_L, x_S, x_LL, x_LS, x_SS, x, wG, wB, k, N, low, kappa = sy.symbols('p pp q s a x_L x_S x_LL x_LS x_SS x wG wB k N low kappa')

# 関数の定義
v_GB = q
v_BG = 1 - q

B = ((1 - q) * wG * v_GB - q * wB * v_BG) * (wB - wG)
C = wG * wB * (v_GB + v_BG)
V = wG * v_GB + wB * v_BG

# シグマの定義
sigma = sy.Sum(((-V) ** k ) * sy.binomial(N, k + 1) * x ** k, (k, 0, N - 1))

# フィットネス関数の定義
fitness = - (a * s * x + 1 - a) * (C * N / V - (B / V) * sigma) + N

# 置換リストの定義
subsses = [
    (wB, sy.Rational(4, 5)),
    (a, sy.Rational(1, 5)),
    (wG, sy.Rational(1, 5)),
]
sulm_round = fitness.subs(subsses)

# 置換後のフィットネス関数の表示
print(sulm_round)

fig = plt.figure(figsize=(12, 10))
ax = []

xLL, xLS, xSS = sy.Rational(1, 5), sy.Rational(2, 5), sy.Rational(4, 5)

# 各プロット用のNの値を定義
Nv = [1, 2, 3, 10, 70, 1000]

for l in range(6):
    ax.append(fig.add_subplot(2, 3, l + 1))     
    sulm_roundN = sulm_round.subs(N, Nv[l]).doit()
    ax[l].set_title("$N=$" + str(Nv[l]))
    ax[l].set_xlabel("$q$")
    ax[l].set_ylabel("$s$")
    ax[l].set_xlim(0, 1)
    ax[l].set_ylim(0, 1)
    LL = sulm_roundN.subs(x, xLL)
    LS = sulm_roundN.subs(x, xLS)
    SS = sulm_roundN.subs(x, xSS)
    qq = np.array([])
    SSLSbox = np.array([])
    LSLLbox = np.array([])
    SSLLbox = np.array([])
    for q_jj in range(0, 101):
        q_j = sy.Rational(q_jj, 100)
        qq = np.append(qq, float(q_j))  # floatに変換
        SS_sub = SS.subs(q, q_j)
        LS_sub = LS.subs(q, q_j)
        LL_sub = LL.subs(q, q_j)
        SSLS = sy.solve(SS_sub - LS_sub, s)[0]  # SS - LS = 0 の解を求める
        LSLL = sy.solve(LS_sub - LL_sub, s)[0]  # LS - LL = 0 の解を求める
        SSLL = sy.solve(SS_sub - LL_sub, s)[0]  # SS - LL = 0 の解を求める
        SSLSbox = np.append(SSLSbox, float(sy.re(SSLS)))
        LSLLbox = np.append(LSLLbox, float(sy.re(LSLL)))
        SSLLbox = np.append(SSLLbox, float(sy.re(SSLL)))

    if l == 0:  # N=1の場合
        ax[l].fill_between(qq, 0, 1, color="b", alpha=0.4)
    else:
        ax[l].fill_between(qq, LSLLbox, 1, color="b", where=LSLLbox > SSLSbox, alpha=0.4)
        # LLの領域
        ax[l].fill_between(qq, SSLSbox, 1, color="b", where=LSLLbox < SSLSbox, alpha=0.4)
        # 双安定領域
        ax[l].fill_between(qq, SSLSbox, LSLLbox, color="w", where=SSLSbox > LSLLbox, hatch="//",alpha=0.6, edgecolor="r")
        # 双安定領域
        ax[l].fill_between(qq, SSLSbox, LSLLbox, color="w", where=SSLSbox > LSLLbox, hatch="\\",alpha=0.6, edgecolor="b")
        # 共存領域
        ax[l].fill_between(qq, SSLSbox, LSLLbox, color="y", edgecolor="0.2", where=SSLSbox < LSLLbox, hatch="..", alpha=0.6)
        # SS領域 
        ax[l].fill_between(qq, 0, SSLSbox, color="r", where=LSLLbox > SSLSbox, alpha=0.4)
        # SS領域
        ax[l].fill_between(qq, 0, LSLLbox, color="r", where=LSLLbox < SSLSbox, alpha=0.4)

plt.tight_layout()
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure4.pdf", format="pdf")
plt.show()
