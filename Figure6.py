import matplotlib.pyplot as plt
import numpy as np
import os

def B(q, v_BG, v_GB):
  return ((1 - q) * wG * v_GB - q * wB * v_BG) * (wB - wG)

def C(v_BG, v_GB):
  return wG * wB * (v_GB + v_BG)

def V(v_BG, v_GB):
  return wG * v_GB + wB * v_BG

def combination(a, b):
    result = 1
    for i in range(b):
        result *= (a-i)/(b-i)
    return result

def sigma(v_BG, v_GB, x):
  ans = 0
  for k in range(N):
    ans += (-V(v_BG, v_GB)) ** k * combination(N, k + 1) * x ** k
  return ans

def hx(low, kappa, x):
  return low * (x / 4/5) ** kappa

def fitness(q, v_BG, v_GB, x):
  return - ( a * s * (x + hx(1 / 2, 5, x)) + 1 - a) * (C(v_BG, v_GB) * N / V(v_BG, v_GB) - (B(q, v_BG, v_GB) / V(v_BG, v_GB)) * sigma(v_BG, v_GB, x)) + N


N = 70
wG = 0.2
wB = 0.8
a = 0.2
xLL = 0.2
xLS = 2/5
xSS = 4/5
s = 0.3

print(fitness)

Sallelefreq = np.array([])
q_list = np.array([])

for q_l in np.arange(1, -0.01, -0.01):
  q_list = np.append(q_list, 1-q_l)
  if fitness(q_l, 1-q_l, q_l, xLS) > fitness(q_l, 1-q_l, q_l, xLL) and fitness(q_l, 1-q_l, q_l, xLS) > fitness(q_l, 1-q_l, q_l, xSS):
    Sallelefreq = np.append(Sallelefreq, (fitness(q_l, 1-q_l, q_l, xLL) - fitness(q_l, 1-q_l, q_l, xLS))/ (fitness(q_l, 1-q_l, q_l, xLL) + fitness(q_l, 1-q_l, q_l, xSS) - 2 * fitness(q_l, 1-q_l, q_l, xLS)))
  elif fitness(q_l, 1-q_l, q_l, xSS) > fitness(q_l, 1-q_l, q_l, xLS) and fitness(q_l, 1-q_l, q_l, xSS) > fitness(q_l, 1-q_l, q_l, xLL):
    Sallelefreq = np.append(Sallelefreq, 1)
  elif fitness(q_l, 1-q_l, q_l, xLL) > fitness(q_l, 1-q_l, q_l, xLS) and fitness(q_l, 1-q_l, q_l, xLL) > fitness(q_l, 1-q_l, q_l, xSS):
    Sallelefreq = np.append(Sallelefreq, 0)

plt.figure(figsize=(6, 10))


max_index = np.argmax(Sallelefreq)
plt.plot(Sallelefreq[:max_index+1], q_list[:max_index+1], color='green', label='S Allele Frequency (Green)')
plt.plot(Sallelefreq[max_index+1:], q_list[max_index+1:], color='blue', label='S Allele Frequency (Blue)')
plt.xlabel('Allele Frequency')
plt.ylabel('q')
plt.title('Allele Frequencies')
plt.legend()
plt.gca().patch.set_alpha(0)
plt.savefig(os.path.dirname(os.path.abspath(__file__)) + "/Figure6_分割_line.pdf", format="pdf", facecolor='none', edgecolor='none')
plt.show()