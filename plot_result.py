import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

# 1. 計算結果の読み込み（1行目のヘッダは無視）
data = np.loadtxt('output/result.dat', comments='#')
x = data[:, 0]
y = data[:, 1]
phi = data[:, 2] # 計算されたポテンシャル（温度）

# 2. 点のデータから三角形の網目を自動生成して描画
triang = mtri.Triangulation(x, y)

plt.figure(figsize=(8, 6))
# 3. 結果を色で塗りつぶす（等高線図）
plt.tricontourf(triang, phi, cmap='jet', levels=20)
plt.colorbar(label='Phi (Potential)')

plt.title('Laplace Equation FEM Result')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal') # 縦横の比率を揃える
plt.savefig('output/result.png') # 画像として保存
print("Image saved to output/result.png")