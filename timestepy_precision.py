import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import cm
from ipywidgets import interact, FloatSlider
import numpy as np
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
#スタイルの設定
sns.set(font="TakaoPGothic", palette="colorblind", style="whitegrid")
#datファイルからデータを読み取る
axisx, axisy = np.loadtxt("./chk_timestepx_precision.dat", comments='!', unpack=True)
#x, yの値の設定
x = axisx
y = axisy
log_x = [np.log(i) for i in x]
log_y = [np.log(i) for i in y]
#回帰直線
a, b = np.polyfit(log_x, log_y, 1)
print('a is :', a)
print('b is :', b)
# 1. Figureのインスタンスを生成
fig = plt.figure()
# 2. Axesのインスタンスを生成
ax = fig.add_subplot(111)
# 3. データを渡してプロット
line1, = ax.plot(x, y, label = 'timestepy_precision', marker='.')
# 4. グラフタイトル, ラベル付け等
ax.set_xlim(0.0125, 0.1)
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel("dt")
ax.text(0.0405, 0.00000043, 'slope : 0.979')
ax.set_title("check the timestepy precision")
ax.legend()
# 5. グラフを描画
plt.grid(which="both")
plt.show()