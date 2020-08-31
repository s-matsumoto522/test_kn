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
axisx1, axisx2, axisx3, axisy1, axisy2, axisy3 = np.loadtxt("./chk_viscosity_precision.dat", comments='!', unpack=True)
#x, yの値の設定
x1 = axisx1
y1 = axisy1
x2 = axisx2
y2 = axisy2
x3 = axisx3
y3 = axisy3
log_x1 = [np.log(i) for i in x1]
log_y1 = [np.log(i) for i in y1]
log_x2 = [np.log(i) for i in x2]
log_y2 = [np.log(i) for i in y2]
log_x3 = [np.log(i) for i in x3]
log_y3 = [np.log(i) for i in y3]
#回帰直線
linear1 = np.polyfit(log_x1, log_y1, 1)
linear2 = np.polyfit(log_x2, log_y2, 1)
linear3 = np.polyfit(log_x3, log_y3, 1)
print(linear1[0], linear2[0], linear3[0])
# 1. Figureのインスタンスを生成
fig1 = plt.figure()
fig2 = plt.figure()
fig3 = plt.figure()
# 2. Axesのインスタンスを生成
ax1 = fig1.add_subplot(111)
ax2 = fig2.add_subplot(111)
ax3 = fig3.add_subplot(111)
# 3. データを渡してプロット
line1, = ax1.plot(x1, y1, label = 'x_viscosity_precision', marker='.')
line2, = ax2.plot(x2, y2, label = 'y_viscosity_precision', marker='.')
line3, = ax3.plot(x3, y3, label = 'z_viscosity_precision', marker='.')
# 4. グラフタイトル, ラベル付け等
ax1.set_xlim(2.0*np.pi/160, 2.0*np.pi/20)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlabel("dX")
ax1.text(0.1, 0.0012, 'slope : 1.968')
ax1.grid(which='both')
ax1.set_title("check the x_viscosity precision")
ax1.legend()
ax2.set_xlim(2.0*np.pi/160, 2.0*np.pi/20)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlabel("dY")
ax2.text(0.1, 0.0012, 'slope : 1.968')
ax2.grid(which='both')
ax2.set_title("check the y_viscosity precision")
ax2.legend()
ax3.set_xlim(2.0*np.pi/160, 2.0*np.pi/20)
ax3.set_xscale("log")
ax3.set_yscale("log")
ax3.set_xlabel("dZ")
ax3.text(0.1, 0.012, 'slope : 1.967')
ax3.grid(which='both')
ax3.set_title("check the z_viscosity precision")
ax3.legend()
# 5. グラフを描画
plt.show()
