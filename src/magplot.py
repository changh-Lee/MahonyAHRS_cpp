from cProfile import label
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# magnetic field strength(nT) at 22.588642995403102N, 113.96927123871706E
Hx = 37697.2
Hy = -2258.6
Hz = 25756.5
Hr = np.sqrt(Hx**2 + Hy**2)
r = np.sqrt(1- Hz**2/Hr**2)
h = Hz/Hr

# 读取CSV文件
mh = pd.read_csv('data/m_horizon_lse.csv', header=None)
data = pd.read_csv('data/m_calibrated.csv', header=None)
m = np.loadtxt('data/m_projection_lse.csv', delimiter=',')
m = np.concatenate((m, np.ones((m.shape[0], 1))), axis=1)

# 绘制数据
plt.plot(np.sqrt(data.iloc[:, 0]**2 + data.iloc[:, 1]**2), label='calibrated')
plt.plot(np.sqrt(mh.iloc[:, 0]**2 + mh.iloc[:, 1]**2), label='original')
plt.legend()
plt.grid(True)
plt.title('Magnetic Norm')
plt.xlabel('Sample')
plt.ylabel('Horizontal Component(nT)')
# plt.show()

# original data
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(m[:, 0], m[:, 1], m[:, 2], c='r', label='Original Data', s=0.5)
ax.legend()
ax.set_box_aspect([1,1,1])  # Set equal aspect ratio for all axes
ax.set_xlabel('mx')
ax.set_ylabel('my')
ax.set_zlabel('mz')
plt.title('Collected Magnetometer Data')

# calibrated data
data.iloc[:,:2] *= r
data.iloc[:,2] = h
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(data.iloc[:, 0], data.iloc[:, 1], data.iloc[:, 2], c='g', label='Processed Result')

# Theoretical Value
theta = np.linspace(0, 2*np.pi, 100)
x = np.cos(theta)
y = np.sin(theta)
z = np.zeros_like(theta)
ax.plot(x, y, z, 'k--')
ax.plot([0, 0], [-1, 1], [0, 0], 'k--')
ax.plot([-1, 1], [0, 0], [0, 0], 'k--')

x *= r
y *= r
z = h * np.ones_like(theta)
ax.plot(x, y, z, 'b-', label='Theoretical Value')
ax.plot([0, 0], [-r, r], [h, h], 'b--')
ax.plot([-r, r], [0, 0], [h, h], 'b--')


# 绘制单位球
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = np.outer(np.cos(u), np.sin(v))
y = np.outer(np.sin(u), np.sin(v))
z = np.outer(np.ones(np.size(u)), np.cos(v))
ax.plot_surface(x, y, z, color='lightgray', alpha=0.5)
ax.scatter([0], [0], [0], color='r')
ax.scatter([0], [0], [h], color='r')
ax.plot([0, 0], [0, 0], [0, h], 'r-', label='Hz/Hr')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
ax.set_box_aspect([1,1,1])
ax.view_init(elev=90., azim=-90)
plt.title('Collaborated Magnetometer Data')
plt.show()