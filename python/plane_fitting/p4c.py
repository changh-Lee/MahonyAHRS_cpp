import numpy as np
import numpy.linalg as la
from svd_solve import svd, svd_solve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fit_plane_LSE import fit_plane_LSE, fit_plane_LSE_RANSAC

# points_file = open('cluttered_table.txt', 'r')
# points_lines = points_file.readlines()
# points_array = [[float(s) for s in line.strip().split(' ')] for line in points_lines]
# m = np.array(points_array, dtype=np.float32)
# # homogeneous coord
# m = np.concatenate((m, np.ones((m.shape[0], 1))), axis=1)

m = np.loadtxt('../../Magnetometer_Calibration_HalII_matlab/ellipsoid_fitting/slope_points.csv', delimiter=',')# mag_out_sample.txt plane_points.csv slope_points.csv
m = np.concatenate((m, np.ones((m.shape[0], 1))), axis=1)
m_projection = np.zeros(m[:,:3].shape)

p, inlier_list = fit_plane_LSE_RANSAC(m)


# x = np.arange(np.min(m[:, 0]), np.max(m[:, 0]), 0.1)
# z = np.arange(np.min(m[:, 2]), np.max(m[:, 2]), 0.1)

# xx, zz = np.meshgrid(x, z)

# yy = (-p[0]*xx -p[2]*zz -p[3])/p[1]

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(m[:, 0], m[:, 1], m[:, 2], c='r', label='Original Data')
# ax.plot_wireframe(xx, yy, zz)

def draw_plane(p, m, ax):
    n = p[:3] / np.linalg.norm(p[:3])
    print('normal vector:', n)
    normal_heading_axis = np.argmax(np.abs(n))
    if normal_heading_axis == 0:
        # grid with yz
        y_min = np.min(m[:, 1])
        y_max = np.max(m[:, 1])
        z_min = np.min(m[:, 2])
        z_max = np.max(m[:, 2])
        y_avg = (y_min + y_max)/2
        z_avg = (z_min + z_max)/2
        
        y = np.arange(y_min, y_max, 0.1)
        z = np.arange(z_min, z_max, 0.1)
        yy, zz = np.meshgrid(y, z)
        xx = (-p[1]*yy - p[2]*zz - p[3])/p[0]
        x_avg = (-p[1]*y_avg - p[2]*z_avg - p[3])/p[0]

        m_projection = magne_projection(m, n, x_avg, y_avg, z_avg)
        ax.scatter(m_projection[:, 0], m_projection[:, 1], m_projection[:, 2], c='g', label='Processed Result')
        m_horizon = magne_horizon(n, m_projection)
        ax.scatter(m_horizon[:, 0], m_horizon[:, 1], m_horizon[:, 2], c='b', label='Horizon Result')
        ax.plot_surface(xx, yy, zz)
    elif normal_heading_axis == 1:
        # grid with xz
        x_min = np.min(m[:, 0])
        x_max = np.max(m[:, 0])
        z_min = np.min(m[:, 2])
        z_max = np.max(m[:, 2])
        x_avg = (x_min + x_max)/2
        z_avg = (z_min + z_max)/2
        
        x = np.arange(x_min, x_max, 0.1)
        z = np.arange(z_min, z_max, 0.1)
        xx, zz = np.meshgrid(x, z)
        yy = (-p[0]*xx - p[2]*zz - p[3])/p[1]
        y_avg = (-p[0]*x_avg - p[2]*z_avg - p[3])/p[1]

        m_projection = magne_projection(m, n, x_avg, y_avg, z_avg)
        ax.scatter(m_projection[:, 0], m_projection[:, 1], m_projection[:, 2], c='g', label='Processed Result')
        m_horizon = magne_horizon(n, m_projection)
        ax.scatter(m_horizon[:, 0], m_horizon[:, 1], m_horizon[:, 2], c='b', label='Horizon Result')
        ax.plot_surface(xx, yy, zz)
    elif normal_heading_axis == 2:
        # grid with xy
        x_min = np.min(m[:, 0])
        x_max = np.max(m[:, 0])
        y_min = np.min(m[:, 1])
        y_max = np.max(m[:, 1])
        x_avg = (x_min + x_max)/2
        y_avg = (y_min + y_max)/2
        
        x = np.arange(x_min, x_max, 0.1)
        y = np.arange(y_min, y_max, 0.1)
        xx, yy = np.meshgrid(x, y)
        zz = (-p[0]*xx - p[1]*yy - p[3])/p[2]
        z_avg = (-p[0]*x_avg - p[1]*y_avg - p[3])/p[2]

        m_projection = magne_projection(m, n, x_avg, y_avg, z_avg)
        ax.scatter(m_projection[:, 0], m_projection[:, 1], m_projection[:, 2], c='g', label='Processed Result')
        m_horizon = magne_horizon(n, m_projection)
        ax.scatter(m_horizon[:, 0], m_horizon[:, 1], m_horizon[:, 2], c='b', label='Horizon Result')
        ax.plot_surface(xx, yy, zz)

# magnetic m projection
def magne_projection(m, n, x, y, z):
    # m: magnetic m, n: n vector
    # return: m_projection
    point = [x, y, z]
    for i in range(m.shape[0]):
        m_projection[i,:] = m[i,:3] - np.dot(n.transpose(), (m[i,:3] - point)) * n
    return m_projection


# projected magnetic points horizontal projection
def magne_horizon(n, m):
    # n: normal vector, m: magnetic points
    # return: m_projection
    z = [0, 0, -1]# positive direction of z is up here, but the z vector in the input must be down
    angle_rad = np.arccos(np.dot(n, z) / (np.linalg.norm(n) * np.linalg.norm(z)))# [0,pi]
    axis = np.cross(n, z)
    if angle_rad > np.pi/2:
        axis = -axis
        angle_rad = np.pi - angle_rad
    angle_deg = np.degrees(angle_rad)
    print('rotation angle:', angle_deg)
    axis = axis / np.linalg.norm(axis)

    # Rodrigues' rotation 
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    t = 1 - c
    x, y, z = axis

    rotation_matrix = np.array([[t*x*x + c, t*x*y - s*z, t*x*z + s*y],
                                [t*x*y + s*z, t*y*y + c, t*y*z - s*x],
                                [t*x*z - s*y, t*y*z + s*x, t*z*z + c]])
    
    m_horizon = (rotation_matrix @ m.transpose()).transpose()
    file_name = 'm_horizon_RANSAC.csv'
    np.savetxt(file_name, m_horizon, fmt='%.6f', delimiter=',')
    return m_horizon

draw_plane(p, m, ax)
ax.legend()
plt.show()