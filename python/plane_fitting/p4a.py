import math
import numpy as np
import numpy.linalg as la
from svd_solve import svd, svd_solve
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from fit_plane_LSE import fit_plane_LSE
import scipy
import time
import sys
import csv
import msvcrt

def read_Xsens(file_path, last_line_number):
    all_data = []
    new_data = []
    with open(file_path, 'r') as file:
        lines = file.readlines()
        # all_data = [line.strip().split(',') for line in lines]
        all_data = [[float(line.strip().split('\t')[i]) for i in range(1, 4)] for line in lines[14:]]

    new_data_lines = all_data[last_line_number:]
    for line in new_data_lines:
        # data = [float(d) for d in line]
        data = [float(line[i]) for i in range(0, 3)]# get magdata columns
        new_data.append(data)
    
    # Update last_line_number
    last_line_number += len(new_data_lines)
        
    return new_data, all_data, last_line_number


def plot_data(data, ax):
    data = np.array(data)
    ax.scatter(data[:, 0], data[:, 1], data[:, 2], s=1, label='Xsens MTi-G-710 MagData')
    # ax.set_aspect('equal')
    plt.draw() 
    plt.legend()
    plt.pause(0.01)

def execute_command(user_input, m):
    if user_input.lower() == 'y':
        m = np.concatenate((m, np.ones((m.shape[0], 1))), axis=1)

        print("Complete collection, plan fitting...")
        plot_fitted_plane(m)

        sys.exit("Plan Fitting Finished")
    elif user_input.lower() == '\x1b':
        sys.exit("Python program terminated by user")
    elif user_input.lower() == 'd':
        print("Sample Count: ", m.shape[0])
    else:
        print("Unknown command:", user_input)

def PlanFitting(file_path):
    if not Simulation:
        print("Start Collecting Magnetometer Dataset: ")
        print("Press 'Y' to start plan fitting, 'esc' to exit:")
        last_line_number = 0
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        plt.title('Collected Magnetometer Data')

        last_print_time = time.time()  # 初始化last_print_time

        while True:
            new_data, all_data, last_line_number = read_Xsens(file_path, last_line_number)
            # Plot new data if available
            if new_data:
                plot_data(new_data, ax)# 可实时绘图
            else:
                waiting_time = time.time()
                if waiting_time - last_print_time >= 10:
                    print("Waiting for new data...")
                    last_print_time = waiting_time

            # windows msvcrt/ linux select?
            if msvcrt.kbhit():
                key = msvcrt.getch().decode('utf-8').lower()
                execute_command(key, np.array(all_data, dtype=float))

            time.sleep(1)  # Adjust the delay as needed

    else:
        file_path = "data/m_slope.csv"  # collected dataset
        m = np.loadtxt(file_path, delimiter=',')# mag_out_sample.txt plane_points.csv slope_points.csv
        m = np.concatenate((m, np.ones((m.shape[0], 1))), axis=1)
        plot_fitted_plane(m)


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
        ax.scatter(m_projection[:, 0], m_projection[:, 1], m_projection[:, 2], c='g', label='Processed Result', s=0.5)
        m_horizon = magne_horizon(n, m_projection)
        ax.scatter(m_horizon[:, 0], m_horizon[:, 1], m_horizon[:, 2], c='b', label='Horizon Result', s=0.5)
        ax.plot_surface(xx, yy, zz, color='lightgray', alpha=0.5)
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
        ax.scatter(m_projection[:, 0], m_projection[:, 1], m_projection[:, 2], c='g', label='Processed Result', s=0.5)
        m_horizon = magne_horizon(n, m_projection)
        ax.scatter(m_horizon[:, 0], m_horizon[:, 1], m_horizon[:, 2], c='b', label='Horizon Result', s=0.5)
        ax.plot_surface(xx, yy, zz, color='lightgray', alpha=0.5)
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
        file_name = 'data/m_projection_lse.csv'
        np.savetxt(file_name, m_projection, fmt='%.6f', delimiter=',')
        ax.scatter(m_projection[:, 0], m_projection[:, 1], m_projection[:, 2], c='g', label='Processed Result', s=0.5)
        m_horizon = magne_horizon(n, m_projection)
        ax.scatter(m_horizon[:, 0], m_horizon[:, 1], m_horizon[:, 2], c='b', label='Horizon Result', s=0.5)
        ax.plot_surface(xx, yy, zz, color='lightgray', alpha=0.5)

# magnetic points projection
def magne_projection(m, n, x, y, z):
    # m: magnetic points, n: n vector
    # return: m_projection
    point = [x, y, z]
    m_projection = np.zeros(m[:,:3].shape)

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
    print('rotation angle:', angle_deg, angle_rad)
    axis = axis / np.linalg.norm(axis)
    print('rotation axis:', axis)
    with open('data/PlaneFitting_param.csv', 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([angle_rad])
        writer.writerow(axis)

    # Rodrigues' rotation 
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    t = 1 - c
    x, y, z = axis

    rotation_matrix = np.array([[t*x*x + c, t*x*y - s*z, t*x*z + s*y],
                                [t*x*y + s*z, t*y*y + c, t*y*z - s*x],
                                [t*x*z - s*y, t*y*z + s*x, t*z*z + c]])
    print("Plane Fitting Rotation Matrix:\n", rotation_matrix)
    
    m_horizon = (rotation_matrix @ m.transpose()).transpose()
    m_horizon[:,2] = ratio
    file_name = 'data/m_horizon_lse.csv'
    np.savetxt(file_name, m_horizon, fmt='%.6f', delimiter=',')
    return m_horizon

def plot_fitted_plane(m):
        p = fit_plane_LSE(m)

        # calculate error:
        dists = np.abs(m @ p) / np.sqrt(p[0]**2 + p[1]**2 + p[2]**2)
        avg_dist = np.mean(dists)
        print('average distance:',avg_dist)

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        ax.scatter(m[:, 0], m[:, 1], m[:, 2], c='r', label='Original Data', s=0.5)

        draw_plane(p, m, ax)
        plt.title('Rotated Magnetometer Data')
        ax.legend()
        plt.show()


if __name__ == "__main__":
    Simulation = False
    Hx = 37697.2
    Hy = -2258.6
    Hz = 25756.5
    r = np.sqrt(Hx**2 + Hy**2)
    ratio = Hz/r# z-axis offset 

    file_path = "data/MT_07782926_007-000.txt" # Xsens dataset(real-time/collected)
    PlanFitting(file_path)


