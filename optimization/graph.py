import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def plot_graph(fn, p_x, p_y, p_z):
    # Plot surface
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    min_pt = min(p_x+p_y)
    max_pt = max(p_x+p_y)

    x = y = np.arange(-abs(min_pt*2), max_pt*2, 1)

    X, Y = np.meshgrid(x, y)
    zs = np.array([fn([x,y]) for x,y in zip(np.ravel(X), np.ravel(Y))])
    Z = zs.reshape(X.shape)
    ax.plot_wireframe(X, Y, Z, rstride=1, cstride=1)

    # Plot points
    ax.scatter(p_x, p_y, p_z, c=p_z, linewidth=2, marker='o', cmap='viridis')

    plt.show()
