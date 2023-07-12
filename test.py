import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot_config(R1, R2, R3, A1, A2, D1):
    # Converter ângulos para radianos
    A1_rad = np.radians(A1)
    A2_rad = np.radians(A2)
    D1_rad = np.radians(D1)

    # Calcular coordenadas dos pontos
    x = np.array([0, R1, R3*np.sin(A2_rad), R3*np.sin(A2_rad) + R2*np.sin(D1_rad)])
    y = np.array([0, 0, R3*np.cos(A2_rad), R2*np.cos(A1_rad) + R1*np.cos(A1_rad)])
    z = np.array([0, R2*np.cos(A1_rad), R2*np.cos(A1_rad) + R1*np.sin(D1_rad), R3*np.cos(A2_rad)])

    # Configurar o plot 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('Eixo X')
    ax.set_ylabel('Eixo Y')
    ax.set_zlabel('Eixo Z')

    # Plotar os pontos
    ax.scatter(x, y, z, c='r', marker='o')

    # Plotar as linhas entre os pontos
    ax.plot([x[0], x[1]], [y[0], y[1]], [z[0], z[1]], c='b')
    ax.plot([x[0], x[2]], [y[0], y[2]], [z[0], z[2]], c='b')
    ax.plot([x[1], x[2]], [y[1], y[2]], [z[1], z[2]], c='b')
    ax.plot([x[0], x[3]], [y[0], y[3]], [z[0], z[3]], c='b')
    ax.plot([x[1], x[3]], [y[1], y[3]], [z[1], z[3]], c='b')
    ax.plot([x[2], x[3]], [y[2], y[3]], [z[2], z[3]], c='b')

    # Ajustar a escala dos eixos para melhor visualização
    ax.auto_scale_xyz([np.min(x), np.max(x)], [np.min(y), np.max(y)], [np.min(z), np.max(z)])

    # Mostrar o plot
    plt.show()

# Valores de exemplo
R1 = 0.9735
R2 = 3.9974
R3 = 0.9220
A1 = 17.51
A2 = 46.76
D1 = 0

R1 = 0.9609
R2 = 0.9609
R3 = 10
A1 = 104.1477
A2 = 300
D1 = 300

# Chamar a função para plotar a configuração
plot_config(R1, R2, R3, A1, A2, D1)
