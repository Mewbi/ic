import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

class Element:

    def __init__(self):
        self.X = 0
        self.Y = 0
        self.Z = 0

def calculate_coords(R1, R2, R3, A1, A2, D1):
    # Converter ângulos para radianos
    A1_rad = np.radians(A1)
    A2_rad = np.radians(A2)
    D1_rad = np.radians(D1)

    # Elements
    H1 = Element()
    H2 = Element()
    O  = Element()
    F  = Element()

    # Calcular coordenadas dos pontos
    O.X = 0
    O.Y = 0
    O.Z = 0

    H1.X = 0
    H1.Y = R1
    H1.Z = 0
    
    H2.X = np.sin(A1_rad) * R2
    H2.Y = np.cos(A1_rad) * R2
    H2.Z = 0

    F.Y = (R1 - R3 * np.cos(A2_rad)) / R1
    # R3 ** 2 = F.X ** 2 + (F.Y - H1.Y) ** 2 + F.Z ** 2
    # np.cos(D1_rad) ** 2 = ( -F.X ** 2 ) / (F.X ** 2 + F.Z ** 2)
    m = np.array([
        [- np.sin(D1_rad) ** 2, np.cos(D1_rad) ** 2],
        [1, 1]
    ])

    m2 = np.array(
        [0, R3 ** 2 - (F.Y -R1) ** 2]
    )

    x = linalg.inv(m) @ m2

    F.X = np.sqrt(x[0])
    F.Z = np.sqrt(x[1])

    return H1, H2, O, F

def plot_config(R1, R2, R3, A1, A2, D1):

    H1, H2, O, F = calculate_coords(R1, R2, R3, A1, A2, D1)

    x = np.array([H1.X, H2.X, O.X, F.X])
    y = np.array([H1.Y, H2.Y, O.Y, F.Y])
    z = np.array([H1.Z, H2.Z, O.Z, F.Z])

    # Configurar o plot 3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel('Eixo X')
    ax.set_ylabel('Eixo Y')
    ax.set_zlabel('Eixo Z')

    # Plotar os pontos
    ax.scatter(x, y, z, c='r', marker='o')

    # Plotar as linhas entre os pontos
    ax.plot([H1.X, O.X], [H1.Y, O.Y], [H1.Z, O.Z], c='b')
    ax.plot([H2.X, O.X], [H2.Y, O.Y], [H2.Z, O.Z], c='b')
    ax.plot([H1.X, F.X], [H1.Y, F.Y], [H1.Z, F.Z], c='b')

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

'''
R1 = 0.9609
R2 = 0.9609
R3 = 10
A1 = 104.1477
A2 = 300
D1 = 300
'''
# Chamar a função para plotar a configuração
plot_config(R1, R2, R3, A1, A2, D1)
