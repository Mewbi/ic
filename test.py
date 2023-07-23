import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

class Element:

    def __init__(self, color="red", size=1, label=""):
        self.X = 0
        self.Y = 0
        self.Z = 0
        self.color = color
        self.size = 20 * 4 * size
        self.label = label

def calculate_coords(R1, R2, R3, A1, A2, D1):
    # Convert to radians
    A1_rad = np.radians(A1)
    A2_rad = np.radians(A2)
    D1_rad = np.radians(D1)

    # Elements
    # Colors: https://en.wikipedia.org/wiki/CPK_coloring
    H1 = Element(color="whitesmoke",  size=1, label="H1")
    H2 = Element(color="whitesmoke",  size=1, label="H2")
    O  = Element(color="red",        size=1, label="O")
    F  = Element(color="lightgreen", size=1, label="F")

    # Calculate points coordenates
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

H1, H2, O, F = calculate_coords(R1, R2, R3, A1, A2, D1)

# Dados dos pontos
x = np.array([H1.X, H2.X, O.X, F.X])
y = np.array([H1.Y, H2.Y, O.Y, F.Y])
z = np.array([H1.Z, H2.Z, O.Z, F.Z])

# Configurações do plot 3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_facecolor('darkgrey')
ax.axis('off')

# Função de inicialização do plot (pontos iniciais)
def init():
    ax.scatter(H1.X, H1.Y, H1.Z, c=H1.color, marker='o', s=H1.size)
    ax.scatter(H2.X, H2.Y, H2.Z, c=H2.color, marker='o', s=H2.size)
    ax.scatter(O.X,  O.Y,  O.Z,  c=O.color,  marker='o', s=O.size)
    ax.scatter(F.X,  F.Y,  F.Z,  c=F.color,  marker='o', s=F.size)

    ax.text(H1.X, H1.Y, H1.Z, H1.label, fontweight='bold')
    ax.text(H2.X, H2.Y, H2.Z, H2.label, fontweight='bold')
    ax.text(O.X,  O.Y,  O.Z,  O.label,  fontweight='bold')
    ax.text(F.X,  F.Y,  F.Z,  F.label,  fontweight='bold')

    ax.plot([H1.X, O.X], [H1.Y, O.Y], [H1.Z, O.Z], c='b')
    ax.plot([H2.X, O.X], [H2.Y, O.Y], [H2.Z, O.Z], c='b')
    ax.plot([H1.X, F.X], [H1.Y, F.Y], [H1.Z, F.Z], c='b')

    ax.auto_scale_xyz([np.min(x), np.max(x)], [np.min(y), np.max(y)], [np.min(z), np.max(z)])

# Função para atualizar o plot a cada frame (pode ser adaptada com variações nos pontos)
def update(frame):
    # Exemplo de variação nos pontos (substitua com a lógica desejada)
    H1.X += H1.X * 0.1
    H1.Y += H1.Y * 0.1
    H1.Z += H1.Z * 0.1

    H2.X += H2.X * 0.1
    H2.Y += H2.Y * 0.1
    H2.Z += H2.Z * 0.1

    O.X -= O.X * 0.1
    O.Y -= O.Y * 0.1
    O.Z -= O.Z * 0.1

    F.X -= F.X * 0.1
    F.Y -= F.Y * 0.1
    F.Z -= F.Z * 0.1

    # Limpar o plot para desenhar os novos pontos
    ax.clear()
    init()

# Criar a animação
animation = FuncAnimation(fig, update, frames=range(10), init_func=init, blit=False)

# Salvar a animação em um arquivo GIF
animation.save('plot_animado.gif', writer='pillow', fps=2)

