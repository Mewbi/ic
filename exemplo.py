import numpy as np
import li_dawes_guo as ldg


# Para inicializar:
#  Isso PRECISA ser chamado antes de qualquer chamada à função.
#  Basicamente ele lê os valores dos arquivos .txt que definem essa função
ldg.init()  

# Geometria da tabela I, artigo Li, Dawes e Guo, JCP 137, 094304 (2012)
f_h2o = np.array([0.9613, 0.9613, 50, 104.25, 0, 0])
R_vdw = np.array([0.9668, 0.9668, 2.2992, 105.07, 70.54, -83.05])
R_ts = np.array([0.9705, 1.0389, 1.3128, 103.00, 121.99, 68.72])
P_vdw = np.array([0.9738, 1.8037, 0.9330, 111.48, 179.16, 0])
ho_hf = np.array([0.9730, 10, 0.9204, 0, 0, 0 ])
P_vdw_ohfh = np.array([0.9778, 3.5351, 0.9199, 17.51, 46.76, 0])

points = [
    np.array([0.9613, 0.9613, 50, 104.25, 0, 0]),
    np.array([0.9668, 0.9668, 2.2992, 105.07, 70.54, -83.05]),
    np.array([0.9705, 1.0389, 1.3128, 103.00, 121.99, 68.72]),
    np.array([0.9738, 1.8037, 0.9330, 111.48, 179.16, 0]),
    np.array([0.9730, 10, 0.9204, 0, 0, 0 ]),
    np.array([0.9778, 3.5351, 0.9199, 17.51, 46.76, 0])
]

for i, p in enumerate(points):
    E_fh2o = ldg.pes(p)
    print(f'Energia do f_h2o: {E_fh2o}')


exit()
