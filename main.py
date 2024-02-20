"""
Implementation of Reiter's model (refer to https://patarnott.com/pdf/SnowCrystalGrowth.pdf) for creating 2d
ice crystals using cellular automata.
"""

from hexalattice.hexalattice import *


def grid_to_color(grid):
    result = []
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            red = grid[i][j]
            green = grid[i][j]
            blue = grid[i][j]
            result.append((red, green, blue))
    return np.float32(result)


alpha = 1.0
beta = 0.4
gamma = 0.001
rows = 51 # must be uneven
columns = rows

s = np.ones((rows, columns), dtype=float)*beta
f = np.zeros((rows, columns), dtype=float)
s[rows // 2][columns // 2] = 1
'''
s[rows//2][columns//2-1] = 0.7
s[rows//2][columns//2+1] = 0.7
s[rows//2-1][columns//2] = 0.7
s[rows//2+1][columns//2] = 0.7
s[rows//2+1][columns//2+1] = 0.7
s[rows//2-1][columns//2+1] = 0.7
'''

for k in range(25):
    r = np.zeros((rows, columns), dtype=float)
    nr = np.zeros((rows, columns), dtype=float)
    for i in range(rows):
        for j in range(columns):
            if i==0 or j==0 or i==rows-1 or j==columns-1:
                continue
            elif s[i][j] >= 1 or s[i][j - 1] >= 1 or s[i][j + 1] >= 1 or s[i - 1][j] >= 1 or s[i + 1][j] >= 1 or s[i + 1][j + 1] >= 1 or s[i - 1][j + 1] >= 1:
                r[i][j] = s[i][j]
            else:
                nr[i][j] = s[i][j]

    for i in range(rows):
        for j in range(columns):
            if i==0 or j==0 or i==rows-1 or j==columns-1:
                continue
            nr[i][j] = nr[i][j]/2 + nr[i][j - 1]/12 + nr[i][j + 1]/12 + nr[i - 1][j]/12 + s[i + 1][j]/12 + s[i + 1][j + 1]/12 + s[i - 1][j + 1]/12

    r[r>0] = r[r>0] + gamma
    s = r + nr
    s[s>1] = 1
    f = s >= 1

    hex_centers, _ = create_hex_grid(nx=rows,
                                     face_color=grid_to_color(f),
                                     ny=columns,
                                     do_plot=True)

plt.show()
