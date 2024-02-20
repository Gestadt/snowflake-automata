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
n = 51 # must be uneven

s = np.ones((n, n), dtype=float)*beta
f = np.zeros((n, n), dtype=float)
s[n // 2][n // 2] = 1
'''
s[rows//2][columns//2-1] = 0.7
s[rows//2][columns//2+1] = 0.7
s[rows//2-1][columns//2] = 0.7
s[rows//2+1][columns//2] = 0.7
s[rows//2+1][columns//2+1] = 0.7
s[rows//2-1][columns//2+1] = 0.7
'''

r = np.zeros((n, n), dtype=float)
nr = np.zeros((n, n), dtype=float)
av = np.zeros((n,n), dtype=float)
for k in range(100):
    r = r * 0
    nr = nr * 0
    av = av * 0
    for i in range(n):
        for j in range(n):
            if s[i][j] >= 1 or s[i][(j - 1)%n] >= 1 or s[i][(j + 1)%n] >= 1 or s[(i - 1)%n][j] >= 1 or s[(i + 1)%n][j] >= 1 or s[(i + 1)%n][(j + 1)%n] >= 1 or s[(i - 1)%n][(j + 1)%n] >= 1:
                r[i][j] = s[i][j] + gamma
            else:
                nr[i][j] = s[i][j]

    for i in range(n):
        for j in range(n):
            av[i][j] = (nr[i][j]/2 +
                        nr[i][(j - 1)%n]/12 +
                        nr[i][(j + 1)%n]/12 +
                        nr[(i - 1)%n][j]/12 +
                        nr[(i + 1)%n][j]/12 +
                        nr[(i + 1)%n][(j + 1)%n]/12 +
                        nr[(i - 1)%n][(j + 1)%n]/12)

    s = r + av
    f = s >= 1

hex_centers, _ = create_hex_grid(nx=n,
                                 face_color=grid_to_color(f),
                                 ny=n,
                                 do_plot=True)

plt.show()
