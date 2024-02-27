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


max_iterations = 200
alpha = 1
beta = 0.5
gamma = 0.001
n = 51

s = np.ones((n, n), dtype=float) * beta
f = np.zeros((n, n), dtype=float)
mid = int((n - 1) / 2)
s[mid][mid] = 1

c = np.zeros((n,n), dtype=bool) # classification (receptive (1) / non-receptive (2))
r = np.zeros((n, n), dtype=float) # receptive grid
nr = np.zeros((n, n), dtype=float) # non-receptive grid
av = np.zeros((n, n), dtype=float) # average (diffusion) of non-receptive grid
for k in range(max_iterations):
    r = r * 0
    nr = nr * 0
    av = av * 0
    for i in range(n):
        for j in range(n):
            if s[i][j] >= 1:
                c[i][j] = 1
                c[(i + 1) % n][j] = 1
                c[(i - 1) % n][j] = 1
                c[i][(j - 1) % n] = 1
                c[i][(j + 1) % n] = 1
                if i % 2 == 0:
                    c[(i + 1) % n][(j - 1) % n] = 1
                    c[(i - 1) % n][(j - 1) % n] = 1
                else:
                    c[(i + 1) % n][(j + 1) % n] = 1
                    c[(i - 1) % n][(j + 1) % n] = 1
    for i in range(n):
        for j in range(n):
            if c[i][j]:
                r[i][j] = s[i][j]+gamma
            else:
                nr[i][j] = s[i][j]

    for i in range(n):
        for j in range(n):
            if i % 2 == 0:
                av[i][j] = (nr[i][j] * (1 - alpha * 0.5) +
                            nr[i][(j - 1) % n] * (alpha / 12) +
                            nr[i][(j + 1) % n] * (alpha / 12) +
                            nr[(i - 1) % n][j] * (alpha / 12) +
                            nr[(i + 1) % n][j] * (alpha / 12) +
                            nr[(i + 1) % n][(j - 1) % n] * (alpha / 12) +
                            nr[(i - 1) % n][(j - 1) % n] * (alpha / 12))
            else:
                av[i][j] = (nr[i][j] * (1 - alpha * 0.5) +
                            nr[i][(j - 1) % n] * (alpha / 12) +
                            nr[i][(j + 1) % n] * (alpha / 12) +
                            nr[(i - 1) % n][j] * (alpha / 12) +
                            nr[(i + 1) % n][j] * (alpha / 12) +
                            nr[(i + 1) % n][(j + 1) % n] * (alpha / 12) +
                            nr[(i - 1) % n][(j + 1) % n] * (alpha / 12))

    s = r + av
    f[s>0.6] = s[s>0.6]

hex_centers, _ = create_hex_grid(nx=n,
                                 face_color=grid_to_color(f/np.max(f)),
                                 ny=n,
                                 do_plot=True)

plt.show()
