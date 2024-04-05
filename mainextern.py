import ctypes
import os

from hexalattice.hexalattice import *
from ctypes import *

max_iterations = 4000
alpha = 1
beta = 0.3
gamma = 0.0035
n = 401


def grid_to_color(grid):
    result = []
    for i in range(grid.shape[0]):
        for j in range(grid.shape[1]):
            red = 1
            green = 1
            blue = 1
            alpha_c = 1
            if grid[i][j] >= 1:
                red = 0
                green = 0
                blue = 0
                alpha_c = 1
            elif grid[i][j] > beta:
                red = 0
                green = 0
                blue = 0
                alpha_c = 0.7

            result.append((red, green, blue, alpha_c))
    return np.float32(result)


cpp = cdll.LoadLibrary(os.path.abspath("libsnowflake_automata_cpp.so"))

# setup functions
cpp.CreateSnowFlake.restype = c_int64
cpp.CreateSnowFlake.argtypes = [c_int, c_double, c_double, c_double, c_int]

cpp.Iterate.restype = c_bool
cpp.Iterate.argtypes = [c_int64]

cpp.GetGrid.restype = ctypes.POINTER(c_double)
cpp.GetGrid.argtypes = [c_int64]

# Run program
snow_flake = cpp.CreateSnowFlake(max_iterations, alpha, beta, gamma, n)

res = cpp.Iterate(snow_flake)
count = 0
while res:
    res = cpp.Iterate(snow_flake)
    count +=1
    if count%100==0:
        print(count)
grid = cpp.GetGrid(snow_flake)

s = np.zeros((n, n), dtype=float)
for i in range(n):
    for j in range(n):
        s[i][j] = grid[i * n + j]

hex_centers, _ = create_hex_grid(nx=n,
                                 face_color=grid_to_color(s),
                                 ny=n,
                                 do_plot=True, line_width=0)

plt.show()
