"""Test D - cellular flow (Maxey 1987): the parameters of the case, used by make_case.py and plot_D.py."""
import math

# air (as all other tests)
rho_f, nu = 1.2, 1.5e-5
mu = rho_f * nu

# cellular flow:  u_x =  U0 sin(pi x / L) cos(pi y / L)
#                 u_y = -U0 cos(pi x / L) sin(pi y / L)
# domain 2L x 2L = 2 x 2 vortices, 1 cell thick in z
L = 0.01                 # vortex size [m]
U0 = 0.1                 # flow speed [m/s]
tau_f = L / U0           # flow time scale [s]
N = 32                   # mesh cells per vortex (in x and in y)
nx = ny = 2 * N
dx = 2 * L / nx          # cell size [m]
dz = dx                  # thickness of the single cell in z [m]

# particles: solid, one size per Stokes number St = tau_p / tau_f
rho_p = 2500.0
ST = [0.1, 0.5, 1.0, 2.0]
n_side = 10              # n_side x n_side particles per size, on a grid, released at rest

# time
dt = 1e-4
end_time = 1.0           # 10 tau_f
write_interval = 50      # write every 50 steps


def diameter(St):
    """d from St: tau_p = rho_p d^2 / (18 mu) = St tau_f"""
    return math.sqrt(18 * mu * St * tau_f / rho_p)


def flow(x, y):
    return (U0 * math.sin(math.pi * x / L) * math.cos(math.pi * y / L),
            -U0 * math.cos(math.pi * x / L) * math.sin(math.pi * y / L))


def start_positions():
    """n_side x n_side grid inside the domain, shifted by 0.3 dx so that no particle starts on a
    cell face (on a face both neighbouring cells are correct, and the test could not decide)"""
    h = 2 * L / n_side
    s = 0.3 * dx
    return [((i + 0.5) * h + s, (j + 0.5) * h + s, 0.5 * dz) for j in range(n_side) for i in range(n_side)]
