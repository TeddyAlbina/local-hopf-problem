"""Mean of the c d^2 cubic coefficient over the equatorial torus (PROOF.md section 6).

Closed form, for every separation and t, with the h_cd of PROOF.md. The
complete-tensor computation is part B of full_metric_check.py.
"""

import sympy as sp

from cd_tensor import cd_profiles
from jets import K, half, matrix, minimized_mixed_cubic, unit, one
from weights_ode import C, S, T, u, rational_in_cosine
from quadratic_coefficient import chart, symprod, t


def coefficient(axis):
    p, q, rotate, flat, g, inverse, d = chart(C, S)
    c = unit(4)
    e = [(1, 0, 0), (0, 1, 0)][axis]
    kz = flat(rotate((0, 0, 1)))
    jz = flat(rotate((0, 0, 1), -1))
    kx = flat(rotate(e))
    xfield = rotate(e, -1)
    jx = flat(xfield)
    norm = sum(v * w for v, w in zip(xfield, jx))
    dot = sum(v * w for v, w in zip(p, q))
    f, g_profile = cd_profiles(t, C)
    f_symbolic = cd_profiles(T, u)[0]
    f_slope = K.from_expr(sp.diff(f_symbolic, u).subs(u, C.as_expr()))
    # The zx prefactor vanishes to first order in heights; zz to second.
    f_jet = (dot - C) * f_slope + f
    zx = q[2] * q[axis] - p[2] * p[axis]
    zz = q[2] * q[2] - p[2] * p[2]

    def tensor(i, j):
        cd = (
            -(zx * f_jet * jz[i] * jz[j] + zz * g_profile * symprod(jx, jz, i, j))
            * half
            * half
        )
        dd = norm * jz[i] * jz[j] * (one / (4 * (1 + 2 * t)))
        return (
            g[i][j]
            + c * kz[i] * kz[j]
            + d * symprod(kx, jz, i, j)
            + d * d * dd
            + 2 * c * d * cd
        )

    value, dd = minimized_mixed_cubic(matrix(tensor), inverse)
    assert not dd
    return value


if __name__ == "__main__":
    mean = (coefficient(0) + coefficient(1)) * half
    assert mean == one / (4 * (1 + 2 * t) ** 2)
    density = sp.factor(rational_in_cosine(mean))
    area = 2 * sp.pi**2 / sp.sqrt(1 + 2 * T)
    integral = sp.factor(area * density)
    assert sp.simplify(integral.subs(T, 1) - sp.pi**2 / (18 * sp.sqrt(3))) == 0
    assert sp.limit(integral, T, 0, dir="+") == sp.pi**2 / 2
    print("Mixed rotation average:", density)
    print("Mixed torus integral:", integral)
    print("Full cubic mean at lambda_c=1/100, lambda_d=1/10:", density / 10000)
