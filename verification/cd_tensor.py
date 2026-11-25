"""The mixed tensor h_cd (PROOF.md section 2; DERIVATIONS.md section 5).

Checks that h_cd cancels the c,d coefficient at every regular point and
every t>0, that the c,c coefficient vanishes, and the t=0 bound on the
coefficient linear in lambda_c. All coefficients refer to the minimized
curvature numerator in the U,V frame.
"""

import sympy as sp

from weights_ode import C, S, T, u
from quadratic_coefficient import K, matrix, one, symprod, t
from normal_corrections import dd, forms, geometry
from first_order_weights import P, Q, regular_on_circle


def cd_profiles(parameter, dot):
    """Coefficients of the two coordinate products before correction."""
    h = 1 + 2 * parameter
    a = 1 + parameter
    d = a**2 - parameter**2 * dot**2
    f = (a**2 + parameter**2 * dot**2) / (h * d**2)
    numerator = h**2 + 6 * h * parameter**2 * dot**2 + parameter**4 * (dot**4 - 1)
    g = -a * numerator / (h**2 * d**3)
    return f, g


def first_variations():
    p, q, rot, flat, *_ = geometry()
    kn, jn = flat(rot((0, 0, 1))), flat(rot((0, 0, 1), -1))
    tensors = forms()  # D1(e1), D1(e2), D2(e1), D2(e2)
    for axis in ((1, 0, 0), (0, 1, 0)):
        kw = flat(rot(axis))
        tensors.append(matrix(lambda i, j: symprod(kn, kw, i, j)))
    dot = sum(x * y for x, y in zip(p, q))
    plus = flat(rot([x + y for x, y in zip(p, q)]))
    minus = flat(rot([y - x for x, y in zip(p, q)]))
    for arg, first, second in ((-dot, plus, minus), (dot, minus, plus)):
        wp, wq = P(arg), Q(arg)
        tensors.append(
            matrix(
                lambda i, j: (
                    wp * symprod(first, kn, i, j) + wq * symprod(second, jn, i, j)
                )
            )
        )
    return tensors


def check_mixed_cd(tensors):
    f, g = cd_profiles(t, C)
    # Delta = p2 tensor p2 - p1 tensor p1, in the in-circle frame.
    delta = ((-(S**2), C * S), (C * S, S**2))
    for i in range(2):
        for j in range(2):
            assert 4 * dd(4 + i, j, tensors) == f * delta[i][j]
            # Both inputs to the second block are z_parallel.
            symmetric = 2 * (dd(4 + i, 2 + j, tensors) + dd(4 + j, 2 + i, tensors))
            assert symmetric == g * delta[i][j]

    # Along U the two points undergo the same planar rotation. Check its
    # second derivative on every basis component, including the diagonal.
    rot90 = sp.Matrix([[0, -1], [1, 0]])
    d = sp.Matrix([[x.as_expr() for x in row] for row in delta])
    first = rot90 * d - d * rot90
    assert (rot90 * first - first * rot90 + 4 * d).applyfunc(sp.cancel) == sp.zeros(2)
    print("Mixed c,d coefficient cancels identically for every regular point and t>0.")


def check_c_bounds(tensors):
    p, q, rot, flat, *_ = geometry()
    dot = sum(x * y for x, y in zip(p, q))
    key = (0, 2, 0, 0, 0, 0)  # coefficient of y^2, where V = d/dy
    h2 = (1 + 2 * t) ** 2
    axes = ((1, 0, 0), (0, 1, 0))
    for i, axis in enumerate(axes):
        for j in range(i, 2):
            norm = sum(x * y for x, y in zip(rot(axis), flat(rot(axes[j]))))
            # Linear curvature of h_cc is -d_yy h_cc(U,U)/2.
            assert 4 * dd(4 + i, 4 + j, tensors) == norm.get(key, K.zero) / h2

    # The two vectors below multiply n_x n_z and n_y n_z respectively,
    # contracted in each case against the projection of e_z into the circle.
    def limit_profile(v):
        return (-168 * v**3 + 154 * v**2 + 148 * v - 99) / 20

    c = C.as_expr()
    for j, sign in ((6, 1), (7, -1)):
        weight = (-8 - 14 * sign * dot + 6 * dot * dot) * (one / 20)
        for i in range(2):
            inherited = weight * (q[i] + sign * p[i])
            value = 4 * dd(4 + i, j, tensors) - 2 * inherited.get(key, K.zero) / h2
            regular_on_circle(value)
            coord = c + sign if i == 0 else S.as_expr()
            assert (
                sp.cancel(value.as_expr().subs(T, 0) - limit_profile(sign * c) * coord)
                == 0
            )

    # A strict bound on the squared Frobenius norm of the limiting pair
    # persists uniformly as t -> 0 by regularity on the compact half-angle line.
    squared = 2 * (1 + u) * limit_profile(u) ** 2 + 2 * (1 - u) * limit_profile(-u) ** 2
    gap = sp.Poly(100 - squared, u)
    assert gap.eval(-1) > 0 and gap.count_roots(-1, 1) == 0
    print(
        "The c,c term vanishes exactly; the combined a,c and b,c loss is <= 10 |lambda_c| rho."
    )


if __name__ == "__main__":
    tensors = first_variations()
    check_mixed_cd(tensors)
    check_c_bounds(tensors)
