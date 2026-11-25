"""Normal corrections M_t, N_t in h_ad and h_bd (DERIVATIONS.md, section 4).

M_t, N_t multiply the J-tensor terms of h_ad and h_bd in PROOF.md section 2.
Checks that the a,d and b,d coefficients vanish at every regular point and
every t, that the d,d blocks have uniform limits of size <= 1/4, and that the
Brendle--Hung form of h_cd gives a nonzero normal slope at the equator.
"""

from functools import lru_cache

import sympy as sp
from jets import unit
from weights_ode import C, S, T
from quadratic_coefficient import (
    K,
    chart,
    curvature,
    half,
    matrix,
    one,
    quadratic,
    symprod,
    t,
)
from first_order_weights import P, Q, gradient, regular_on_circle

# Each row lists a polynomial in u, in ascending powers. Rows are in ascending powers of t.
M_ROWS = [
    [-134, 32, 124, -48],
    [-36, 116, 192, -24, -64],
    [-74, 159, 617, -19, -175, -36],
    [-212, -322, 306, 394, -94, -72],
]
N_ROWS = [
    [-242, -16, 268, 336],
    [-812, -1566, 1146, 1660, 456],
    [-1056, -3013, -543, 2449, 1423, 124],
    [-618, -1957, -1158, 792, 812, 165, -36],
    [-212, -534, -16, 700, 300, -166, -72],
]


def profiles(parameter, dot):
    """Return M_t(dot), N_t(dot) for the unnormalized mixed coefficient."""
    a = 1 + parameter * (1 + dot)

    def poly(rows):
        return sum(
            parameter**i * sum(k * dot**j for j, k in enumerate(row))
            for i, row in enumerate(rows)
        )

    return (
        poly(M_ROWS) / (80 * (1 + 2 * parameter) ** 2 * a),
        poly(N_ROWS) / (80 * (1 + 2 * parameter) ** 2 * a * a),
    )


def solve(H, g):
    rows = [row[:] + [z] for row, z in zip(H, g)]
    for i in range(4):
        v = rows[i][i]
        rows[i] = [x / v for x in rows[i]]
        for j in range(4):
            if j != i:
                v = rows[j][i]
                rows[j] = [a - v * b for a, b in zip(rows[j], rows[i])]
    return [row[-1] for row in rows]


@lru_cache
def geometry():
    p, q, rot, flat, g, inv, s = chart(C, S)
    return p, q, rot, flat, g, inv, s, quadratic(g, inv).hessian


def ad(axis, kind=1, channel="a"):
    p, q, rot, flat, g, inv, s, H = geometry()
    dot = sum(x * y for x, y in zip(p, q))
    minus = [y - x for x, y in zip(p, q)]
    plus = [y + x for x, y in zip(p, q)]
    km, kp = flat(rot(minus)), flat(rot(plus))
    kn, jn = flat(rot((0, 0, 1))), flat(rot((0, 0, 1), -1))
    if channel == "a":
        weight, weight2 = P(-dot), Q(-dot)
        v1, v2 = kp, km
        ex = minus
        pol = -92
    else:
        weight, weight2 = P(dot), Q(dot)
        v1, v2 = km, kp
        ex = plus
        pol = 92
    kw, jw = flat(rot(axis)), flat(rot(axis, -1))
    f = ((pol * dot - 406) + 18 * (2 * dot * dot - 1)) * (one / 400)
    g2 = ((-pol * 77 * one / 46) * dot - 823 + 9 * (2 * dot * dot - 1)) * (one / 200)
    proj = sum(x * y for x, y in zip(ex, axis))
    c = unit(4)

    def metric(i, j):
        A = weight * symprod(v1, kn, i, j) + weight2 * symprod(v2, jn, i, j)
        if kind == 1:
            D = symprod(kw, jn, i, j)
            coeff = f
        else:
            D = symprod(kn, jw, i, j)
            coeff = g2
        return g[i][j] + s * (D + c * A) + 2 * c * s * s * coeff * proj * kn[i] * kn[j]

    R = curvature(matrix(metric), inv)
    gradD = gradient(R, (0, 0, 0, 0, 0, 1))
    gradA = gradient(R, (0, 0, 0, 0, 1, 1))
    return R[0][0].get((0, 0, 0, 0, 1, 2), K.zero) - sum(
        x * y for x, y in zip(gradA, solve(H, gradD))
    )


def forms():
    _p, _q, rot, flat, _g, _inv, _s, _H = geometry()
    R0 = flat(rot((0, 0, 1)))
    J0 = flat(rot((0, 0, 1), -1))
    out = []
    for k in (1, 2):
        for ax in [(1, 0, 0), (0, 1, 0)]:
            kw, jw = flat(rot(ax)), flat(rot(ax, -1))
            v, w = (kw, J0) if k == 1 else (R0, jw)
            out.append(matrix(lambda i, j, v=v, w=w: symprod(v, w, i, j)))
    return out


def dd(i, j, fs):
    _p, _q, _rot, _flat, g, inv, s, H = geometry()
    F, G = fs[i], fs[j]
    c = unit(4)
    R = curvature(matrix(lambda a, b: g[a][b] + s * (F[a][b] + c * G[a][b])), inv)
    r0 = gradient(R, (0, 0, 0, 0, 0, 1))
    r1 = gradient(R, (0, 0, 0, 0, 1, 1))
    return (
        R[0][0].get((0, 0, 0, 0, 1, 2), K.zero)
        - sum(x * y for x, y in zip(r0, solve(H, r1)))
    ) * half


def cd():
    _p, _q, rot, flat, g, inv, s = chart(K.zero, one)
    eps = unit(4)
    kx = flat(rot((1, 0, eps)))
    kz = flat(rot((-eps, 0, 1)))
    jz = flat(rot((-eps, 0, 1), -1))

    def metric(i, j):
        return g[i][j] + s * (kz[i] * kz[j] + symprod(kx, jz, i, j))

    R = curvature(matrix(metric), inv)
    g0 = gradient(R, (0, 0, 0, 0, 0, 1))
    g1 = gradient(R, (0, 0, 0, 0, 1, 1))
    H = quadratic(g, inv).hessian
    raw = R[0][0].get((0, 0, 0, 0, 1, 2), K.zero) - sum(
        x * y for x, y in zip(g1, solve(H, g0))
    )
    # On the torus the Brendle--Hung h_cd contributes -(F(0)-F_double_prime(0))/h^2,
    # with F(u)=G(sqrt(1-u^2)/sqrt(3)), where F(0)-F_double_prime(0)=3/4.
    return raw - 3 * one / (4 * (1 + 2 * t) ** 2)


def check():
    for ch, sgn in [("a", -1), ("b", 1)]:
        reference = profiles(t, C if ch == "a" else -C)
        for k, profile in enumerate(reference, 1):
            # Two in-circle test directions suffice by linearity. The
            # vectors p2-p1, p1+p2 have components (C+sgn,S).
            assert ad((1, 0, 0), k, ch) == (C + sgn) * profile
            assert ad((0, 1, 0), k, ch) == S * profile
            regular_on_circle(profile)
    print("The a,d and b,d coefficients vanish with the normal corrections.")
    c, s = C.as_expr(), S.as_expr()
    limits = {
        (0, 2): -c * c / 4,
        (0, 3): -c * s / 4,
        (1, 2): -c * s / 4,
        (1, 3): c * c / 4,
        (2, 2): c * c / 4,
        (2, 3): c * s / 4,
        (3, 3): -c * c / 4,
    }
    fs = forms()
    for (i, j), target in limits.items():
        val = dd(i, j, fs)
        regular_on_circle(val)
        assert sp.cancel(val.as_expr().subs(T, 0) - target) == 0
    print("All remaining dd matrix entries have uniform limits of magnitude <= 1/4.")
    witness = sp.factor(cd().as_expr())
    assert (
        sp.cancel(
            witness - (1 - T) * (3 * T + 1) / (4 * (1 + T) ** 2 * (1 + 2 * T) ** 2)
        )
        == 0
    )
    print("With the Brendle--Hung h_cd, the mixed normal slope is", witness)


if __name__ == "__main__":
    check()
