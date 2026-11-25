"""Chart, curvature jets and minimized quadratic coefficient (PROOF.md sections 3-5).

chart() places p1=(1,0,0), p2=(cos delta, sin delta, 0) with the great-circle
normal along the third axis. quadratic() computes the fixed-plane and
minimized quadratic coefficients in the graph chart, and asserts that the
plane is flat with zero first variation. The script checks that the d,d
weight 1/(4h) makes the equatorial coefficient vanish. It also reproduces the
Brendle--Hung constant 0.36893... at t=1, and shows that their tensors give
t^3 V^(2) -> -16/25 as t -> 0. The first four jet slots are x, y, u, v; slot
six records the metric parameter s.
"""

from dataclasses import dataclass
from itertools import product

import jets
from jets import (
    K,
    curvature,
    half,
    indices,
    inverse_jet,
    jet,
    matrix,
    one,
    r,
    unit,
    wedge,
)
from jets import (
    cheeger as t,
)
from sympy import Rational as Q
from sympy import factor, simplify, sqrt


def chart(cosine, sine):
    x, y, u, v, _, s = map(unit, range(6))
    a, b = x - y, x + y
    ca, sa = 1 - a * a * half, a
    cb, sb = (
        (1 - b * b * half) * cosine - b * sine,
        (1 - b * b * half) * sine + b * cosine,
    )
    ra, rb = 1 - u * u * half, 1 - v * v * half
    p = [ra * ca, ra * sa, u]
    q = [rb * cb, rb * sb, v]

    def rot(axis, sign=1):
        ax, ay, az = axis
        a1 = az - u * (ca * ax + sa * ay)
        a2 = (az - v * (cb * ax + sb * ay)) * sign
        return [
            (a1 + a2) * half,
            (a2 - a1) * half,
            ra * (sa * ax - ca * ay),
            rb * (sb * ax - cb * ay) * sign,
        ]

    e = [[int(i == j) for j in range(3)] for i in range(3)]
    rotations = list(map(rot, e))
    A = matrix(lambda i, j: jet(0))
    A[0][0] = A[1][1] = (2 + u * u + v * v) * half * half
    A[1][0] = A[0][1] = (v * v - u * u) * half * half
    A[2][2], A[3][3] = 1 - u * u, 1 - v * v
    inv = matrix(lambda i, j: 2 * (A[i][j] + sum(w[i] * w[j] for w in rotations) * t))
    G = matrix(lambda i, j: jet(0))
    G[0][0], G[1][1] = jet(1 / (1 + 2 * t)), jet(1)
    D = 2 * ((1 + t) ** 2 - (t * cosine) ** 2)
    G[2][2] = G[3][3] = jet((1 + t) / D)
    G[2][3] = G[3][2] = jet(-t * cosine / D)
    g = inverse_jet(inv, G, 2)
    leading = matrix(lambda i, j: jet(inv[i][j].get(jets.origin, 0)))

    def flat(w):
        return [sum(g[i][j] * w[j] for j in indices) for i in indices]

    identity = jets.matmul(G, leading)
    assert all(
        identity[i][j] == jet(int(i == j)) for i, j in product(indices, repeat=2)
    )
    return p, q, rot, flat, g, leading, s


def symprod(a, b, i, j):
    return (a[i] * b[j] + a[j] * b[i]) * half


@dataclass
class Quadratic:
    value: object
    tilt_gradient: list
    hessian: list


def quadratic(metric, inv):
    R = curvature(metric, inv)
    E = [[int(i == j) for j in indices] for i in indices]
    U, V, N, M = E
    W = wedge(U, V)
    Z = [0] * 4
    tilts = [(N, Z), (M, Z), (Z, N), (Z, M)]

    def add(A, B):
        return [a + b for a, b in zip(A, B)]

    linear = [add(wedge(a, V), wedge(U, b)) for a, b in tilts]

    def contraction(k, A, B):
        order = (0, 0, 0, 0, 0, k)
        return sum(
            (
                a * b * R[i][j].get(order, 0)
                for i, a in enumerate(A)
                if a
                for j, b in enumerate(B)
                if b
            ),
            K.zero,
        )

    assert contraction(0, W, W) == contraction(1, W, W) == 0

    def hess(i, j):
        a, b = tilts[i]
        e, f = tilts[j]
        return 2 * (
            contraction(0, linear[i], linear[j])
            + contraction(0, W, add(wedge(a, f), wedge(e, b)))
        )

    grad = [2 * contraction(1, W, h) for h in linear]
    H = matrix(hess)
    rows = [row + [v] for row, v in zip(H, grad)]
    for i in indices:
        pivot = rows[i][i]
        rows[i] = [a / pivot for a in rows[i]]
        for j in indices:
            if i != j:
                w = rows[j][i]
                rows[j] = [a - w * b for a, b in zip(rows[j], rows[i])]
    value = contraction(2, W, W) - sum(w * row[-1] for w, row in zip(grad, rows)) * half
    return Quadratic(value, grad, H)


def torus(axis, weight):
    cosine = (1 - r * r) / (1 + r * r)
    sine = 2 * r / (1 + r * r)
    _p, _q, rot, flat, g, inv, s = chart(cosine, sine)
    kx = flat(rot(axis))
    jz = flat(rot((0, 0, 1), -1))
    jx = rot(axis, -1)
    norm = sum(a * b for a, b in zip(jx, flat(jx)))
    return quadratic(
        matrix(
            lambda i, j: (
                g[i][j]
                + s * symprod(kx, jz, i, j)
                + s * s * norm * jz[i] * jz[j] * weight
            )
        ),
        inv,
    )


def meridian(cosine, sine, lc=0, ld=0, *, b_intercept=one * 6 / 5, b_slope=0):
    p, q, rot, flat, g, inv, s = chart(cosine, sine)
    # Global (e_z,e_x,e_y) is the local positively oriented coordinate basis.
    # At separation pi/2 the point is (global e_z, global e_x).
    plus = [a + b for a, b in zip(p, q)]
    minus = [b - a for a, b in zip(p, q)]
    dot = sum(a * b for a, b in zip(p, q))
    kp, km = flat(rot(plus)), flat(rot(minus))
    kx, jx = flat(rot((0, 1, 0))), flat(rot((0, 1, 0), -1))
    ky, jy = flat(rot((0, 0, 1))), flat(rot((0, 0, 1), -1))
    kz, jz = flat(rot((1, 0, 0))), flat(rot((1, 0, 0), -1))

    def metric(i, j):
        a = (2 - dot) * symprod(kp, kx, i, j) * (one * 6 / 5) - 2 * (1 - dot) * (
            2 + dot
        ) * symprod(km, jx, i, j)
        b = (2 + dot) * symprod(km, ky, i, j) * (b_intercept + dot * b_slope) - 2 * (
            1 + dot
        ) * (2 - dot) * symprod(kp, jy, i, j)
        c = kz[i] * kz[j]
        d = symprod(kx, jz, i, j)
        return g[i][j] + s * (a + b + lc * c + ld * d)

    # Only the b,b term in the second-order tensor has nonzero pullback
    # to this great-circle torus. Its linearized curvature is added below.
    return quadratic(matrix(metric), inv)


if __name__ == "__main__":
    ts = t.as_expr()
    weight = 1 / (4 * (1 + 2 * t))
    cosine = (1 - r * r) / (1 + r * r)
    D = (1 + t) ** 2 - (t * cosine) ** 2
    for axis in [(1, 0, 0), (0, 1, 0), (1, 1, 0)]:
        # Polarization of a homogeneous quadratic proves the identity
        # for every horizontal axis at every nondegenerate separation.
        assert torus(axis, weight).value == 0
    uncorrected = torus((1, 0, 0), one / 12).value
    assert uncorrected == (1 - t) * cosine**2 / (6 * D)
    print("With the weight 1/(4h) the equatorial d,d coefficient vanishes identically.")

    sample = meridian(K.zero, one)
    base = sample.value
    # These six pairs determine every polynomial of total degree <= 2
    # in lambda_c,lambda_d. Check the full first-order tensor in each case.
    for lc, ld in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1)]:
        assert meridian(K.zero, one, lc, ld).value == base

    # At delta=pi/2 odd harmonics do not contribute to f_bb''(delta).
    # The geometric sum from the nonpolynomial tail is -k^2/(1+k^2).
    k = 4 - sqrt(15)
    fpp = simplify(
        4 * Q(5579141, 153600)
        - 16 * Q(48211, 307200)
        + 36 * Q(25, 165888)
        - Q(42025, 6) * sqrt(Q(5, 3)) * k**2 / (1 + k**2)
    )
    V = base.as_expr() - 18 * fpp / (1 + 2 * ts) ** 2
    at_one = simplify(V.subs(ts, 1))
    at_half = simplify(V.subs(ts, Q(1, 2)))
    assert simplify(at_one - (-Q(3472117, 384) + Q(42025, 18) * sqrt(15))) == 0
    assert at_half.is_negative
    assert simplify((ts**3 * V).limit(ts, 0)) == -Q(16, 25)
    print(
        "Metric quadratic and minimizing-plane contribution =", factor(base.as_expr())
    )
    print("f_bb''(pi/2) (longitude separation derivative) =", fpp)
    print("V at 1 =", at_one, "~", at_one.evalf(12))
    print("V at 1/2 =", at_half, "~", at_half.evalf(12))
    print("limit of t^3 V as t -> 0+ =", (ts**3 * V).limit(ts, 0))

    # Tilting the longitude direction on sphere 1 alone moves through
    # product-flat mixed planes at t=0 at every point of the meridian torus.
    tilt = (1, 0, -1, 0)
    gradient = sum(a * b for a, b in zip(sample.tilt_gradient, tilt))
    hessian = sum(
        sample.hessian[i][j] * tilt[i] * tilt[j] for i, j in product(indices, repeat=2)
    )
    assert gradient.as_expr().subs(ts, 0) == Q(8, 5)
    assert hessian.as_expr().subs(ts, 0) == 0
    print("Gradient on the product-flat tilt at t=0 = 8/5; Hessian = 0.")
