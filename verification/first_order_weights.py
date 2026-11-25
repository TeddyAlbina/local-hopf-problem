"""First-order weights P_t, Q (PROOF.md section 2; DERIVATIONS.md section 2).

Checks the structure of the tilt Hessian as t -> 0, that the soft tilt
gradients of the a, b and d channels vanish to the required order, and that
the limiting coefficient F_0 satisfies F_0 > 1 on [-1,1] (exact root count).
Also evaluates the a,d slope at the equator without the normal corrections
(DERIVATIONS.md, section 4).
"""

import sympy as sp
from single_separation_weights import fbb_second_derivative
from jets import half, unit
from weights_ode import C, S, T, rational_in_cosine, u
from quadratic_coefficient import (
    K,
    chart,
    curvature,
    matrix,
    one,
    quadratic,
    r,
    symprod,
    t,
)


def geom():
    *_, g, inv, _s = chart(C, S)
    H = quadratic(g, inv).hessian
    H = sp.Matrix([[rational_in_cosine(z) for z in row] for row in H])
    return H


def blocks(H, sg):
    # hard / soft x,y by sign under interchange, normalized unit basis
    v = sp.Matrix([1, sg, 1, -sg]) / 2
    w = sp.Matrix([1, sg, -1, sg]) / 2
    A = sp.factor(v.dot(H * v))
    B = sp.factor(v.dot(H * w))
    D = sp.factor(w.dot(H * w))
    return A, B, D


def P(v):
    return (one * 16 / 5 + v * half - v * v * (2 * one / 5)) * (
        1 - (1 + v) * (1 + v) * t * t
    )


def Q(v):
    return 2 * v * v - v * (7 * one / 10) - 4


def gradient(R, key):
    return [
        2 * sign * R[0][ix].get(key, K.zero)
        for ix, sign in [(3, -1), (4, -1), (1, 1), (2, 1)]
    ]


def axis_data(d_axis=None):
    p, q, rot, flat, g, inv, s = chart(C, S)
    dot = sum(x * y for x, y in zip(p, q))
    normal = (0, 0, 1)
    kn, jn = flat(rot(normal)), flat(rot(normal, -1))
    plus = flat(rot([x + y for x, y in zip(p, q)]))
    minus = flat(rot([y - x for x, y in zip(p, q)]))
    first, second = P(dot), Q(dot)
    if d_axis:
        jz = flat(rot(d_axis, -1))

    def metric(i, j):
        if d_axis:
            b = symprod(kn, jz, i, j)
        else:
            b = first * symprod(minus, kn, i, j) + second * symprod(plus, jn, i, j)
        return g[i][j] + s * b

    R = curvature(matrix(metric), inv)
    return gradient(R, (0, 0, 0, 0, 0, 1)), R[0][0].get((0, 0, 0, 0, 0, 2), K.zero)


def same(expr):
    assert sp.cancel(expr) == 0


def unit_vectors(sig):
    return sp.Matrix([1, sig, 1, -sig]) / 2, sp.Matrix([1, sig, -1, sig]) / 2


def regular_on_circle(x):
    """No poles at t=0 on the compactified real half-angle line."""
    if not x:
        return
    n, d = x.numer, x.denom
    rr = r.as_expr()
    denom_zero = sp.Poly(d.as_expr().subs(T, 0), rr)
    assert not denom_zero.is_zero and denom_zero.count_roots(-sp.oo, sp.oo) == 0
    assert n.degree(1) <= d.degree(1) == denom_zero.degree()


def certify_schur():
    H = geom()
    plus = sp.Matrix.hstack(*unit_vectors(1))
    minus = sp.Matrix.hstack(*unit_vectors(-1))
    assert (plus.T * H * minus).applyfunc(sp.cancel) == sp.zeros(2)
    b_grad, raw = axis_data()
    d_data = [axis_data(j)[0] for j in [(1, 0, 0), (0, 1, 0)]]
    loss0 = 0
    for sig in [1, -1]:
        A, B, D = blocks(H, sig)
        ratio = sp.factor(B / A)
        schur = sp.factor(D - B * B / A)
        same(A.subs(T, 0) - 2)
        same(ratio.subs(T, 0))
        same(ratio.diff(T).subs(T, 0) - (1 - sig * u))
        scaled = sp.factor(schur / (T**3 * (1 - sig * u)))
        same(scaled.subs(T, 0) - 2)
        for expr in (scaled, A):
            denom = sp.Poly(sp.denom(expr).subs(T, 0), u)
            assert not denom.is_zero and denom.count_roots(-1, 1) == 0
        k_ratio = K.from_expr(sp.cancel(ratio.subs(u, C.as_expr())))
        for grad in [b_grad, *d_data]:
            hard = (grad[0] + sig * grad[1] + grad[2] - sig * grad[3]) / 2
            soft = (grad[0] + sig * grad[1] - grad[2] + sig * grad[3]) / 2
            residual = soft - k_ratio * hard
            regular_on_circle(residual**2 / (t**4 * (1 - sig * C)))
            regular_on_circle(hard)
            if grad is b_grad:
                hard_sq = rational_in_cosine(hard**2)
                loss0 += sp.cancel(hard_sq.subs(T, 0) / 4)
    regular_on_circle(raw)
    geometric = sp.factor(
        rational_in_cosine(K.from_expr(raw.as_expr().subs(T, 0))) - loss0
    )
    expected = (
        1152 * u**6
        - 400 * u**5
        - 10400 * u**4
        - 2956 * u**3
        + 18544 * u**2
        + 6920 * u
        - 4643
    ) / 400
    same(geometric - expected)
    limit = sp.apart(
        geometric - 18 * fbb_second_derivative(u), u, extension=sp.sqrt(15)
    )
    lower_sqrt = sp.Rational(38729833462, 10**10)
    assert lower_sqrt**2 < 15
    same(limit.coeff(sp.sqrt(15)) - sp.Rational(42025, 2))
    lower = sp.Poly(
        sp.cancel((4 - u) * (limit.xreplace({sp.sqrt(15): lower_sqrt}) - 1)), u
    )
    assert lower.eval(-1) > 0 and lower.count_roots(-1, 1) == 0
    print("F_0(u)>1 on [-1,1]. Geometric part:", geometric)


def cross():
    ca, sa = K.zero, one
    p, q, rot, flat, g, inv, s = chart(ca, sa)
    # pair (e1,e2); orient global x along local e1 (or e2), z normal.
    eps = unit(4)
    xa, za = (1, 0, eps), (-eps, 0, 1)
    dot = sum(x * y for x, y in zip(p, q))
    kp = flat(rot([x + y for x, y in zip(p, q)]))
    km = flat(rot([y - x for x, y in zip(p, q)]))
    kx, jx = flat(rot(xa)), flat(rot(xa, -1))
    kz, jz = flat(rot(za)), flat(rot(za, -1))
    # symmetric trial
    da = P(-dot)
    qa = Q(-dot)
    px = sum((y - x) * z for x, y, z in zip(p, q, xa))
    pz = sum((y - x) * z for x, y, z in zip(p, q, za))
    # Paper defines h_ad, enters with factor two
    cad = px * (-406 - 92 * dot + 18 * (2 * dot * dot - 1)) * (one / 400)
    cad2 = pz * (-823 + 154 * dot + 9 * (2 * dot * dot - 1)) * (one / 200)

    def metric(i, j):
        ha = da * symprod(kp, kx, i, j) + qa * symprod(km, jx, i, j)
        hd = symprod(kx, jz, i, j)
        had = cad * symprod(kx, kz, i, j) + cad2 * kx[i] * kx[j]
        return g[i][j] + s * (ha + hd) + 2 * s * s * had

    R = curvature(matrix(metric), inv)
    # For a^2,d^2: eta derivative zero by reflection. Hessian independent of eta.
    grads = [
        [
            2 * sgn * R[0][ix].get((0, 0, 0, 0, i, 1), K.zero)
            for ix, sgn in [(3, -1), (4, -1), (1, 1), (2, 1)]
        ]
        for i in range(2)
    ]
    G = [sp.Matrix([v.as_expr() for v in gr]) for gr in grads]
    H = sp.Matrix([[v.as_expr() for v in row] for row in quadratic(g, inv).hessian])
    X = R[0][0].get((0, 0, 0, 0, 1, 2), K.zero).as_expr() - (G[0].T * H.inv() * G[1])[0]
    slope = sp.factor(X)
    same(
        slope
        - (106 * T**3 + 37 * T**2 + 18 * T + 67) / (40 * (1 + T) * (1 + 2 * T) ** 2)
    )
    print("Linear a,d coefficient at the equator:", slope)


if __name__ == "__main__":
    certify_schur()
    cross()
