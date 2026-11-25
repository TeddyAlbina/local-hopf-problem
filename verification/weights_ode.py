"""Product-flat tilt equation for the b weights (DERIVATIONS.md, section 2).

Derives the ODE that the weights P, Q must satisfy at t=0 and its regular
perturbation solution, and evaluates the residual soft-tilt gradient of the
Brendle--Hung weights. The variable u denotes p1.p2; derivatives in u hold
t fixed.
"""

from quadratic_coefficient import (
    K,
    chart,
    curvature,
    half,
    matrix,
    meridian,
    one,
    r,
    symprod,
    t,
)
from sympy import Matrix, Poly, Rational, cancel, diff, expand, factor, symbols

u = symbols("u")
T = t.as_expr()
C = (1 - r * r) / (1 + r * r)
S = 2 * r / (1 + r * r)
a = 1 + T * (1 + u)
m = 1 + T * (1 - u)
h = 1 + 2 * T
q = -2 * (1 + u) * (2 - u)
drift = (3 - 5 * u) / 2 - T**2 * (1 - u * u) * ((3 * T + 2) * u - T - 1) / (
    2 * h * m * a
)
potential = (
    (4 * T + 3) / (4 * a)
    - (9 * T + 6) / (4 * h)
    - (6 * T * T + 3 * T - 1) / (4 * h * m)
)
source = (1 - u) * diff(q, u) / (2 * a) - (
    2 + T * (1 - u) * (4 + 3 * T) + T**3 * (1 - u * u)
) * q / (4 * m * m * a)


def gradients(j0, j1, j2, with_source):
    """Curvature jet for Y=(2+u)*W/a on the great circle with normal e_y."""
    p, p2, rot, flat, g, inv, s = chart(C, S)
    dot = sum(v * w for v, w in zip(p, p2))
    du = dot - C
    weight = (1 + (1 + dot) * t) * (j0 + du * j1 + du * du * (half * j2))
    kp = flat(rot([v + w for v, w in zip(p, p2)]))
    km = flat(rot([w - v for v, w in zip(p, p2)]))
    ky = flat(rot((0, 0, 1)))
    jy = flat(rot((0, 0, 1), -1))

    def metric(i, j):
        perturbation = weight * symprod(km, ky, i, j)
        if with_source:
            perturbation -= 2 * (1 + dot) * (2 - dot) * symprod(kp, jy, i, j)
        return g[i][j] + s * perturbation

    R = curvature(matrix(metric), inv)
    key = (0, 0, 0, 0, 0, 1)
    return [
        2 * sign * R[0][i].get(key, K.zero)
        for sign, i in [(-1, 3), (-1, 4), (1, 1), (1, 2)]
    ]


def rational_in_cosine(expr):
    """After dividing by sin(delta), all powers of r=tan(delta/2) are even."""
    n, d = expr.as_expr().as_numer_denom()

    def convert(poly):
        terms = Poly(poly, r.as_expr()).terms()
        assert all(i % 2 == 0 for (i,), _ in terms)
        return sum(v * ((1 - u) / (1 + u)) ** (i // 2) for (i,), v in terms)

    return cancel(convert(n) / convert(d))


def series_solution(order):
    """Compute Y_n by inverting the leading operator on polynomials."""
    import sympy as sp

    # Clear rational denominators before extracting powers.
    coeffs = [1 - u * u, drift, potential, source]
    denom = 4 * h * m * m * a
    coeffs = [Poly(cancel(c * denom), T) for c in coeffs]

    def apply(k, y):
        return sum(
            c.nth(k) * v for c, v in zip(coeffs[:3], (diff(y, u, 2), diff(y, u), y))
        )

    out = []
    for n in range(order + 1):
        residual = expand(
            coeffs[3].nth(n) + sum(apply(n - i, y) for i, y in enumerate(out))
        )
        degree = Poly(residual, u).degree()
        names = sp.symbols(f"c:{degree + 1}")
        trial = sum(c * u**i for i, c in enumerate(names))
        eq = Poly(expand(apply(0, trial) + residual), u)
        solution = sp.solve(eq.all_coeffs(), names, dict=True)
        assert len(solution) == 1
        out.append(trial.subs(solution[0]))
    return out


if __name__ == "__main__":
    jets = [(1, 0, 0, False), (0, 1, 0, False), (0, 0, 1, False), (0, 0, 0, True)]
    columns = [gradients(*args) for args in jets]
    for gr, coef in zip(columns, [potential, drift, 1 - u * u, source]):
        assert gr[0] == gr[1] and gr[2] == -gr[3]
        assert cancel(rational_in_cosine((gr[0] - gr[2]) * (1 + 2 * t) / S) - coef) == 0
    assert cancel(drift.subs(u, -1) - 4) == cancel(drift.subs(u, 1) + 1) == 0
    y0, y1, y2 = series_solution(2)
    print("Y_0 =", factor(y0))
    print("Y_1 =", factor(y1))
    print("Y_2 =", factor(y2))
    w0 = factor(y0 / (2 + u))
    w1 = factor((y1 + (1 + u) * y0) / (2 + u))
    print("W_0 =", w0)
    print("W_1 =", w1)
    # Evaluate the remaining tilt component at the calibrated point. On this
    # meridian the a,c,d first-order pieces are reflection-even, tangentially
    # zero, and have zero tilt gradients. The b tensor is reflection-odd.
    data = [y0.subs(u, 0), diff(y0, u).subs(u, 0), diff(y0, u, 2).subs(u, 0), 1]
    gamma = sum(
        g[2].as_expr().subs({T: 0, r.as_expr(): 1}) * x for g, x in zip(columns, data)
    )
    assert gamma == -Rational(62, 15)
    H = Matrix([[x.as_expr() for x in row] for row in meridian(K.zero, one).hessian])
    e = Matrix([1, 1, 1, -1])
    penalty = factor((e.T * H.inv() * e)[0] / 2)
    assert cancel(penalty - (T + 1) * (T * T + 8 * T + 4) / (4 * T)) == 0
    print("Residual tilt component at t=0, u=0 =", gamma)
    print("Inverse-Hessian factor =", penalty)
    print("Limit of t * V^(2) at u=0 =", -(gamma**2))
    # Six amplitude pairs check the parameter independence of the tilt gradient.
    reference = meridian(K.zero, one)
    for lc, ld in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1)]:
        other = meridian(K.zero, one, lc, ld)
        assert other.tilt_gradient == reference.tilt_gradient
