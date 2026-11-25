"""Affine b weights tuned at one separation (DERIVATIONS.md, section 1).

Choose the b weight (A + B u)(2 + u) so that the tilt gradient vanishes at
u = p1.p2 = 0. This repairs the quadratic coefficient at that separation
(positive limit) but not at others: at u = 4/5 the product-flat tilt
derivative tends to 186/125 and t^3 V^(2) -> -8649/3125. So the weight must
depend on u as in PROOF.md section 2. Curvature is computed from metric jets.
"""

from quadratic_coefficient import K, meridian, one, t
from sympy import Rational, chebyshevt, factor, simplify, sqrt


def solve_for_weights():
    """The meridian tilt gradient is affine in the new intercept and slope."""

    def gradient(a, b):
        return meridian(K.zero, one, b_intercept=a, b_slope=b).tilt_gradient

    original = gradient(6 * one / 5, 0)
    assert (original[0] - original[2]).as_expr().subs(t.as_expr(), 0) == Rational(8, 5)
    g = gradient(0, 0)
    a = [x - y for x, y in zip(gradient(1, 0), g)]
    b = [x - y for x, y in zip(gradient(0, 1), g)]
    determinant = a[0] * b[2] - a[2] * b[0]
    intercept = (b[0] * g[2] - b[2] * g[0]) / determinant
    slope = (a[2] * g[0] - a[0] * g[2]) / determinant
    assert gradient(intercept, slope) == [K.zero] * 4
    return intercept, slope


def fbb_second_derivative(cosine):
    # delta is longitude separation. The Fourier tail, differentiated
    # twice, is a geometric series with real sum (k*cos(delta)-k^2)/D.
    harmonics = (
        Rational(2203841, 1920),
        Rational(5579141, 153600),
        Rational(214597, 103680),
        Rational(48211, 307200),
        Rational(251, 19200),
        Rational(25, 165888),
    )
    k = 4 - sqrt(15)
    tail = (k * cosine - k**2) / (1 - 2 * k * cosine + k**2)
    return simplify(
        -sum(n * n * h * chebyshevt(n, cosine) for n, h in enumerate(harmonics, 1))
        + Rational(42025, 6) * sqrt(Rational(5, 3)) * tail
    )


def full_quadratic(cosine, sine, intercept, slope):
    def result(lc, ld):
        return meridian(cosine, sine, lc, ld, b_intercept=intercept, b_slope=slope)

    sample = result(0, 0)
    # Quadratic dependence on the two free amplitudes. Six unisolvent
    # evaluations over Q(t) suffice. Only h_bb^(2) restricts nontrivially
    # to this great-circle product, since its normal is global e_y.
    for lc, ld in [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1)]:
        other = result(lc, ld)
        assert other.value == sample.value
        assert other.tilt_gradient == sample.tilt_gradient
    ts = t.as_expr()
    fpp = fbb_second_derivative(cosine.as_expr())
    return sample.tilt_gradient, sample.value.as_expr() - 18 * fpp / (1 + 2 * ts) ** 2


if __name__ == "__main__":
    A, B = solve_for_weights()
    ts = t.as_expr()
    print("A =", factor(A.as_expr()))
    print("B =", factor(B.as_expr()))
    for u, v in [(K.zero, one), (4 * one / 5, 3 * one / 5)]:
        grad, V = full_quadratic(u, v, A, B)
        tilt_limit = (grad[0] - grad[2]).as_expr().limit(ts, 0)
        print("p1 dot p2 =", u)
        print("  product-flat tilt derivative at t -> 0:", tilt_limit)
        if not u:
            limit = simplify(V.limit(ts, 0))
            assert grad == [K.zero] * 4
            assert limit.is_positive
            print("  limit of full V:", limit, "~", limit.evalf(12))
        else:
            leading = simplify((ts**3 * V).limit(ts, 0))
            assert leading == -Rational(8649, 3125)
            assert tilt_limit == Rational(186, 125)
            print("  limit of t^3 V:", leading)
        print("  V at t=1/10:", V.subs(ts, Rational(1, 10)).evalf(12))
