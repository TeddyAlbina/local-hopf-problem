# The Local Hopf Problem on $S^2\times S^2$

> **Local Hopf Problem I** (Bettiol, 2015). Is there a sequence $g_n$ of metrics on $S^2\times S^2$ with ${\rm sec}_{g_n}>0$ that converges to the standard product metric?

This repository contains a computer-assisted proof that the answer is yes.

> **Theorem.** There is a sequence of smooth Riemannian metrics on $S^2\times S^2$ with strictly positive sectional curvature that converges to the product of round metrics in $C^\infty$.

The product metric $g_{\rm prod}$ has flat planes: take one tangent direction from each sphere. No first-order change of $g_{\rm prod}$ makes all of them positive (Berger, Strake). Bourguignon, Deschamps and Sentenac showed that analytic deformations cannot succeed for products without infinitesimal isometries. The round product has a six-dimensional isometry group, so their result does not apply to it. Bettiol posed the local question and recorded it as open, alongside the Hopf Problem itself, which Brendle and Hung have since answered with a metric far from the product. The theorem shows that arbitrarily small smooth changes of $g_{\rm prod}$ make every plane positive.

The full argument is in [**proof/PROOF.md**](proof/PROOF.md); section 9 there compares it with the known obstructions. The same proof as a LaTeX article is in [paper/](paper) ([PDF](paper/local-hopf.pdf)). [proof/DERIVATIONS.md](proof/DERIVATIONS.md) explains how the tensors were found and why each modification is needed, and [verification/](verification) certifies every algebraic step.

## How the proof works

1. **Background.** For a Cheeger parameter $t>0$, deform $g_{\rm prod}$ by the diagonal $SO(3)$-action. This gives a nonnegatively curved metric $g_t$, whose flat planes form a compact 4-manifold $Z$. As $t\to0$, $2g_t\to g_{\rm prod}$.
2. **Perturbation.** Perturb to $G_{t,s}=g_t+s\,h^{(1)}+s^2h^{(2)}+s^3h^{(3)}$. This follows Brendle and Hung, who used it at $t=1$ ([arXiv:2608.19068](https://arxiv.org/abs/2608.19068)). Their tensors fail as $t\to0$, because the tilt stiffness of the flat planes degenerates like $t^3$. We use first-order weights, normal corrections and a mixed tensor that work for all small $t$.
3. **Reduction.** A perturbation theorem on the bundle of 2-planes reduces positivity for small $s$ to three conditions:
   - quadratic growth of the background curvature away from $Z$, from Ziller's curvature inequality for Cheeger deformations, for every $t>0$;
   - a lower bound $V^{(2)}\ge\frac3{10}\rho$ on the second-order coefficient for $0<t\le\frac1{100}$, where $\rho$ vanishes only on an equatorial torus;
   - a positive third-order coefficient on that torus, with exact mean $1/(40000(1+2t)^2)$, made constant by a conformal $h^{(3)}$.
4. **Limit.** Each fixed $t\in(0,\frac1{100}]$ gives positive curvature for small $s>0$. Letting $t\to0$ proves the theorem.

## Verification

- **Exact, in general position.** The structure of the quadratic coefficient is checked in general position for every $t$ and every separation. The check is exact over $\mathbb Q(t,r)$, at 165 rotations on which polynomials of degree $\le4$ on $SO(3)$ are determined.
- **Rigorous bound.** The bound $V^{(2)}\ge\frac3{10}\rho$ for $t\le\frac1{100}$ is proved by ball arithmetic over the whole parameter box. The box includes the coincident and antipodal pairs.
- **Exact cubic.** The cubic mean on the torus is computed exactly from the complete tensors.
- **Two backends.** All scripts run on SymPy or on FLINT (`HOPF_BACKEND=flint`), which implement the same exact field. Logs, produced with FLINT, are in [verification/logs](verification/logs).

```sh
python3 -m venv .venv && . .venv/bin/activate
pip install -r requirements.txt
cd verification
HOPF_BACKEND=flint python full_metric_check.py
python t0_certificate.py
```

See [verification/README.md](verification/README.md) for what each script certifies.

## References

- R. G. Bettiol, *On different notions of positivity of curvature*, PhD thesis, University of Notre Dame (2015).
- J.-P. Bourguignon, A. Deschamps, P. Sentenac, *Conjecture de H. Hopf sur les produits de variétés*, Ann. Sci. ÉNS (4) 5 (1972); *Quelques variations particulières d'un produit de métriques*, Ann. Sci. ÉNS (4) 6 (1973).
- M. Strake, *Curvature increasing metric variations*, Math. Ann. 276 (1987).
- S. Brendle, P.-K. Hung, *A metric on $S^2\times S^2$ with positive sectional curvature*, [arXiv:2608.19068](https://arxiv.org/abs/2608.19068).
- W. Ziller, *Examples of Riemannian manifolds with non-negative sectional curvature*, Surveys in Differential Geometry XI (2007).
- M. Müter, *Krümmungserhöhende Deformationen mittels Gruppenaktionen*, Dissertation, Münster (1987).

## License

The verification code is released under the MIT License, and the proof, paper and documentation under CC BY 4.0. See [LICENSE](LICENSE).
