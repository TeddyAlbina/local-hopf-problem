# Lean formalization

This directory is a Lean 4 / Mathlib formalization of the proof in [`../proof/PROOF.md`](../proof/PROOF.md). It is **not** a complete formalization from first principles. Mathlib does not yet have several things the proof needs:

- sectional curvature of Riemannian metrics on manifolds;
- the bundle of tangent 2-planes;
- Cheeger deformations;
- the $C^\infty$ topology on metrics.

The computer-certified identities of [`../verification/`](../verification) are also not re-run inside Lean.

So the formalization does two things:

1. **Proves** every step of the argument that can be stated with what Mathlib has: the logic of the proof, the analysis of section 3, the deduction of Proposition 5.2, the diagonal argument of section 7, and the closed-form facts about the profiles.
2. **States as explicit hypotheses** the geometric, literature and computational inputs. These are the fields of the structure `LocalHopfInput` in [`LocalHopf/Main.lean`](LocalHopf/Main.lean). Each field's docstring names its source.

The main theorem is

```lean
theorem local_hopf (I : LocalHopfInput) :
    ∃ x : ℕ → I.Met, (∀ n, I.PosSec (x n)) ∧ I.dist.TendstoCInf x I.gprod
```

There are no `axiom` declarations and no `sorry`. `#print axioms local_hopf` reports only Lean's standard `propext`, `Classical.choice` and `Quot.sound` (see [`PrintAxioms.lean`](PrintAxioms.lean)).

## Building

```sh
cd lean
lake exe cache get   # download prebuilt Mathlib
lake build
lake env lean PrintAxioms.lean
```

Pinned versions: Lean `v4.36.0`, Mathlib tag `v4.36.0`.

## What is proved

| File | PROOF.md | Proved in Lean |
|---|---|---|
| `Vectors.lean` | §1, §4 | $g_{\rm prod}(K_\omega,K_{\omega'})=\langle P_p\omega,\omega'\rangle$; $J_b$ is horizontal for $b\perp p_1,p_2$; $P_pb=2b$, hence $C_tK_b=K_b/(1+2t)$; the moment map $g_{\rm prod}(\xi,K_\omega)=\langle a(\xi),\omega\rangle$; the plane $\pi_{b,p}$ is spanned by $K_b,J_b$ and its spanning vectors are orthonormal. None of these needs $p_1\neq\pm p_2$. |
| `Profiles.lean` | §2, §5, §6, DERIVATIONS | The ODE for $P_0,Q$; $D\ge h>0$ and the positive denominators of $M_t,N_t,\varphi_t,\psi_t$; $\varphi_t>0$; the $d,d$ weight $1/(4h)=1/12$ at $t=1$; $\lambda_dM_t(0)\ne0$; the [BH] $c,d$ defect is $\ne0$ on $0<t<1$; the closed form of $\lvert A_0\rvert^2+\lvert B_0\rvert^2$ and the bound $<100$; $\mu_t=\lambda_c\lambda_d^2/(4h^2)=1/(40000h^2)>0$. |
| `Minimization.lean` | Lemma 3.1 | The lower bound $F(z,s)\ge s^2W_2+s^3W_3-C's^4$, proved from a coercivity estimate and the rescaled Taylor bound. This includes "minimizers are $O(s)$", the inequality $s^2\lambda a^2-Cs^3a\ge -C^2s^4/(4\lambda)$, and completing the square for $W_2$, $\bar z=-H^{-1}g$. |
| `Perturbation.lean` | Theorem 3.4 | Full compactness argument: from the far/near estimates and (H2), $u(\cdot,s)>0$ for $0<s<s_0$. |
| `QuadraticBound.lean` | Prop. 5.2 | From the normal form of Prop. 5.1 to $V^{(2)}\ge\rho\,\mathcal G$, including all Cauchy–Schwarz and Frobenius-norm estimates; then $V^{(2)}\ge\kappa_2\rho$ from $\mathcal G\ge\kappa_2$. |
| `TorusMean.lean` | §6 | For a trigonometric polynomial of degree $\le3$, both the four-rotation average and $\frac1{2\pi}\int_0^{2\pi}$ equal the constant term. |
| `Diagonal.lean` | §7 | The diagonal sequence $t_j\downarrow0$, $s_j\to0$ converges in every $C^j$. |
| `Main.lean` | §5–§7 | Extension of $U_2\ge h\kappa_2\rho$ from the dense set $\Omega$ to all of $Z$; $U_3=h\mu_t$ on all of $\Sigma$; (H2) with $\kappa=h\min(\kappa_2,\mu_t)$; positivity for small $s$; the main theorem. |
| `Consistency.lean` | — | A one-point model of `LocalHopfInput`. It shows that the hypotheses are not contradictory, so `local_hopf` is not vacuous. It does not show that the real geometry satisfies them. |

## What is assumed

These are the fields of `LocalHopfInput`. Formalizing them is the remaining work.

| Field | Source in PROOF.md | Kind |
|---|---|---|
| `g_tendsto`, `G_tendsto` | §1, §2, §7: $2g_t\to g_{\rm prod}$, and the tensors are smooth at fixed $t$ | geometry |
| `Z`, compact; `Ω_dense`, `Sig_dense` | §1, §3: $Z$ is a compact 4-manifold, and regular planes are dense in $Z$ and in $\Sigma$ | geometry |
| `frame` | §5: unit normal, $\rho=n_x^2+n_y^2=\lvert z'\rvert^2$, $\lvert x'\rvert\le1$ | geometry |
| `sign` | §3, "How it is applied": the sign of $u$ is the sign of ${\rm sec}$ | geometry |
| `U₂_cont`, `U₃_cont`, `tube` | Def. 3.2, Lemma 3.1 (Taylor step), Prop. 4.2 (Ziller's inequality, Lemma 4.1), Lemma 4.3 | analysis and literature |
| `U₂_normal_form` | Prop. 5.1 with Lemma 3.3 | `full_metric_check.py`, part A |
| `certificate` | Prop. 5.2: $\mathcal G\ge\frac3{10}$ on $(0,\frac1{100}]\times[-1,1]$ | `t0_certificate.py` (ball arithmetic) |
| `U₃_on_Sig` | Prop. 6.1, the $h^{(3)}$ correction, Lemma 3.3 | `full_metric_check.py`, part B, and §6 |

## Possible next steps

- **Lemma 4.1.** It is pure linear algebra in $\mathbb R^3$ and could be proved with the tools in `Vectors.lean`. That would reduce the `tube` hypothesis to Ziller's inequality and the Taylor step.
- **Taylor step.** Derive hypotheses (a) and (b) of `minimization_lower_bound` from $C^4$ bounds, using Mathlib's `taylor_mean_remainder`.
- **Certificates.** Export the exact coefficient functions in `verification/data/quadratic_coefficients.json` to Lean, and check $\mathcal G\ge\frac3{10}$ with verified interval arithmetic. The profile tails $f_{aa},f_{bb}$ need an analytic enclosure.
- **Geometry.** Model curvature of metrics on $S^2\times S^2$ once Mathlib's Riemannian geometry supports it.
