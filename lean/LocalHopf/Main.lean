/-
Copyright (c) 2026. Released under the MIT License (see ../../LICENSE).

# The main theorem, assembled from its inputs

Mathlib does not yet have sectional curvature of Riemannian metrics on manifolds, the bundle of
tangent 2-planes, Cheeger deformations, or the `C^∞` topology on metrics. The geometric objects
are therefore taken as *data* in the structure `LocalHopfInput`, and every statement that PROOF.md
obtains from differential geometry, from the literature, or from a verification script is a
*field* of that structure, with its source in the docstring. Nothing is postulated with `axiom`,
and no proof uses `sorry`.

`local_hopf` then derives, with every step checked by Lean:

1. `U₂ ≥ h·(3/10)·ρ` on all of `Z`
   (Prop. 5.2 deduction `V2_ge_kappa_rho`, then density of the regular planes);
2. `U₃ = h·μ_t` on all of `Σ = {ρ = 0}` (density in `Σ`), with `μ_t = 1/(40000 h²) > 0`;
3. hypothesis (H2) of Theorem 3.4 with `κ = h·min(3/10, μ_t)`;
4. positivity of `G_{t,s}` for `0 < s < s₀(t)` (Theorem 3.4, `eventually_positive`);
5. a sequence of positively curved metrics converging to `g_prod` in `C^∞`
   (section 7, `diagonal_sequence`).
-/
import LocalHopf.Profiles
import LocalHopf.Perturbation
import LocalHopf.QuadraticBound
import LocalHopf.Diagonal

namespace LocalHopf

open Set Filter Topology Profiles Matrix

/-- The threshold `t₀ = 1/100` of Proposition 5.2. -/
noncomputable def t₀ : ℝ := 1 / 100

/-- The inputs of the proof. Each field is labelled with where PROOF.md gets it:

* **[geometry]** standard differential geometry not yet in Mathlib;
* **[Ziller]** the curvature inequality of Ziller / Müter for Cheeger deformations (§4);
* **[analysis]** Taylor expansion and the tubular-neighbourhood estimates of §3–4;
* **[script: …]** a computation certified by the named script in `verification/`. -/
structure LocalHopfInput where
  /-- Smooth Riemannian metrics on `S² × S²`. -/
  Met : Type
  /-- The `C^j` distances, defining `C^∞` convergence. -/
  dist : CkDistances Met
  /-- `PosSec g` : `g` has strictly positive sectional curvature. -/
  PosSec : Met → Prop
  /-- The product of unit round metrics. -/
  gprod : Met
  /-- `g t = 2 g_t`, twice the background (scaling does not change the sign of `sec`). -/
  g : ℝ → Met
  /-- `G t s = 2 G_{t,s}`. -/
  G : ℝ → ℝ → Met
  /-- **[geometry]** `2 g_t → g_prod` in `C^∞` as `t → 0⁺` (§1). -/
  g_tendsto : ∀ j, Tendsto (fun t => dist.d j (g t) gprod) (𝓝[>] 0) (𝓝 0)
  /-- **[geometry]** At fixed `t` the tensors `h⁽ⁱ⁾` are smooth, so `G_{t,s} → g_t` in each `C^j`
  as `s → 0⁺` (§2, §7). -/
  G_tendsto : ∀ t ∈ Ioc 0 t₀, ∀ j,
    Tendsto (fun s => dist.d j (G t s) (g t)) (𝓝[>] 0) (𝓝 0)
  /-- The bundle `𝒩` of tangent 2-planes. -/
  𝒩 : Type
  /-- The flat locus `Z ⊂ 𝒩`, a compact embedded 4-manifold (§1). -/
  Z : Type
  [topZ : TopologicalSpace Z]
  [compactZ : CompactSpace Z]
  /-- `ρ = 1 - ⟨n, e_z⟩²` on `Z`. -/
  ρ : Z → ℝ
  ρ_cont : Continuous ρ
  ρ_nonneg : ∀ x, 0 ≤ ρ x
  /-- The planes over regular pairs `p₁ ≠ ± p₂`. -/
  Ω : Set Z
  /-- **[geometry]** `Ω` is dense in `Z` (§3). -/
  Ω_dense : Dense Ω
  /-- **[geometry]** `Σ ∩ Ω` is dense in `Σ = {ρ = 0}` (regular pairs are dense in the equatorial
  torus). -/
  Sig_dense : ∀ x, ρ x = 0 → x ∈ closure {y | y ∈ Ω ∧ ρ y = 0}
  /-- Frame data at a regular plane: the components `n_x, n_y, n_z` of the unit normal, the
  projections `x', z'` of `e_x, e_z` to the plane of `Γ` (in an orthonormal basis of that plane),
  and `u = ⟨p₁, p₂⟩`. -/
  nx : Z → ℝ
  ny : Z → ℝ
  nz : Z → ℝ
  xp : Z → Fin 2 → ℝ
  zp : Z → Fin 2 → ℝ
  sep : Z → ℝ
  /-- **[geometry]** Elementary facts about the frame (§5): `n` is a unit vector,
  `ρ = n_x² + n_y² = |z'|²`, `|x'| ≤ 1`, and `u ∈ [-1, 1]`. -/
  frame : ∀ x ∈ Ω, nx x ^ 2 + ny x ^ 2 = ρ x ∧ nx x ^ 2 + ny x ^ 2 + nz x ^ 2 = 1 ∧
    zp x ⬝ᵥ zp x = ρ x ∧ xp x ⬝ᵥ xp x ≤ 1 ∧ sep x ∈ Icc (-1 : ℝ) 1
  /-- The coefficient functions `F^x_t, F^y_t, A_t, B_t, M¹_t, M²_t` of Proposition 5.1,
  functions of `(t, u)`. -/
  Fx : ℝ → ℝ → ℝ
  Fy : ℝ → ℝ → ℝ
  A : ℝ → ℝ → Fin 2 → ℝ
  B : ℝ → ℝ → Fin 2 → ℝ
  M1 : ℝ → ℝ → Matrix (Fin 2) (Fin 2) ℝ
  M2 : ℝ → ℝ → Matrix (Fin 2) (Fin 2) ℝ
  /-- `u t π s = R_{G_{t,s}}(e₁, e₂, e₁, e₂)` for a `g_t`-orthonormal basis of `π`. -/
  u : ℝ → 𝒩 → ℝ → ℝ
  /-- The coefficients `U₂, U₃` of Definition 3.2. -/
  U₂ : ℝ → Z → ℝ
  U₃ : ℝ → Z → ℝ
  /-- **[analysis]** `U₂, U₃` are continuous on `Z` (Definition 3.2). -/
  U₂_cont : ∀ t ∈ Ioc 0 t₀, Continuous (U₂ t)
  U₃_cont : ∀ t ∈ Ioc 0 t₀, Continuous (U₃ t)
  /-- **[geometry]** The sign of `u` is the sign of `sec_{G_{t,s}}` (§3, "How it is applied"),
  and `sec_{2G} = sec_G / 2`. -/
  sign : ∀ t ∈ Ioc 0 t₀, ∀ s, 0 < s → (∀ π, 0 < u t π s) → PosSec (G t s)
  /-- **[Ziller] + [analysis]** The estimates used in the proof of Theorem 3.4:
  away from `Z`, `u ≥ σδ² - Cs` (Prop. 4.2, from Ziller's inequality and Lemma 4.1);
  in the tube, `u ≥ s²U₂ + s³U₃ - Cs⁴` (Lemma 3.1 with (H0) and Lemma 4.3 (H1), for the
  `C⁴` expansion of `u`; see `minimization_lower_bound` for the step from the Taylor bounds). -/
  tube : ∀ t ∈ Ioc 0 t₀, ∃ (far : 𝒩 → Prop) (proj : 𝒩 → Z) (a C s₁ : ℝ), 0 < a ∧ 0 < s₁ ∧
    (∀ π s, 0 < s → s ≤ s₁ → far π → a - C * s ≤ u t π s) ∧
    (∀ π s, 0 < s → s ≤ s₁ → ¬ far π →
      s ^ 2 * U₂ t (proj π) + s ^ 3 * U₃ t (proj π) - C * s ^ 4 ≤ u t π s)
  /-- **[script: `full_metric_check.py`, part A] + [analysis: Lemma 3.3]**
  Proposition 5.1 in the graph chart, transported to `U₂ = h V⁽²⁾` by Lemma 3.3. -/
  U₂_normal_form : ∀ t ∈ Ioc 0 t₀, ∀ x ∈ Ω,
    U₂ t x = hh t * V2form lamC lamD (Fx t (sep x)) (Fy t (sep x)) (nx x) (ny x) (nz x)
      (A t (sep x)) (B t (sep x)) (xp x) (zp x) (M1 t (sep x)) (M2 t (sep x))
  /-- **[script: `t0_certificate.py`]** The ball-arithmetic bound `𝒢 ≥ 3/10` on
  `(0, 1/100] × [-1, 1]` (Proposition 5.2). -/
  certificate : ∀ t ∈ Ioc 0 t₀, ∀ v ∈ Icc (-1 : ℝ) 1,
    3 / 10 ≤ Gcal lamC lamD (Fx t v) (Fy t v) (A t v) (B t v) (M1 t v) (M2 t v)
  /-- **[script: `full_metric_check.py`, part B] + [analysis: §6, Lemma 3.3]**
  After the conformal correction `h⁽³⁾ = 2hχ_t g_t`, the cubic coefficient on `Σ ∩ Ω` equals
  `h μ_t`, where `μ_t = λ_c λ_d² / (4h²)` is the exact rotation mean of Proposition 6.1. -/
  U₃_on_Sig : ∀ t ∈ Ioc 0 t₀, ∀ x ∈ Ω, ρ x = 0 → U₃ t x = hh t * mu t

attribute [instance] LocalHopfInput.topZ LocalHopfInput.compactZ

namespace LocalHopfInput

variable (I : LocalHopfInput)

theorem t_nonneg {t : ℝ} (ht : t ∈ Ioc 0 t₀) : 0 ≤ t := ht.1.le

/-- **Proposition 5.2 on all of `Z`.** `U₂ ≥ h κ₂ ρ` with `κ₂ = 3/10`. -/
theorem U₂_lower {t : ℝ} (ht : t ∈ Ioc 0 t₀) (x : I.Z) :
    hh t * (3 / 10) * I.ρ x ≤ I.U₂ t x := by
  have hS : IsClosed {x | hh t * (3 / 10) * I.ρ x ≤ I.U₂ t x} :=
    isClosed_le (continuous_const.mul I.ρ_cont) (I.U₂_cont t ht)
  have hΩ : I.Ω ⊆ {x | hh t * (3 / 10) * I.ρ x ≤ I.U₂ t x} := by
    intro y hy
    obtain ⟨h1, h2, h3, h4, h5⟩ := I.frame y hy
    have hb := V2_ge_kappa_rho lamC lamD _ _ (I.nx y) (I.ny y) (I.nz y) (I.ρ y) (3 / 10)
      _ _ (I.xp y) (I.zp y) _ _ (by norm_num [lamC]) h1 h2 h3 h4 (I.certificate t ht _ h5)
    show hh t * (3 / 10) * I.ρ y ≤ I.U₂ t y
    rw [I.U₂_normal_form t ht y hy, mul_assoc]
    exact mul_le_mul_of_nonneg_left hb (hh_pos (t_nonneg ht)).le
  have := hS.closure_subset_iff.2 hΩ
  rw [I.Ω_dense.closure_eq] at this
  exact this (mem_univ x)

/-- `U₃ = h μ_t` on all of `Σ`. -/
theorem U₃_on_Sig_all {t : ℝ} (ht : t ∈ Ioc 0 t₀) (x : I.Z) (hx : I.ρ x = 0) :
    I.U₃ t x = hh t * mu t := by
  have hS : IsClosed {x | I.U₃ t x = hh t * mu t} :=
    isClosed_eq (I.U₃_cont t ht) continuous_const
  have hsub : {y | y ∈ I.Ω ∧ I.ρ y = 0} ⊆ {x | I.U₃ t x = hh t * mu t} :=
    fun y hy => I.U₃_on_Sig t ht y hy.1 hy.2
  exact hS.closure_subset_iff.2 hsub (I.Sig_dense x hx)

/-- **Hypothesis (H2)** with `κ = h·min(3/10, μ_t) > 0`. -/
theorem H2 {t : ℝ} (ht : t ∈ Ioc 0 t₀) :
    let κ := hh t * min (3 / 10) (mu t)
    0 < κ ∧ (∀ x, κ * I.ρ x ≤ I.U₂ t x) ∧ (∀ x, I.ρ x = 0 → κ ≤ I.U₃ t x) := by
  have hh0 := hh_pos (t_nonneg ht)
  have hμ := mu_pos (t_nonneg ht)
  refine ⟨mul_pos hh0 (lt_min (by norm_num) hμ), fun x => ?_, fun x hx => ?_⟩
  · calc hh t * min (3 / 10) (mu t) * I.ρ x ≤ hh t * (3 / 10) * I.ρ x := by
          apply mul_le_mul_of_nonneg_right _ (I.ρ_nonneg x)
          exact mul_le_mul_of_nonneg_left (min_le_left _ _) hh0.le
      _ ≤ I.U₂ t x := I.U₂_lower ht x
  · rw [I.U₃_on_Sig_all ht x hx]
    exact mul_le_mul_of_nonneg_left (min_le_right _ _) hh0.le

/-- **Positivity at fixed `t`** (section 7, first paragraph): Theorem 3.4 gives `s₀(t) > 0`. -/
theorem pos_small_s {t : ℝ} (ht : t ∈ Ioc 0 t₀) :
    ∃ s₀ > 0, ∀ s, 0 < s → s < s₀ → I.PosSec (I.G t s) := by
  obtain ⟨hκ, hU₂, hU₃⟩ := I.H2 ht
  obtain ⟨far, proj, a, C, s₁, ha, hs₁, hfar, hnear⟩ := I.tube t ht
  obtain ⟨s₀, hs₀, hpos⟩ := eventually_positive (I.u t) (I.U₂ t) (I.U₃ t) I.ρ (I.U₃_cont t ht)
    I.ρ_cont _ hκ I.ρ_nonneg hU₂ hU₃ far proj a C s₁ ha hs₁ hfar hnear
  exact ⟨s₀, hs₀, fun s hs hss => I.sign t ht s hs fun π => hpos π s hs hss⟩

end LocalHopfInput

/-- **Main theorem** (Local Hopf Problem I). Given the inputs listed in `LocalHopfInput`, there is
a sequence of metrics on `S² × S²` with strictly positive sectional curvature converging to the
product of round metrics in `C^∞`. -/
theorem local_hopf (I : LocalHopfInput) :
    ∃ x : ℕ → I.Met, (∀ n, I.PosSec (x n)) ∧ I.dist.TendstoCInf x I.gprod :=
  diagonal_sequence I.dist I.PosSec I.G I.g I.gprod t₀ (by norm_num [t₀])
    (fun t ht => I.pos_small_s ht) I.G_tendsto I.g_tendsto

end LocalHopf
