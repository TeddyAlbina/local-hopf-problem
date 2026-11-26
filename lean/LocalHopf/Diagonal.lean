/-
Copyright (c) 2026. Released under the MIT License (see ../../LICENSE).

# Section 7: the diagonal argument

For each `t ∈ (0, t₀]` there is `s₀(t) > 0` such that `G_{t,s}` has positive curvature for
`0 < s < s₀(t)`; at fixed `t` the tensors are smooth, so `G_{t,s} → g_t` in every `C^j` as `s → 0`;
and `g_t → g_prod` in `C^∞` as `t → 0`. PROOF.md, section 7, chooses `t_j ↓ 0` and `s_j` to get a
sequence converging to `g_prod` in `C^∞`.

`C^∞` convergence is modelled by a family of distances `d_j` (the `C^j` distances), increasing in
`j`, with the triangle inequality: `x_n → g` in `C^∞` iff `d_j(x_n, g) → 0` for every `j`.
-/
import Mathlib.Order.Filter.AtTopBot.Basic
import Mathlib.Topology.Order.Basic
import Mathlib.Topology.Algebra.Order.Field
import Mathlib.Analysis.SpecificLimits.Basic
import Mathlib.Tactic

namespace LocalHopf

open Filter Topology Set

/-- The `C^j` distances on a space of metrics. -/
structure CkDistances (Met : Type*) where
  /-- `d j a b` is the `C^j` distance. -/
  d : ℕ → Met → Met → ℝ
  nonneg : ∀ j a b, 0 ≤ d j a b
  triangle : ∀ j a b c, d j a c ≤ d j a b + d j b c
  mono : ∀ a b, Monotone fun j => d j a b

/-- `C^∞` convergence: convergence in every `C^j`. -/
def CkDistances.TendstoCInf {Met : Type*} (D : CkDistances Met) (x : ℕ → Met) (g : Met) : Prop :=
  ∀ j, Tendsto (fun n => D.d j (x n) g) atTop (𝓝 0)

/-- **Section 7.** A two-parameter family with positivity for small `s` at each `t`, converging
to `g_t` as `s → 0` and with `g_t → g_prod` as `t → 0`, yields a sequence of positively curved
metrics converging to `g_prod` in `C^∞`. -/
theorem diagonal_sequence {Met : Type*} (D : CkDistances Met) (Pos : Met → Prop)
    (G : ℝ → ℝ → Met) (g : ℝ → Met) (gprod : Met) (t₀ : ℝ) (ht₀ : 0 < t₀)
    (hpos : ∀ t ∈ Ioc 0 t₀, ∃ s₀ > 0, ∀ s, 0 < s → s < s₀ → Pos (G t s))
    (hG : ∀ t ∈ Ioc 0 t₀, ∀ j, Tendsto (fun s => D.d j (G t s) (g t)) (𝓝[>] 0) (𝓝 0))
    (hg : ∀ j, Tendsto (fun t => D.d j (g t) gprod) (𝓝[>] 0) (𝓝 0)) :
    ∃ x : ℕ → Met, (∀ n, Pos (x n)) ∧ D.TendstoCInf x gprod := by
  -- `t_n = t₀ / (n+1)`
  set tn : ℕ → ℝ := fun n => t₀ / (n + 1)
  have htn : ∀ n, tn n ∈ Ioc 0 t₀ := fun n => by
    refine ⟨by positivity, ?_⟩
    simp only [tn]
    rw [div_le_iff₀ (by positivity)]
    nlinarith [(Nat.cast_nonneg n : (0:ℝ) ≤ n)]
  -- choose `s_n` with positivity and `d_n(G, g) < 1/(n+1)`
  have hchoose : ∀ n, ∃ s, Pos (G (tn n) s) ∧ D.d n (G (tn n) s) (g (tn n)) < 1 / (n + 1) := by
    intro n
    obtain ⟨s₀, hs₀, hP⟩ := hpos (tn n) (htn n)
    have e1 : ∀ᶠ s in 𝓝[>] (0:ℝ), D.d n (G (tn n) s) (g (tn n)) < 1 / (n + 1) :=
      (hG (tn n) (htn n) n).eventually (gt_mem_nhds (by positivity))
    have e2 : ∀ᶠ s in 𝓝[>] (0:ℝ), 0 < s ∧ s < s₀ := by
      have h1 : ∀ᶠ s in 𝓝[>] (0:ℝ), 0 < s := self_mem_nhdsWithin
      have h2 : ∀ᶠ s in 𝓝[>] (0:ℝ), s < s₀ :=
        nhdsWithin_le_nhds (gt_mem_nhds hs₀)
      exact h1.and h2
    obtain ⟨s, hs1, hs2⟩ := (e1.and e2).exists
    exact ⟨s, hP s hs2.1 hs2.2, hs1⟩
  choose sn hsnP hsnd using hchoose
  refine ⟨fun n => G (tn n) (sn n), hsnP, fun j => ?_⟩
  -- `t_n → 0⁺`
  have htn0 : Tendsto tn atTop (𝓝[>] 0) := by
    apply tendsto_nhdsWithin_iff.2
    refine ⟨?_, Eventually.of_forall fun n => (htn n).1⟩
    have := (tendsto_one_div_add_atTop_nhds_zero_nat).const_mul t₀
    rw [mul_zero] at this
    exact Tendsto.congr (fun n => by simp only [tn]; ring) this
  have hgn := (hg j).comp htn0
  -- squeeze: `d_j(x_n, g_prod) ≤ 1/(n+1) + d_j(g_{t_n}, g_prod)` for `n ≥ j`
  have hbound : ∀ᶠ n in atTop, D.d j (G (tn n) (sn n)) gprod ≤
      1 / ((n : ℝ) + 1) + D.d j (g (tn n)) gprod := by
    filter_upwards [eventually_ge_atTop j] with n hn
    have h1 := D.triangle j (G (tn n) (sn n)) (g (tn n)) gprod
    have h2 : D.d j (G (tn n) (sn n)) (g (tn n)) ≤ D.d n (G (tn n) (sn n)) (g (tn n)) :=
      D.mono _ _ hn
    have h3 := hsnd n
    linarith
  have hlim : Tendsto (fun n : ℕ => 1 / ((n : ℝ) + 1) + D.d j (g (tn n)) gprod) atTop (𝓝 0) := by
    simpa using tendsto_one_div_add_atTop_nhds_zero_nat.add hgn
  exact squeeze_zero' (Eventually.of_forall fun n => D.nonneg _ _ _) hbound hlim

end LocalHopf
