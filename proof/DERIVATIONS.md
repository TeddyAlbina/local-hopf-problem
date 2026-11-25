# How the tensors were found

This note explains where the profiles in section 2 of [PROOF.md](PROOF.md) come from, and why each one is necessary. The proof does not depend on it: every property used in the proof is certified directly by the scripts in [`verification/`](../verification).

The notation is that of PROOF.md. On the flat torus $\Gamma\times\Gamma$ of a regular pair, $U$ is the simultaneous rotation and $V$ the opposite rotation. The fields satisfy $J^\flat_\omega(U)=0$ and $J^\flat_\omega(V)=-\langle\omega,n\rangle$. $V^{(2)}$ is the quadratic coefficient of the curvature numerator after minimizing over plane tilts:
$$V^{(2)}=(\text{fixed-plane coefficient})-\tfrac12\,g^{T}H^{-1}g .$$
Here $g$ is the tilt gradient of the first variation, and $H$ is the tilt Hessian of the background.

## 1. Why small $t$ is hard

The tilt Hessian $H$ of $g_t$ splits into two 2×2 sectors, spanned by $\tfrac12(1,\sigma,1,-\sigma)$ and $\tfrac12(1,\sigma,-1,\sigma)$ for $\sigma=\pm1$. In each sector, write $A,B,D$ for the entries and $S=D-B^2/A$. Then, uniformly in $u$,
$$A\to2,\qquad B/A=t(1-\sigma u)+O(t^2),\qquad S/\big[t^3(1-\sigma u)\big]\to2.$$
The *soft* direction has stiffness of order $t^3$. Its penalty in $V^{(2)}$ is
$$\frac{(g_k-Bg_h/A)^2}{2S},$$
which is $O(t^{-3})$ unless the soft component of the tilt gradient vanishes to order $t^2$.

The first-order tensors of Brendle and Hung [BH] do not have this property. At $u=0$ the soft gradient tends to $8/5$, and $t^3V^{(2)}\to-16/25$ (`quadratic_coefficient.py`). Tuning a weight at a single separation does not help either. An affine weight that kills the gradient at $u=0$ still gives $t^3V^{(2)}\to-8649/3125$ at $u=4/5$ (`single_separation_weights.py`). The weights must depend on $u$.

## 2. First-order weights $P_t$, $Q$

Consider the $b$ tensor with general weights,
$$P(u)\,K^\flat_{p_2-p_1}\odot K^\flat_{e_y}+Q(u)\,K^\flat_{p_1+p_2}\odot J^\flat_{e_y},$$
on the torus whose normal is $e_y$. At $t=0$, its tilt gradient along the product-flat direction vanishes exactly when
$$2(1-u^2)P''+(3-5u)P'-P+\big((1-u)Q\big)'=0 .$$
For a polynomial $P$, integrating from $u=1$ determines $Q$. A search over low-degree $P$ for which the limiting coefficient $F_0$ is positive gives
$$P_0(u)=\tfrac{16}5+\tfrac u2-\tfrac{2u^2}5,\qquad Q(u)=2u^2-\tfrac{7u}{10}-4 .$$
The factor $1-t^2(1+u)^2$ in $P_t$ removes the remaining $O(t)$ part of the soft gradient. The $a$ tensor uses the same weights at $-u$, which is its image under the antipodal map on the first factor.

The soft penalty is then $O(t)$, so
$$F_0(u)=\tfrac1{400}\big(1152u^6-400u^5-10400u^4-2956u^3+18544u^2+6920u-4643\big)-18\,\partial_\delta^2 f_{bb}(\delta).$$
Here $\delta=\arccos u$ is the longitude separation. Replace $\sqrt{15}$ by the rational lower bound $38729833462\cdot10^{-10}$, which is valid since its coefficient is positive. Then $(4-u)(F_0-1)$ becomes a rational polynomial with no roots in $[-1,1]$. So $F_0>1$ (`first_order_weights.py`).

## 3. The $d,d$ weight $1/(4h)$

On the equatorial torus the quadratic coefficient must vanish identically. The $d,d$ contribution there is $(1-t)\cos^2\delta/(6D)$ with the [BH] weight $\tfrac1{12}$, and exactly $0$ with the weight $\tfrac1{4h}$ (`quadratic_coefficient.py`). At $t=1$ the two weights agree.

## 4. Normal corrections $M_t$, $N_t$ in $h_{ad}$, $h_{bd}$

**Which terms survive.** Split each fixed axis into its components normal to and inside the plane of $\Gamma$.
- Reflection in that plane fixes $\Gamma\times\Gamma$.
- Normal–normal and inside–inside components are even under the reflection. They pull back to constants on the torus and contribute nothing to the quadratic coefficient.
- Products of even and odd components vanish by parity.

The odd part of the $a$ tensor is $n_x$ times the axial $a$ variation. The odd parts of the $d$ tensor are
$$n_z\,D_1(x')+n_x\,D_2(z'),\qquad D_1(w)=K^\flat_w\odot J^\flat_n,\quad D_2(w)=K^\flat_n\odot J^\flat_w .$$
With the [BH] $h_{ad}$, the $a,d$ coefficient is
$$M_t(u)\,w_xn_xn_z+N_t(u)\,w_zn_x^2,\qquad w=p_2-p_1 .$$
This is the raw curvature jet, minus the pairing of tilt gradients through $H^{-1}$, plus the linear term of $h_{ad}$.

**The correction.** Under $U$, each coordinate of $w$ is a first harmonic, and $u$ is constant. The tensor
$$-M_t(u)\,w_x\,J^\flat_{e_x}\odot J^\flat_{e_z}-N_t(u)\,w_z\,J^\flat_{e_x}\odot J^\flat_{e_x}$$
pulls back only to the $VV$ entry. Its linearized curvature $-U^2h(V,V)/2$, with the factor 2 in front of $h_{ad}$, cancels both terms exactly. The same holds for $b,d$ with $v=p_1+p_2$ and $M_t(-u)$, $N_t(-u)$. Since the formula has no division by $|p_1\times p_2|$, the correction is smooth on $M$. The coefficient tables are in `normal_corrections.py`.

**Necessity.** Correct only $h_{ad}$, and tilt the normal of the pair $(e_x,e_y)$ toward $e_y$. The part of the slope that is odd under simultaneous antipodes is
$$\lambda_dM_t(0)=-\lambda_d\frac{106t^3+37t^2+18t+67}{40(1+t)(1+2t)^2}\neq0 .$$
The quadratic coefficient vanishes on the equator, so this pair or its antipode has negative quadratic coefficient on one side. With neither correction, the $a,d$ slope is the same expression with the opposite sign (`first_order_weights.py`).

## 5. The mixed tensor $h_{cd}$

**The uncorrected term.** Let $C(w)=K^\flat_n\odot K^\flat_w$ and $\Delta=p_2p_2^T-p_1p_1^T$. In the frame $p_1=(1,0)$, $p_2=(\cos\delta,\sin\delta)$, the bilinear form $B$ of the minimized quadratic coefficient satisfies, entrywise,
$$4B\big(C(\cdot),D_1(\cdot)\big)=\varphi_t\,\Delta,\qquad 4\,{\rm sym}\,B\big(C(\cdot),D_2(\cdot)\big)=\psi_t\,\Delta .$$
Without an $h_{cd}$ term, the $c,d$ coefficient is therefore
$$n_z^2\varphi_t\Delta_{zx}+n_xn_z\psi_t\Delta_{zz}.$$

**The correction.** Under $U$ the traceless matrix $\Delta$ rotates at twice the angular speed, so $U^2\Delta=-4\Delta$. The tensor
$$h_{cd}=-\tfrac14\big[\varphi_t\Delta_{zx}J^\flat_{e_z}\odot J^\flat_{e_z}+\psi_t\Delta_{zz}J^\flat_{e_x}\odot J^\flat_{e_z}\big]$$
pulls back only to the $VV$ entry. With the factor 2 in the metric, its linear curvature cancels the term above at every regular point and every $t>0$ (`cd_tensor.py`). Since $D=(1+t)^2-t^2u^2\ge1+2t$, it is smooth.

**Necessity.** The [BH] tensor $h_{cd}$ contains $G(\sin2\theta/\sqrt3)$ with $G(x)=\arctan x/x$. At $(e_x,e_y)$, tilt $n$ toward $e_x$. The $c,d$ coefficient of the normal derivative is then
$$\frac{(1-t)(3t+1)}{4(1+t)^2(1+2t)^2},$$
which is nonzero for $0<t<1$. So that tensor gives a negative quadratic coefficient on one side of the equator (`normal_corrections.py`).

## 6. The remaining channels

- **$c,c$:** vanishes with $h_{cc}$, for all three entries of the in-plane matrix.
- **Linear in $\lambda_c$** (from $h_{ac}$, $h_{bc}$ and the cross terms): of the form $n_z[n_xA_t+n_yB_t]\cdot z'$. The limits at $t=0$ are
  $$A_0=f(u)(p_2+p_1),\qquad B_0=f(-u)(p_2-p_1),\qquad f(v)=\tfrac1{20}(-168v^3+154v^2+148v-99).$$
- **$d,d$:** after cancelling against $h_{dd}$, it has the form $2n_zn_x\,x'^TM^1_tz'+n_x^2\,z'^TM^2_tz'$. The limits at $t=0$ have entries of absolute value at most $\tfrac14$.

These forms, in general position and for every $t>0$, are Proposition 5.1 of PROOF.md. The bounds used are Proposition 5.2.

## Reference

**[BH]** S. Brendle, P.-K. Hung, *A metric on $S^2\times S^2$ with positive sectional curvature*, arXiv:2608.19068.
