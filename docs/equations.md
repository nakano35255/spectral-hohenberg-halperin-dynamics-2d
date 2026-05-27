# 基礎方程式

このページでは、本ソルバーが時間発展を計算する基礎方程式をまとめます。


## 1. 基本変数と空間離散化

### 1.1 基本場

本ソルバーは、以下の場を基本変数として扱います。

- $\rho(\boldsymbol{x},t)$
- $\boldsymbol{j}(\boldsymbol{x},t) = (j_x(\boldsymbol{x},t),j_y(\boldsymbol{x},t))$
- $\psi_\alpha(\boldsymbol{x},t) \quad (\alpha=0,1,\dots,N_\psi-1)$

速度場 $\boldsymbol{v}(\boldsymbol{x},t)$ は、$\boldsymbol{j}$ と $\rho$ を用いて以下のように定義されます。

```math
\boldsymbol{v}(\boldsymbol{x},t) = \frac{\boldsymbol{j}(\boldsymbol{x},t)}{\rho(\boldsymbol{x},t)}
```

> [!NOTE]
> スカラー場の成分数 $N_\psi$ は入力スクリプトで指定します。 $N_\psi = 0$ に設定した場合、スカラー場の方程式は除外され、等温の揺らぐNavier-Stokesダイナミクスに帰着します。

## 1.2 Fourier 展開と微分演算

計算領域は2次元周期境界 $\Omega = [0,L_x) \times [0,L_y)$ とします。
空間方向の添字を $a,b,c,d \in \{x,y\}$ とし、同一の添字が繰り返される場合はEinsteinの縮約記法に従います。

任意の場 $f(\boldsymbol{x})$ のFourier展開、および空間微分のFourier空間における定義は以下の通りです。

```math
f(\boldsymbol{x}) = \sum_{\boldsymbol{k}} \hat{f}_{\boldsymbol{k}} \exp(i\boldsymbol{k}\cdot\boldsymbol{x})
```
```math
\partial_a f \longleftrightarrow i k_a \hat{f}_{\boldsymbol{k}}
```
```math
\nabla^2 f \longleftrightarrow -|\boldsymbol{k}|^2 \hat{f}_{\boldsymbol{k}}
```


## 2. 支配方程式（動作モード）

対象とする物理系に応じ、3つの動作モード（圧縮性・非圧縮性・静止流体）を選択可能です。
モードの選択については

- [モード選択](WIKI.md)

を参照してください。


### 2.1 圧縮性モード

密度、運動量、スカラー場の完全なダイナミクスを解きます。支配方程式は以下の通りです。
移流項は入力スクリプトから有効化・無効化できます。

```math
\partial_t \rho = - \nabla\cdot\boldsymbol{j}
```

```math
\partial_t j_a = - \partial_b(j_a v_b) - \nabla_a p + \eta \nabla^2 v_a + \zeta \nabla_a (\nabla \cdot \boldsymbol{v}) + \partial_b\widetilde{\Sigma}_{ab}
```

```math
\partial_t \psi_\alpha = - \nabla\cdot(\psi_\alpha\boldsymbol{v}) + M_\alpha\nabla^2\mu_\alpha + \nabla\cdot\widetilde{\boldsymbol{F}}_\alpha
```


各物性値（移動度 $M_\alpha$、せん断粘性 $\eta$、体積粘性 $\zeta$）は空間的に一様な定数とします。
また、異なるスカラー成分間のOnsager交差係数は $0$ とし、$M_{\alpha\beta} = M_\alpha \delta_{\alpha\beta}$ と仮定しています。

熱揺らぎを表す確率フラックス $\widetilde{\boldsymbol{F}}_\alpha$ と確率応力テンソル $\widetilde{\Sigma}_{ab}$ は、平均0のガウシアン・ホワイトノイズであり、以下の共分散を持ちます。

```math
\left\langle
\widetilde{F}_{\alpha,a}(\boldsymbol{x},t)
\widetilde{F}_{\beta,b}(\boldsymbol{x}',t')
\right\rangle
=
2 k_B T M_\alpha
\delta_{\alpha\beta}\delta_{ab}
\delta(\boldsymbol{x}-\boldsymbol{x}')
\delta(t-t')
```

```math
\left\langle
\widetilde{\Sigma}_{ab}(\boldsymbol{x},t)
\widetilde{\Sigma}_{cd}(\boldsymbol{x}',t')
\right\rangle
=
2 k_B T
\left[
\eta(\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc})
+(\zeta-\eta)\delta_{ab}\delta_{cd}
\right]
\delta(\boldsymbol{x}-\boldsymbol{x}')
\delta(t-t')
```


### 2.2 非圧縮性モード

密度の空間平均 $\rho_0$ を基準密度として固定し（$\rho = \rho_0$）、連続の式を非圧縮条件に置き換えます。

```math
\nabla\cdot\boldsymbol{v} = 0 \quad \left( \boldsymbol{v} = \frac{\boldsymbol{j}}{\rho_0} \right)
```

このモードにおいて圧力 $p$ は、非圧縮条件を満たすためのLagrange乗数として決定されます。支配方程式は以下の通りです。

```math
\partial_t j_a = - \partial_b(j_a v_b) - \nabla_a p + \eta \nabla^2 v_a + \partial_b\widetilde{\Sigma}_{ab}
```

```math
\partial_t \psi_\alpha = - \nabla\cdot(\psi_\alpha\boldsymbol{v}) + M_\alpha\nabla^2\mu_\alpha + \nabla\cdot\widetilde{\boldsymbol{F}}_\alpha
```

実装では、圧力に対するPoisson方程式を明示的に解く代わりに、Fourier空間で運動量の時間増分を横波成分へ射影します。
波数 $\boldsymbol{k}\neq\boldsymbol{0}$ に対して射影演算子を

```math
P_{ab}(\boldsymbol{k}) = \delta_{ab} - \frac{k_a k_b}{|\boldsymbol{k}|^2}
```

と定義すると、運動量方程式の右辺 $R_a$ に対して実際に更新される成分は

```math
\widehat{R}^{\perp}_a(\boldsymbol{k})
=
P_{ab}(\boldsymbol{k})\widehat{R}_b(\boldsymbol{k})
```

です。これにより、圧力勾配のような縦波成分は更新から取り除かれ、$\boldsymbol{k}\cdot\widehat{\boldsymbol{j}}=0$ が保たれます。

### 2.3 静止流体モード

流体のマクロな運動場をゼロ（$\boldsymbol{v}=\boldsymbol{0}$）とし、スカラー場の時間発展のみを解きます。

```math
\partial_t \psi_\alpha = M_\alpha\nabla^2\mu_\alpha + \nabla\cdot\widetilde{\boldsymbol{F}}_\alpha
```

これは保存型の確率的Cahn-Hilliard方程式（Hohenberg-Halperin Model B）に相当します。



## 3. 熱力学量モデル

圧力 $p$ および化学ポテンシャル $\mu_\alpha$ の具体的な関数形は、ユーザーが任意に定義・拡張することが可能です。（参考：機能の追加手順）
すでに実装されているものについては

- [圧力の選択](docs/equations.md)
- [化学ポテンシャルの選択](docs/build.md)

も参照してください。

### 3.1 圧力

等温音速 $c_T$ を用いた、密度 $\rho$ に対する線形な状態方程式を利用可能です。

```math
p(\rho) = c_T^2 \rho
```

また、ユーザー定義により局所的な非線形モデル $p(\rho)$ を組み込むことも可能です。
`非圧縮性モード`では圧力の選択は不要です。


### 3.2 化学ポテンシャル 

化学ポテンシャルは、空間微分の線形項と、局所的な非線形ポテンシャル項 $\mu^{\mathrm{phys}}_\alpha$ の和として定式化されます。

```math
\mu_\alpha = k_{0,\alpha}\psi_\alpha - k_{2,\alpha}\nabla^2\psi_\alpha + k_{4,\alpha}\nabla^4\psi_\alpha + \mu^{\mathrm{phys}}_\alpha(\psi_0,\dots,\psi_{N_\psi-1})
```

Fourier空間における評価式は以下の通りです。線形項は波数ベクトル $\boldsymbol{k}$ に対する代数演算として厳密に評価されます。

```math
\widehat{\mu_\alpha}(\boldsymbol{k}) = \left( k_{0,\alpha} + k_{2,\alpha}|\boldsymbol{k}|^2 + k_{4,\alpha}|\boldsymbol{k}|^4 \right) \widehat{\psi_\alpha}(\boldsymbol{k}) + \widehat{\mu^{\mathrm{phys}}_\alpha}(\boldsymbol{k})
```
