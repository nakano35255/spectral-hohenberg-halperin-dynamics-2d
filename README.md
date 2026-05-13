# Spectral Fluctuating Isothermal Fluid

本リポジトリは、2次元・等温・多成分混合流体の「揺らぐ流体力学（Fluctuating Hydrodynamics, FHD）方程式」を解くための数値計算ソルバーです。熱揺らぎを含む流体の挙動を、周期境界条件を持つ長方形領域 $\Omega = [0, L_x) \times [0, L_y)$ において疑似スペクトル法を用いて高精度にシミュレーションすることを目的としています。

本ソルバーのアーキテクチャは、理想気体から、強い非理想性を持つ実在液体混合物まで、幅広い非線形な物性モデルをシームレスに拡張できるよう設計されています。


## 支配方程式 (Governing Equations)

本ソルバーは、$N$ 成分の質量密度 $\rho_k(\boldsymbol{x}, t)$ および、全運動量密度 $\boldsymbol{j}(\boldsymbol{x}, t)$ を基礎変数として時間発展を行います。
全密度を $\rho = \sum_{k=1}^N \rho_k$、流速を $\boldsymbol{v} = \frac{\boldsymbol{j}}{\rho}$ と定義します。

### 1. 連続の式（成分質量の保存）

各成分 $k$ の密度は、移流と拡散（決定論的および確率的）によって変化します。

```math
\frac{\partial \rho_k}{\partial t} = -\nabla \cdot (\rho_k \boldsymbol{v}) + \nabla \cdot \boldsymbol{F}_k + \nabla \cdot \widetilde{\boldsymbol{F}}_k
```

### 2. 運動量方程式

全運動量密度は、移流、圧力勾配、および粘性応力（決定論的および確率的）によって変化します。

```math
\frac{\partial \boldsymbol{j}}{\partial t} = -\nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) - \nabla p + \nabla \cdot \boldsymbol{\tau} + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```

## 物理モデルの詳細

### 熱力学と状態方程式 (Thermodynamics & EOS)

すべての熱力学変数は、体積あたりのヘルムホルツ自由エネルギー密度 $f(\rho_1, \dots, \rho_N)$ から一貫して導出されます。これにより、ギブス・デュエムの式 $\nabla p = \sum_{k=1}^N \rho_k \nabla \mu_k$ が満たされ、熱力学的整合性が保証されます。

- 化学ポテンシャル: $\mu_k = \frac{\partial f}{\partial \rho_k}$
- 熱力学的圧力: $p = -f + \sum_{k=1}^N \rho_k \mu_k$

**【ベースラインモデル：理想気体混合物】**

等温の理想気体混合物（$k_B$: ボルツマン定数, $T$: 絶対温度, $m_k$: 分子質量, $\rho_{k0}$: 基準密度）をデフォルトとしています。

```math
f = k_B T \sum_{k=1}^N \frac{\rho_k}{m_k} \left[ \ln\left(\frac{\rho_k}{\rho_{k0}}\right) - 1 \right]
```

対応する化学ポテンシャルと圧力は次のように与えられます。

```math
\mu_k = \frac{k_B T}{m_k} \ln\left(\frac{\rho_k}{\rho_{k0}}\right)
```

```math
p = k_B T \sum_{k=1}^N \frac{\rho_k}{m_k}
```

**【拡張性：非理想液体への対応】**

アルコール水溶液などの実在液体を扱う場合、自由エネルギーに Redlich-Kister 多項式などの「過剰自由エネルギー（Excess Free Energy）」を追加することで、活量係数を伴う非線形な熱力学モデルへと容易に拡張可能です。


### 成分拡散フラックス (Species Diffusion)

化学ポテンシャル勾配を駆動力とし、オンサーガーの移動度行列 $L_{kl}(\boldsymbol{\rho})$ を介して拡散フラックスが決まります。質量保存則を満たすため、$\sum_{k=1}^N L_{kl} = 0$ が課されます。

```math
\boldsymbol{F}_k = \sum_{l=1}^N L_{kl} \nabla \mu_l
```

ここでは、連続の式に現れる $+\nabla \cdot \boldsymbol{F}_k$ に合わせて
$\boldsymbol{F}_k$ の符号を定義しています。

**【ベースラインモデル：2成分理想拡散】**

相互拡散係数 $D$ を用いたベースラインは以下の通りです。

```math
L_{11} = L_{22} = -L_{12} = -L_{21} = \frac{m_1 m_2 \rho_1 \rho_2}{k_B T (m_1 \rho_2 + m_2 \rho_1)} D
```

【拡張性：実在液体の複雑な拡散】

非理想液体では、熱力学的因子（活量係数の対数微分）や、Vignes則に基づく指数関数的な拡散係数の変化が現れます。本ソルバーは物性値評価を独立させているため、これらの複雑な対数・指数・べき乗モデルをそのまま組み込むことができます。

### 粘性応力 (Viscous Stress)

せん断粘性を $\eta(\boldsymbol{\rho})$、体積粘性を $\zeta(\boldsymbol{\rho})$ とします。2次元における決定論的な粘性応力テンソルは以下の通りです。

```math
\tau_{xx} = 2 \eta \frac{\partial v_x}{\partial x} + (\zeta - \eta) (\nabla \cdot \boldsymbol{v})
```

```math
\tau_{yy} = 2 \eta \frac{\partial v_y}{\partial y} + (\zeta - \eta) (\nabla \cdot \boldsymbol{v})
```

```math
\tau_{xy} = \tau_{yx} = \eta \left( \frac{\partial v_x}{\partial y} + \frac{\partial v_y}{\partial x} \right)
```

**【ベースラインモデル：単原子理想気体】**

気体分子運動論に基づき、等温理想気体におけるせん断粘性は密度に依存しない定数 $\eta_0$ とし、体積粘性はゼロ（ストークスの仮説）とします。$$\eta(\boldsymbol{\rho}) = \eta_0, \quad \zeta(\boldsymbol{\rho}) = 0$$

**【拡張性：実在液体の非線形粘性】**

ガラス形成液体やアルコール水溶液などでは、粘性が密度や濃度に対して指数関数的（Doolittle式）や対数的（Grunberg-Nissan式）に増大します。本アーキテクチャでは、これらの極めて強い非線形性も物理空間で安定して評価可能です。

### 確率的フラックスと応力ノイズ (Stochastic Terms)

熱力学第二法則（半正定値性）を満たすよう、独立な標準正規分布ノイズ $\boldsymbol{W}_m$ および $\xi$ を用いて構成します。スケーリング係数を $C = \sqrt{\frac{2 k_B T}{\Delta V \Delta t}}$ とします（$\Delta V$: セル体積, $\Delta t$: 時間刻み）。

**質量フラックスノイズ:**

コレスキー分解行列 $B$ （$L_{kl} = \sum_m B_{km} B_{lm}$）を用いて構成します。

```math
\widetilde{\boldsymbol{F}}_k = C \sum_{m=1}^{N-1} B_{km} \boldsymbol{W}_m
```

**応力テンソルノイズ:**

$\eta \ge 0, \zeta \ge 0$ の熱力学的制約下で半正定値性が明示されるよう、等方成分と偏差成分に分解して構成します。
 
```math
\widetilde{\Sigma}_{xx} = C \left( \sqrt{\zeta}\,\xi_1 + \sqrt{\eta}\,\xi_2 \right)
```

```math
\widetilde{\Sigma}_{yy} = C \left( \sqrt{\zeta}\,\xi_1 - \sqrt{\eta}\,\xi_2 \right)
```

```math
\widetilde{\Sigma}_{xy} = \widetilde{\Sigma}_{yx} = C \sqrt{\eta}\,\xi_3
```
