# 基礎方程式 (Governing Equations)

本ソルバーは、2次元・等温環境下における $N$ 成分混合流体の揺らぐ流体力学（Fluctuating Hydrodynamics, FHD）方程式を、疑似スペクトル法を用いて高精度に解く汎用ソルバーです。
本ページでは、ソルバー内で実際に時間発展が行われる基礎方程式とその制約条件について定義します。

## 1. 基本変数

本ソルバーは以下の2つの場を独立変数として時間発展させます。

**成分 $k$ の質量密度:**

```math
\rho_k(\boldsymbol{x}, t) \quad (k = 1, \dots, N)
```

**全運動量密度:**

```math
\boldsymbol{j}(\boldsymbol{x}, t)
```

これらを用いて、全密度 $\rho$ と流速（質量中心速度）$\boldsymbol{v}$ を次のように定義します。

```math
\rho := \sum_{k=1}^N \rho_k, \quad \boldsymbol{v} := \frac{\boldsymbol{j}}{\rho}
```


## 2. 支配方程式

### 2.1 保存則

成分質量 $\rho_k(\boldsymbol{x}, t)$ および全運動量 $\boldsymbol{j}(\boldsymbol{x}, t)$ は、以下の保存則に従って時間発展します。


**成分質量の保存:**

```math
\frac{\partial \rho_k}{\partial t} = -\nabla \cdot (\rho_k \boldsymbol{v}) + \nabla \cdot \boldsymbol{F}_k + \nabla \cdot \widetilde{\boldsymbol{F}}_k
```

**全運動量の保存:**

```math
\frac{\partial \boldsymbol{j}}{\partial t} = -\nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) - \nabla p + \nabla \cdot \boldsymbol{\tau} + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```

ここで、方程式の右辺に現れる各項は以下の物理量を表します。
- $\boldsymbol{F}_k, \widetilde{\boldsymbol{F}}_k$: 成分 $k$ の決定論的および確率的拡散フラックス
- $p$: 熱力学的圧力
- $\boldsymbol{\tau}, \widetilde{\boldsymbol{\Sigma}}$: 決定論的および確率的粘性応力テンソル


### 2.2 熱力学と拡散フラックス

すべての熱力学変数は、ヘルムホルツ自由エネルギー密度 $f(\rho_1, \dots, \rho_N)$ から一貫して導出されます。

**化学ポテンシャル:**

```math
\mu_k = \frac{\partial f}{\partial \rho_k}
```

**熱力学的圧力:**

```math
p = -f + \sum_{k=1}^N \rho_k \mu_k
```

**決定論的拡散:**

```math
\boldsymbol{F}_k = \sum_{l=1}^N L_{kl} \nabla \mu_l
```

**【重要】** ユーザーは自由エネルギー $f(\{\rho_k\})$ の表式を指定するだけで済みます。圧力と化学ポテンシャルは、ユーザーが定義した $f$ からソルバー内部で自動生成されるため、熱力学的な整合性（ギブス・デュエムの式）が確実に保証されます。

**【重要】** 本ソルバーでは自由エネルギーに勾配項（例：$\frac{\kappa}{2}|\nabla \rho|^2$）が含まれる場合にも対応するため、圧力項を

```math
    \nabla p = \sum_{k=1}^N \rho_k \nabla \mu_k
```

というギブス・デュエムの式に基づき、化学ポテンシャルから計算しています。これにより、自由エネルギーに勾配項が含まれ、スカラー圧力記述が破綻する（圧力テンソルが導入される）場合にも、正しい方程式を解くことができます。
自由エネルギーに勾配項が含まれる場合、化学ポテンシャルは自動的に変分微分（$\mu_k = \frac{\delta \mathcal{F}}{\delta \rho_k}$）へと拡張されます。ユーザーは複雑なテンソル計算を一切意識することなく、気液・液液相分離などの高度な界面ダイナミクスをシミュレーション可能です。（詳細はSec. 5.2参照）


### 2.3 応力

せん断粘性を $\eta$、体積粘性を $\zeta$ とします。空間のインデックスを $a, b \in \{x, y\}$ とすると、2次元空間における決定論的な粘性応力テンソルの各成分 $\tau_{ab}$ は、以下の一つの式で記述されます。

```math
\tau_{ab} = \eta \left( \frac{\partial v_a}{\partial x_b} + \frac{\partial v_b}{\partial x_a} \right) + (\zeta - \eta) \delta_{ab} (\nabla \cdot \boldsymbol{v})
```


### 2.4　ノイズ

熱揺らぎのモデル化には、揺動散逸定理（Fluctuation-Dissipation Theorem）を満たす平均ゼロのガウシアン・ホワイトノイズを用います。

**【数学的定義】**

確率的フラックス $\widetilde{\boldsymbol{F}}_k$ と確率的応力テンソル $\widetilde{\boldsymbol{\Sigma}}$ の空間的・時間的な相関（共分散構造）は、オンサーガー係数と粘性係数を用いて次のように定義されます（空間インデックス $a,b,c,d \in \{x, y\}$）。

質量フラックスノイズの共分散：

```math
\langle \widetilde{F}_{k,a}(\boldsymbol{x}, t) \widetilde{F}_{l,b}(\boldsymbol{x}', t') \rangle = 2 k_B T L_{kl} \delta_{ab} \delta(\boldsymbol{x} - \boldsymbol{x}') \delta(t - t')
```

応力テンソルノイズの共分散：
```math
\langle \widetilde{\Sigma}_{ab}(\boldsymbol{x}, t) \widetilde{\Sigma}_{cd}(\boldsymbol{x}', t') \rangle = 2 k_B T \left[ \eta (\delta_{ac}\delta_{bd} + \delta_{ad}\delta_{bc}) + (\zeta - \eta) \delta_{ab}\delta_{cd} \right] \delta(\boldsymbol{x} - \boldsymbol{x}') \delta(t - t')
```

ここで、$k_B$ はボルツマン定数、$T$ は絶対温度、$\delta(\cdot)$ はディラックのデルタ関数です。


**【数値計算のための離散化と構成手法】**

実際の数値シミュレーションでは、セル体積を $\Delta V$、時間刻みを $\Delta t$ とし、熱揺らぎの強さを決定するスケーリング係数を $C = \sqrt{\frac{2 k_B T}{\Delta V \Delta t}}$ と定義します。この係数と、独立な標準正規分布に従う乱数を用いて、上記の共分散構造を満たすノイズを構成します。

**質量フラックスノイズの構成:**

移動度行列 $L_{kl}$ に従う相関を持たせるため、これをコレスキー分解した行列を $B$ を導入します。

```math
L_{kl} = \sum_{m=1}^{N-1} B_{km} B_{lm}
```

各成分が独立な標準正規乱数からなるベクトル $\boldsymbol{W}_m$ を用いて、次のように構成します。

```math
\widetilde{\boldsymbol{F}}_k = C \sum_{m=1}^{N-1} B_{km} \boldsymbol{W}_m
```

**応力テンソルノイズの構成:**

テンソルの対称性（$\widetilde{\Sigma}_{xy} = \widetilde{\Sigma}_{yx}$）と半正定値性を担保しつつ上記の共分散を正しく再現するため、3つの独立な標準正規乱数 $\xi_1, \xi_2, \xi_3$ を用いて、等方成分と偏差成分の線形結合として構成します。

```math
\widetilde{\Sigma}_{xx} = C \left( \sqrt{\zeta}\,\xi_1 + \sqrt{\eta}\,\xi_2 \right)
```
```math
\widetilde{\Sigma}_{yy} = C \left( \sqrt{\zeta}\,\xi_1 - \sqrt{\eta}\,\xi_2 \right)
```
```math
\widetilde{\Sigma}_{xy} = \widetilde{\Sigma}_{yx} = C \sqrt{\eta}\,\xi_3
```


### 2.5 質量保存の制約

本方程式系において、各成分の連続の式を足し合わせると、系全体の連続の式（全質量の保存）が得られなければなりません。

```math
\frac{\partial \rho}{\partial t} + \nabla \cdot \boldsymbol{j} = 0
```

この全質量保存が破綻なく成立するには、拡散フラックスの総和が厳密にゼロにならなければなりません。この制約は数式として以下のように記述されます。

```math
\sum_{k=1}^N \boldsymbol{F}_k = \boldsymbol{0} \quad \Rightarrow \quad \sum_{k=1}^N L_{kl} = 0 \quad (\text{すべての } l \text{ について})
```
```math
\sum_{k=1}^N \widetilde{\boldsymbol{F}}_k = \boldsymbol{0}
```

**【重要】** これはモデル化の都合ではなく、流速を「質量中心速度」として定義したことに起因する流体力学的な絶対要請です。本ソルバーを利用する際、行列 $L_{kl}$ に対するこの制約を満たさないパラメータ設定が検出された場合、計算の物理的整合性が保てないためソルバーはエラーを出力して停止します。


### 2.6 パラメータへの追加制約（Additiveノイズ近似）

現在のバージョンの本ソルバーでは、移動度行列 $L_{kl}$、せん断粘性 $\eta$、体積粘性 $\zeta$ はすべて成分密度に依存しない定数であると仮定します。この制約により、2.4節のノイズ項が状態変数（密度）に依存しない「加法的（Additive）ノイズ」となります。※ 将来的に、密度依存性を許容するスキームへの拡張を予定しています）


## 3. 静止流体モード

本ソルバーでは、設定フラグを切り替えることで、計算領域全体の流速場を常にゼロに固定（凍結）して計算を行う「静止流体モード」を利用できます。

```math
    \boldsymbol{v} = \boldsymbol{0}
```

このモードを有効にした場合、方程式系およびソルバーの内部処理は以下のように振る舞います。

* **1. 運動量方程式のバイパス:** ナビエ・ストークス方程式の計算、および各成分の移流項 $\nabla \cdot (\rho_k \boldsymbol{v})$ の計算が完全にスキップされます。これにより、計算コストが大幅に削減されます。
* **2. 拡散方程式系への帰着:** 系のダイナミクスは、熱力学的な化学ポテンシャルの勾配のみを駆動力とする、純粋な多成分拡散方程式へと完全に帰着します。各成分の連続の式は、実質的に以下の形で解かれることになります。

```math
\frac{\partial \rho_k}{\partial t} = \nabla \cdot \left( \sum_{l=1}^N L_{kl} \nabla \mu_l \right) + \nabla \cdot \widetilde{\boldsymbol{F}}_k
```

**【主な利用用途】**
本ソルバーは多成分流体ソルバーとして設計されていますが、このフラグを一つ切り替えるだけで、汎用的な Phase-field（フェーズフィールド）ソルバー として機能します。
流体の流動（移流）を無視できる固体合金のスピノーダル分解（Cahn-Hilliard方程式など、詳細は5.3節参照）などのシミュレーションに利用できます。

## 4. 非圧縮性モード

本ソルバーでは、フラグの切り替えにより「圧縮性流体」と「非圧縮性流体」をシームレスに切り替えてシミュレーションすることが可能です。非圧縮性モードを有効にした場合、系全体の全密度 $\rho$ が空間的・時間的に一定（$\rho = \text{const.}$）であると仮定されます。このとき、全密度の連続の式から、流速場に対する非圧縮条件（発散ゼロ）が要求されます。

```math
    \nabla \cdot \boldsymbol{v} = 0
```

この条件が課されることで、方程式系は以下のように振る舞います。

* **1. 圧力の役割の変化:** 熱力学的圧力 $p$ は状態方程式から決まる変数ではなく、流速の非圧縮条件 $\nabla \cdot \boldsymbol{v} = 0$ を満たすための「ラグランジュの未定乗数」として働きます。
* **2. 応力テンソルの簡略化:** $\nabla \cdot \boldsymbol{v} = 0$ となるため、2.3節の体積粘性 $\zeta$ に依存する項が完全に消滅し、純粋なせん断応力のみが残ります。

**【重要】** 非圧縮モードが選ばれた場合、 体積粘性 $\zeta$ の値は何に設定しても使用されません。


## 5. 例

### 5.1 単成分・（非）圧縮性ナビエ・ストークス方程式 (Single-Component Navier-Stokes)

成分が1種類のみ（$N=1$）の、一般的な以下の（非）圧縮性流体のダイナミクスを記述します。

```math
\frac{\partial \rho}{\partial t} = -\nabla \cdot (\rho \boldsymbol{v})
```
```math
\frac{\partial \boldsymbol{j}}{\partial t} = -\nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) - \nabla p + \eta \nabla^2 \boldsymbol{v} + \zeta \nabla (\nabla \cdot \boldsymbol{v}) + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```

* **設定:** $N=1$ とし、流速場をON。
* **パラメータ:** 2.5節の制約 $\sum L_{kl} = 0$ により、自動的に $L_{11} = 0$ となります。
* **物理的振る舞い:** 拡散フラックスと質量ノイズが自動的にゼロ（$\boldsymbol{F}_1 = \widetilde{\boldsymbol{F}}_1 = \boldsymbol{0}$）となり、方程式系は標準的な「単成分の揺らぐ圧縮性ナビエ・ストークス方程式」へと厳密に一致します。


なお、自由エネルギーを等温理想気体の自由エネルギー
```math
f(\rho) = c_T^2 \rho \left[\log\left(\frac{\rho}{\rho_0}\right) - 1\right]
```
に選ぶことで、圧力を

```math
p = c_T^2 \rho
```

とできます。ここで、$c_T^2 = k_B T / m$ は等温音速の二乗、$\rho_0$ は化学ポテンシャルの基準密度です。


### 5.2 Navier-Stokes-Korteweg（NSK）方程式

Kortewegが1901年に書き下した「気液相共存系（キャビテーションや沸騰など）」を扱う流体方程式の、揺らぐ流体力学に対する拡張です。界面における密度の連続的な変化（Diffuse Interface）から生じる非等方的な応力を考慮します。

```math
    \frac{\partial \rho}{\partial t} + \nabla \cdot (\rho \boldsymbol{v}) = 0
```
```math
\frac{\partial \boldsymbol{j}}{\partial t} = -\nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) - \nabla \cdot \boldsymbol{P}_{\text{NSK}} +  \eta \nabla^2 \boldsymbol{v} + \zeta \nabla (\nabla \cdot \boldsymbol{v}) + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```

* **設定:** $N=1$ とし、流速場をON。
* **パラメータ:** 2.5節の制約 $\sum L_{kl} = 0$ により、自動的に $L_{11} = 0$ となります。
* **熱力学:** 自由エネルギー $f(\rho)$ として、バルクの気液相分離を表すポテンシャルと、界面ペナルティとなる勾配項の和を設定します。

```math
f(\rho) = \frac{A}{2} (\rho - \rho_{\text{gas}})^2 (\rho - \rho_{\text{liquid}})^2 + \frac{\kappa}{2} |\nabla \rho|^2
```

このとき、物理学的に厳密な圧力テンソル $\boldsymbol{P}_{\text{NSK}}$ の形は以下の通りです。

```math
\boldsymbol{P}_{\text{NSK}} = \underbrace{ p_{\text{bulk}} \boldsymbol{I} }_{\text{等方バルク圧力}} - \underbrace{ \left[ \kappa \left( \rho \nabla^2 \rho + \frac{1}{2} |\nabla \rho|^2 \right) \boldsymbol{I} - \kappa (\nabla \rho \otimes \nabla \rho) \right] }_{\boldsymbol{K} : \text{Korteweg応力テンソル（異方性を含む界面張力）}}
```

**【重要】** 本ソルバーでは拡張されたギブス・デュエムの式に基づき、圧力テンソルを導入することなく、化学ポテンシャルを用いて上記の式を計算しています。ユーザーは自由エネルギーを正しく設定するだけNSK方程式の数値計算が可能です。



### 5.3 揺らぐ Cahn-Hilliard 方程式（Hohenberg-Halperin Model B ）

流動を伴わない、固体合金などの純粋な相分離（スピノーダル分解）ダイナミクスを記述する標準モデルです。

```math
    \frac{\partial \phi}{\partial t} = \nabla \cdot \left[ 4M \nabla \frac{\delta \mathcal{F}}{\delta \phi} \right] + 2\nabla \cdot \widetilde{\boldsymbol{F}}
```

ここで、$\mathcal{F}$ は系全体の自由エネルギー汎関数であり、以下のように与えられます。

```math
    \mathcal{F}[\phi] = \int \left[ \frac{A}{4} (\phi^2 - \phi_0^2)^2 + \frac{\kappa}{2} |\nabla \phi|^2 \right] d\boldsymbol{x}
```

このモデルは以下の2成分流体モデルと数学的に等価です。なので、本ソルバーでは以下を解きます。

```math
    \frac{\partial \rho_1}{\partial t} = \nabla \cdot \left[ M \nabla (\mu_1 - \mu_2) \right] + \nabla \cdot \widetilde{\boldsymbol{F}}_1
```
```math
    \frac{\partial \rho_2}{\partial t} = -\nabla \cdot \left[ M \nabla (\mu_1 - \mu_2) \right] - \nabla \cdot \widetilde{\boldsymbol{F}}_1
```

ここで、熱力学を駆動する自由エネルギー密度 $f$ は、$\phi = \rho_1 - \rho_2$ を用いて以下のように設定します。


```math
    f(\rho_1, \rho_2) = \frac{A}{4} (\phi^2 - \phi_0^2)^2 + \frac{\kappa}{2} |\nabla \phi|^2
```

* **設定:** $N=2$ とし、流速場をOFF。これにより移流項と運動量方程式が無視され、純粋な拡散方程式系となります。
* **パラメータ:** 2.5節の質量保存の制約（$\sum L_{kl} = 0$）により、独立な移動度パラメータは1つに決まります。これを $M$ とおくと、自動的に $L_{11} = L_{22} = M$、$L_{12} = L_{21} = -M$ と設定されます。

**【解説】** 上記の設定を行うと、ソルバー内部で化学ポテンシャルが $\mu_1 = \frac{\delta \mathcal{F}}{\delta \rho_1} = \frac{\delta \mathcal{F}}{\delta \phi}$、$\mu_2 = \frac{\delta \mathcal{F}}{\delta \rho_2} = -\frac{\delta \mathcal{F}}{\delta \phi}$ として生成されます。
したがって、駆動力は $\mu_1 - \mu_2 = 2 \frac{\delta \mathcal{F}}{\delta \phi}$ となります。さらに、流速場が OFF であるため全密度 $\rho_1 + \rho_2$ は完全に保存（凍結）され、オーダーパラメータ $\phi = \rho_1 - \rho_2$ の時間発展（$\partial_t \phi = \partial_t \rho_1 - \partial_t \rho_2 = 2\partial_t \rho_1$）を計算すると、係数に $4$ が掛かった教科書通りの Cahn-Hilliard 方程式が厳密に再現されます。


### 5.4 Hohenberg-Halperin Model H

ポリマーブレンドやエマルションなど、流体の流動と相分離（スピノーダル分解）が相互作用するダイナミクスを記述する標準モデルです。オーダーパラメータ $\phi(\boldsymbol{x}, t)$ と流速場 $\boldsymbol{v}(\boldsymbol{x}, t)$ の連成方程式として、教科書では一般に以下のように書き下されます。

```math
    \frac{\partial \phi}{\partial t} + \nabla \cdot (\phi \boldsymbol{v}) = \nabla \cdot \left[ 4M \nabla \frac{\delta \mathcal{F}}{\delta \phi} \right] + 2\nabla \cdot \widetilde{\boldsymbol{F}}
```
```math
    \frac{\partial \boldsymbol{j}}{\partial t} + \nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) = -\nabla p_{\text{bulk}} + \underbrace{ \phi \nabla \frac{\delta \mathcal{F}}{\delta \phi} }_{\text{界面張力}} + \eta \nabla^2 \boldsymbol{v} + \zeta \nabla (\nabla \cdot \boldsymbol{v}) + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```

（※ $\mathcal{F}$ は 5.3 節 Model B と同じ自由エネルギー汎関数）

本ソルバーにおいて以下の2成分の揺らぐ流体モデルとして計算可能です。

```math
    \frac{\partial \rho_1}{\partial t} + \nabla \cdot (\rho_1 \boldsymbol{v}) = \nabla \cdot \left[ M \nabla (\mu_1 - \mu_2) \right] + \nabla \cdot \widetilde{\boldsymbol{F}}_1
```
```math
    \frac{\partial \rho_2}{\partial t} + \nabla \cdot (\rho_2 \boldsymbol{v}) = -\nabla \cdot \left[ M \nabla (\mu_1 - \mu_2) \right] - \nabla \cdot \widetilde{\boldsymbol{F}}_1
```
```math
    \frac{\partial \boldsymbol{j}}{\partial t} = -\nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) - \nabla \cdot \boldsymbol{P}_{\text{ModelH}} + \eta \nabla^2 \boldsymbol{v} + \zeta \nabla (\nabla \cdot \boldsymbol{v}) + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```

* **設定:** $N=2$ とし、流速場をONにします。
* **パラメータ:** 5.3 節（Model B）と全く同じです（$\sum L_{kl} = 0$ により独立パラメータは $M$ のみ）。
* **熱力学:** 全密度 $\rho = \rho_1 + \rho_2$ の揺らぎ（音波）を抑制して非圧縮性に近づけるため、5.3 節の自由エネルギーにベースとなる等方圧力の項を追加して設定します。

```math
    f(\rho_1, \rho_2) = \underbrace{ c_T^2 \rho \ln \rho }_{\text{全密度の安定化}} + \underbrace{ \frac{A}{4} (\phi^2 - \phi_0^2)^2 + \frac{\kappa}{2} |\nabla \phi|^2 }_{\text{相分離と界面のペナルティ}} \quad \text{（ただし } \phi = \rho_1 - \rho_2 \text{）}
```

このとき、物理学的に厳密な圧力テンソル $\boldsymbol{P}_{\text{ModelH}}$ の形は以下の通りです。

```math
\boldsymbol{P}_{\text{ModelH}} = \underbrace{ p_{\text{bulk}} \boldsymbol{I} }_{\text{等方バルク圧力}} - \underbrace{ \left[ \kappa \left( \phi \nabla^2 \phi + \frac{1}{2} |\nabla \phi|^2 \right) \boldsymbol{I} - \kappa (\nabla \phi \otimes \nabla \phi) \right] }_{\text{毛管応力（界面張力）テンソル}}
```

**【重要】**  本ソルバーでは拡張されたギブス・デュエムの式に基づき、複雑な圧力テンソルをプログラム上で導入することなく、化学ポテンシャル（$-\rho_1 \nabla \mu_1 - \rho_2 \nabla \mu_2$）を用いて自動的に上記の力を計算しています。ユーザーは自由エネルギーを正しく設定するだけで Model H の数値計算が可能です。

**【解説：非圧縮性の取り扱い】**
一般に、教科書的な Model H は全密度 $\rho$ が空間的・時間的に一定である「非圧縮性流体（$\nabla \cdot \boldsymbol{v} = 0$）」として定式化されており、流速場を非圧縮に保つための圧力（ポアソン方程式）を解く必要があります。
本ソルバーの非圧縮モードでは、射影法（Fractional Step法等）を用いてこのポアソン方程式を解きます。この際、右辺にまるごと代入された駆動力（$- \rho_1 \nabla \mu_1 - \rho_2 \nabla \mu_2$）のうち、純粋な勾配成分である等方的なバルク圧力（$-\nabla p_{\text{bulk}}$）は、ヘルムホルツ分解の性質により自動的にポアソン方程式の擬似圧力に吸収・相殺されます。
