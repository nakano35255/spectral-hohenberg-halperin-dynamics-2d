# IMEX Midpoint Integrator

このページでは、`integrator` コマンドで指定できる数値積分スキームのうち、`imex_midpoint` について説明します。

## 目的

`integrator imex_midpoint` を指定すると、時間発展スキームとして **陰的中点予測子・修正子法（Implicit Midpoint Predictor-Corrector method; $w_1=1/2$ ）** が使用されます。

本実装の IMEX は、Delong et al. (2013) で提案された「陰的中点予測子・修正子法（ $w_1=1/2$ ）」に基づいています。音響項や粘性項などの時間刻み制約（CFL条件）を厳しくする線形項を「陰的」に解き、移流などの非線形項を「陽的」に評価することで、計算の安定性を保ちつつ大きな時間刻み幅 $\Delta t$ をとることを目的としています。強い揺らぎの下でも正しい平衡分布のスペクトルを再現できるロバストな手法です。

## 形式

```txt
integrator imex_midpoint
timestep <dt>
```

## 例

```sh
integrator imex_midpoint
timestep   0.01
```

## 引数

`integrator imex_midpoint` 自体は追加の引数を取りません。時間刻み幅は `timestep` で指定します。

## アルゴリズム

時刻 $t_n$ における状態を $u^{(0)} = u_n$ とします。
決定論的な右辺（deterministic RHS）を、陰的に扱う線形項 $\hat{F}_{\mathrm{lin}}(u) = \hat{L}\hat{u}$ と、陽的に扱う非線形項 $\hat{F}_{\mathrm{nonlin}}(u)$ に分割します。
また、1ステップの中で使用する2つの独立な確率論的な右辺（stochastic RHS）を $F_{\mathrm{sto,1}}, F_{\mathrm{sto,2}}$ とします。

本スキームは、予測子（Predictor）と修正子（Corrector）の2つのステージから構成されます。擬スペクトル法を用いているため、線形演算子 $L$ は波数空間において対角化されており、陰的ステップは各波数 $k$ ごとの代数方程式として直接解かれます。

### 1. 予測子（Predictor）ステージ:

時間を $\Delta t/2$ 進めた中間状態 $u^{(1/2)}$ を予測します。線形項には後退オイラー法（完全陰的）、非線形項には前進オイラー法（陽的）が適用されます。

```math
u^{(1/2)}
=
u^{(0)}
+ \frac{\Delta t}{2} F_{\mathrm{lin}}(u^{(1/2)})
+ \frac{\Delta t}{2} F_{\mathrm{nl}}(u^{(0)})
+ \sqrt{\frac{\Delta t}{2}} F_{\mathrm{sto,1}}
```

実装上は、未知数 $u^{(1/2)}$ を左辺に集め、以下の連立方程式（波数空間では除算）として解きます。

```math
\left( I - \frac{\Delta t}{2}L \right) u^{(1/2)}
=
u^{(0)}
+ \frac{\Delta t}{2} F_{\mathrm{nl}}(u^{(0)})
+ \sqrt{\frac{\Delta t}{2}} F_{\mathrm{sto,1}}
```

### 2. 修正子（Corrector）ステージ:

時間を $\Delta t$ 進めた最終状態 $u^{(1)}$ を求めます。線形項にはCrank-Nicolson法（陰的中点則）、非線形項には中点則（陽的）が適用されます。

```math
u^{(1)}
=
u^{(0)}
+ \frac{\Delta t}{2} \left[ F_{\mathrm{lin}}(u^{(0)}) + F_{\mathrm{lin}}(u^{(1)}) \right]
+ \Delta t\, F_{\mathrm{nl}}(u^{(1/2)})
+ \sqrt{\frac{\Delta t}{2}} \left( F_{\mathrm{sto,1}} + F_{\mathrm{sto,2}} \right)
```

これも未知数 $u^{(1)}$ を左辺に集め、以下のように解きます。

```math
\left( I - \frac{\Delta t}{2}L \right) u^{(1)}
=
\left( I + \frac{\Delta t}{2}L \right) u^{(0)}
+ \Delta t\, F_{\mathrm{nl}}(u^{(1/2)})
+ \sqrt{\frac{\Delta t}{2}} \left( F_{\mathrm{sto,1}} + F_{\mathrm{sto,2}} \right)
```

最終的に、次のステップの値を $u_{n+1} = u^{(1)}$ とします。


## 1ステップ全体の更新量

上記の予測子・修正子ステージをまとめると、1ステップ全体での更新量 $\Delta u_n := u_{n+1} - u_n$ は次のように記述できます。

```math
\Delta u_n
=
\Delta t
\left[
\frac{1}{2} F_{\mathrm{lin}}(u^{(0)})
+ \frac{1}{2} F_{\mathrm{lin}}(u^{(1)})
+ F_{\mathrm{nl}}(u^{(1/2)})
\right]
+ \sqrt{\Delta t}
\left[
\frac{1}{\sqrt{2}} F_{\mathrm{sto,1}}
+ \frac{1}{\sqrt{2}} F_{\mathrm{sto,2}}
\right]
```

したがって、Fourier空間における成分 $\alpha$、波数 $k$ についての1ステップの更新量は、次のようにまとまります。

```math
\Delta \hat{\phi}^{\alpha}_{n}(k)
=
\Delta t
\left[
\frac{1}{2} \hat{F}_{\mathrm{lin}}^{\alpha,(0)}(k)
+ \frac{1}{2} \hat{F}_{\mathrm{lin}}^{\alpha,(1)}(k)
+ \hat{F}_{\mathrm{nl}}^{\alpha,(1/2)}(k)
\right]
+ \sqrt{\Delta t} \left[ \frac{1}{\sqrt{2}} \hat{F}_{\mathrm{sto,1}}^{\alpha,n}(k) + \frac{1}{\sqrt{2}} \hat{F}_{\mathrm{sto,2}}^{\alpha,n}(k) \right]
```


## 関連ページ

- [Integrator Commands](./integrator.md)
- [Euler Integrator](./euler.md)