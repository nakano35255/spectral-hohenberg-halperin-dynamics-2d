# SRK3 Integrator

このページでは、`integrator` コマンドで指定できる数値積分スキームのうち、`srk3` について説明します。

## 目的

`integrator srk3` を指定すると、時間発展スキームとして **3段の確率的Runge-Kutta法（Stochastic Runge-Kutta method; SRK3）** が使用されます。

本実装の SRK3 は、決定論的項を3つのステージで評価し、Garciaらによる係数に基づいて確率論的項を組み合わせる手法です。通常のシミュレーション（プロダクションラン）では、精度と安定性の観点からこの `srk3` の使用を強く推奨します（省略時のデフォルト設定でもあります）。

## 形式

```txt
integrator srk3
timestep <dt>
```

## 例

```sh
integrator srk3
timestep   0.01
```

## 引数

`integrator srk3` 自体は追加の引数を取りません。時間刻み幅は `timestep` で指定します。

## アルゴリズム

時刻 $t_n$ における状態を $u^{(0)} = u_n$ とします。決定論的な右辺（deterministic RHS）を $F_{\mathrm{det}}(u,t)$ とします。また、1ステップの中で使用する2つの独立な確率論的な右辺（stochastic RHS）を $F_{\mathrm{sto,A}}, F_{\mathrm{sto,B}}$ とします。

各ステージで評価されるノイズ項 $S_i$ は以下の形で与えられます。

```math
S_i
=
F_{\mathrm{sto,A}}
+ \beta_i F_{\mathrm{sto,B}}
```

実装に使用しているGarcia係数 $\beta_i$ は以下の通りです。

```math
\beta_1 = \frac{2\sqrt{2}+\sqrt{3}}{5}, \quad
\beta_2 = \frac{-4\sqrt{2}+3\sqrt{3}}{5}, \quad
\beta_3 = \frac{\sqrt{2}-2\sqrt{3}}{10}
```

各ステージの更新式（stage update）は以下のように計算されます。

```math
u^{(1)}
=
u^{(0)}
+ \Delta t\,F_{\mathrm{det}}(u^{(0)},t_n)
+ \sqrt{\Delta t}\,S_1
```
```math
u^{(2)}
=
\frac{3}{4}u^{(0)}
+ \frac{1}{4}
\left[
u^{(1)}
+ \Delta t\,F_{\mathrm{det}}(u^{(1)},t_n+\Delta t)
+ \sqrt{\Delta t}\,S_2
\right]
```
```math
u^{(3)}
=
\frac{1}{3}u^{(0)}
+ \frac{2}{3}
\left[
u^{(2)}
+ \Delta t\,F_{\mathrm{det}}\left(u^{(2)},t_n+\frac{1}{2}\Delta t\right)
+ \sqrt{\Delta t}\,S_3
\right]
```

最終的に、次のステップの値を $u_{n+1} = u^{(3)}$ とします。


## 1ステップ全体の更新量

上記の3つのステージをまとめると、1ステップ全体での更新量 $\Delta u_n := u_{n+1} - u_n$ は次のように記述できます。

```math
\Delta u_n
=
\Delta t
\left[
\frac{1}{6}F_{\mathrm{det}}(u^{(0)},t_n)
+ \frac{1}{6}F_{\mathrm{det}}(u^{(1)},t_n+\Delta t)
+ \frac{2}{3}F_{\mathrm{det}}\left(u^{(2)},t_n+\frac{1}{2}\Delta t\right)
\right]
+ \sqrt{\Delta t}
\left[
\frac{1}{6}S_1
+ \frac{1}{6}S_2
+ \frac{2}{3}S_3
\right]
```

Garcia係数は $\frac{1}{6}\beta_1 + \frac{1}{6}\beta_2 + \frac{2}{3}\beta_3 = 0$ を満たすように設計されているため、確率論的項については以下の関係が成り立ちます。

```math
\frac{1}{6}S_1
+ \frac{1}{6}S_2
+ \frac{2}{3}S_3
=
F_{\mathrm{sto,A}}
```

したがって、Fourier空間における成分 $\alpha$、波数 $k$ についての1ステップの更新量は、次のようにまとまります。

```math
\Delta \hat{\phi}^{\alpha}_{n}(k)
=
\Delta t
\left[
\frac{1}{6}\hat{F}_{\mathrm{det}}^{\alpha,(0)}(k)
+ \frac{1}{6}\hat{F}_{\mathrm{det}}^{\alpha,(1)}(k)
+ \frac{2}{3}\hat{F}_{\mathrm{det}}^{\alpha,(2)}(k)
\right]
+ \sqrt{\Delta t}\,\hat{F}_{\mathrm{sto,A}}^{\alpha,n}(k)
```




## 関連ページ

- [Integrator Commands](./integrator.md)
- [Euler Integrator](./euler.md)
