# Euler Integrator

このページでは、`integrator` コマンドで指定できる数値積分スキームのうち、`euler` について説明します。

## 目的

`integrator euler` を指定すると、時間発展スキームとして **Euler-Maruyama法** が使用されます。

## 形式

```txt
integrator euler
timestep <dt>
```

## 例

```sh
integrator euler
timestep   0.005
```

## 引数

`integrator euler` 自体は追加の引数を取りません。時間刻み幅は `timestep` で指定します。

## アルゴリズム

時刻 $t_n = n \Delta t$ における状態を $u_n$ とします。
決定論的な右辺（deterministic RHS）を $F_{\mathrm{det}}(u_n,t_n)$、確率論的な右辺（stochastic RHS）を $F_{\mathrm{sto},n}$ と定義すると、Euler-Maruyama法の更新式（update）は以下のようになります。

```math
u_{n+1}
=
u_n
+ \Delta t\,F_{\mathrm{det}}(u_n,t_n)
+ \sqrt{\Delta t}\,F_{\mathrm{sto},n}
```

これをFourier空間における成分 $\alpha$、波数 $k$ について記述すると、以下のようになります

```math
\hat{\phi}^{\alpha}_{n+1}(k)
=
\hat{\phi}^{\alpha}_{n}(k)
+ \Delta t\,\hat{F}^{\alpha,n}_{\mathrm{det}}(k)
+ \sqrt{\Delta t}\,\hat{F}^{\alpha,n}_{\mathrm{sto}}(k)
```


## 関連ページ

- [Integrator Commands](./integrator.md)
- [SRK3 Integrator](./srk3.md)
