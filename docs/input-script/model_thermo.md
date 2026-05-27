# `model thermo`

## 目的

密度場から圧力項を決める熱力学モデルを指定します。

## 形式

```txt
model thermo <type> [key value ...]
```

## 例

```sh
model thermo linear_eos cT 10.0
model thermo linear_eos sound_speed 10.0
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `linear_eos`
- `cT <value>`
  - 型: double
  - 等温音速を指定します。
  - `sound_speed <value>` も同じ意味で使えます。

## デフォルト値

`model thermo` を省略した場合、圧力項を持たない熱力学モデルが使われます。

`model thermo linear_eos` を指定し、`cT` を省略した場合は以下が使われます。

```txt
cT 0.0
```

## 詳細

`linear_eos` は、圧力の線形係数として $c_T^2$ を使う状態方程式です。
実装上は、Fourier 空間で

```math
\widehat{p}_{\mathrm{lin}}(\boldsymbol{k})
= c_T^2 \widehat{\rho}(\boldsymbol{k})
```

に対応する圧力勾配項を運動量方程式へ加えます。

```math
- \nabla_a p_{\mathrm{lin}}
```

`incompressible` モードでは、圧力の勾配成分は横波射影によって除去されます。
そのため、非圧縮性流体の拘束圧力を `model thermo` で明示的に指定する必要はありません。

## 制限・注意

- 現在指定できる熱力学モデルは `linear_eos` のみです。
- `cT` は `0` 以上である必要があります。
- 非線形の密度依存圧力は、現在の実装では未実装です。

## 関連コマンド

- [`time_evolution`](./time_evolution.md)
- [`model ... coeff`](./model_coeff.md)
- [`model transport`](./model_transport.md)
