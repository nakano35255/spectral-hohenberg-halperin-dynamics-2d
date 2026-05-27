# `model transport`

## 目的

粘性係数とスカラー場の移動度を指定します。

## 形式

```txt
model transport <type> [key value ...]
```

## 例

```sh
model transport constant eta 1.0 zeta 1.0
model transport constant eta 1.0 zeta 1.0 M[0,0] 0.5
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `constant`
- `eta <value>`
  - 型: double
  - せん断粘性係数を指定します。
  - `shear_viscosity <value>` も同じ意味で使えます。
- `zeta <value>`
  - 型: double
  - 体積粘性係数を指定します。
  - `bulk_viscosity <value>` も同じ意味で使えます。
- `M[q,q] <value>`
  - 型: double
  - スカラー場成分 `q` の移動度を指定します。

## デフォルト値

`model transport` は必須です。
指定しない場合、入力スクリプトの読み込み時にエラーになります。

`model transport constant` を指定し、係数を省略した場合は以下が使われます。

```txt
eta 0.0
zeta 0.0
M[q,q] 0.0
```

## 詳細

`constant` は、空間的に一定の輸送係数を使うモデルです。

流体モードでは、`eta` と `zeta` が運動量方程式の粘性項および運動量ノイズの強度に使われます。

スカラー場が存在する場合、`M[q,q]` は成分 `q` の保存型緩和項とスカラー場ノイズの強度に使われます。
概念的には、[`model free_energy`](./model_free_energy.md) で指定した化学ポテンシャルと組み合わせて

```math
\nabla_a M_{qq} \nabla_a \mu_q
```

に対応する項を評価します。

## 制限・注意

- 現在指定できる輸送係数モデルは `constant` のみです。
- `eta`, `zeta`, `M[q,q]` は `0` 以上である必要があります。
- `order_parameters` は、`model transport` より前に指定してください。
- 係数 index `q` は `0 <= q < order_parameters` の範囲である必要があります。
- 現在の実装では、移動度行列は対角成分のみ対応しています。非対角成分 `M[q,r]` は、値が `0.0` の場合だけ受理されます。
- `model transport <type> ...` を再度指定すると、そのカテゴリの設定はデフォルト値から作り直されます。既存係数だけを更新したい場合は [`model ... coeff`](./model_coeff.md) を使ってください。

## 関連コマンド

- [`order_parameters`](./order_parameters.md)
- [`model free_energy`](./model_free_energy.md)
- [`model ... coeff`](./model_coeff.md)
- [`fix ... noise`](./fix_noise.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
