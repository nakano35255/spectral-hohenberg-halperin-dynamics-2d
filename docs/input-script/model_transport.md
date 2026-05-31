# `model transport`

## 目的

粘性係数とスカラー場の移動度を指定します。
このコマンドで指定できる輸送係数モデルは `constant` のみです。
輸送係数は空間的に一様な定数として扱います。

## 形式

```txt
model transport <type> [key value ...]
```

## 例

```sh
model transport constant eta 1.0 zeta 1.0
model transport constant eta 1.0 zeta 1.0 M[0,0] 0.5
```

```sh
model transport constant eta 1.0
model transport coeff zeta 1.0 M[0,0] 0.5
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `constant`
- `[key value ...]`
  - 輸送係数モデルごとの係数を key-value 形式で指定します。

以下の表で `alpha` はスカラー場成分のインデックスを表します。

### `constant`

空間的に一様な粘性係数と移動度を与えるモデルです。

```math
\eta = \mathrm{const.}, \quad
\zeta = \mathrm{const.}, \quad
M_\alpha = \mathrm{const.}
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `eta` | double | no | せん断粘性係数 |
| `shear_viscosity` | double | no | `eta` と同じ意味 |
| `zeta` | double | no | 体積粘性係数 |
| `bulk_viscosity` | double | no | `zeta` と同じ意味 |
| `M[alpha,alpha]` | double | no | 成分 `alpha` の移動度 |

## デフォルト値

`model transport` は必須です。
指定しない場合、入力スクリプトの読み込み時にエラーになります。

各係数を省略した場合、その係数は `0.0` になります。

## coeff コマンド

`model transport coeff` を使用すると、すでに指定済みの輸送係数モデルに対して係数だけを追加・更新できます。
指定した key の値だけが上書きされ、指定しなかった係数は現在の値を保持します。

```sh
model transport constant eta 1.0
model transport coeff zeta 1.0 M[0,0] 0.5
```

この例では、最初の行で `eta` を設定し、次の行で同じ `constant` モデルに `zeta` と `M[0,0]` を追加設定します。
詳細は [`model ... coeff`](./model_coeff.md) を参照してください。

## 詳細

`model transport` は、流体場とスカラー場の散逸項、およびノイズを有効にした場合のノイズ強度に使う輸送係数を定義します。

流体を時間発展する `compressible` / `incompressible` モードでは、`eta` と `zeta` が運動量方程式の粘性項に使われます。

```math
\eta \nabla^2 v_a + \zeta \partial_a(\nabla\cdot\boldsymbol{v})
```

スカラー場が存在する場合、`M[alpha,alpha]` は成分 `alpha` の保存型緩和項に使われます。
[`model free_energy`](./model_free_energy.md) で指定した化学ポテンシャル $\mu_\alpha$ と組み合わせて、

```math
M_\alpha \nabla^2 \mu_\alpha
```

に対応する項が評価されます。

`fix ... noise` を有効にした場合、`eta`, `zeta`, `M[alpha,alpha]` は対応する保存型ノイズの強度にも使われます。スカラー場ノイズには、`fix ... noise` の `chi` で指定する倍率も掛かります。

## 制限・注意

- このコマンドで指定できる輸送係数モデルは `constant` のみです。
- `eta`, `zeta`, `M[alpha,alpha]` は `0` 以上である必要があります。
- `order_parameters` は、`model transport` より前に指定してください。
- 係数 index `alpha` は `0 <= alpha < order_parameters` の範囲である必要があります。
- 現在の実装では、移動度行列は対角成分のみ対応しています。非対角成分 `M[alpha,beta]` は、値が `0.0` の場合だけ受理されます。
- `model transport <type> ...` を再度指定すると、そのカテゴリの設定はデフォルト値から作り直されます。既存係数だけを更新したい場合は [`model ... coeff`](./model_coeff.md) を使ってください。

## 関連コマンド

- [`order_parameters`](./order_parameters.md)
- [`model free_energy`](./model_free_energy.md)
- [`model ... coeff`](./model_coeff.md)
- [`fix ... noise`](./fix_noise.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
