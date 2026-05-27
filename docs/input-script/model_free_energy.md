# `model free_energy`

## 目的

スカラー場（オーダーパラメータ）に対する自由エネルギーモデルを指定します。

## 形式

```txt
model free_energy <type> [key value ...]
```

## 例

```sh
model free_energy quadratic a[0] 1.0
model free_energy phi4 a[0] -3.0 b[0] 1.0 u[0] 1.0
```

```sh
model free_energy phi4 a[0] -3.0 b[0] 1.0
model free_energy coeff u[0] 1.0
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `quadratic`
    - `phi4`
- `[key value ...]`
  - 自由エネルギーモデルごとの係数を key-value 形式で指定します。

以下の表で `alpha` はスカラー場成分のインデックスを表します。

### `quadratic`

線形の化学ポテンシャルを与えるモデルです。

```math
\mu_\alpha
= a_\alpha \psi_\alpha
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `a[alpha]` | double | no | 成分 `alpha` の2次係数 |

### `phi4`

勾配項と4次項を含む化学ポテンシャルを与えるモデルです。

```math
\mu_\alpha
= a_\alpha \psi_\alpha
- b_\alpha \nabla^2 \psi_\alpha
+ u_\alpha \psi_\alpha^3
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `a[alpha]` | double | no | 成分 `alpha` の局所2次係数 |
| `b[alpha]` | double | no | 成分 `alpha` の勾配係数 |
| `u[alpha]` | double | no | 成分 `alpha` の4次係数 |

## デフォルト値

`model free_energy` を省略した場合、化学ポテンシャルを持たない自由エネルギーモデルが使われます。

各係数を省略した場合、その係数は `0.0` になります。

## coeff コマンド

`model free_energy coeff` を使用すると、すでに指定済みの自由エネルギーモデルに対して係数だけを追加・更新できます。
指定した key の値だけが上書きされ、指定しなかった係数は現在の値を保持します。

```sh
model free_energy phi4 a[0] -3.0 b[0] 1.0
model free_energy coeff u[0] 1.0
```

この例では、最初の行で `a[0]` と `b[0]` を設定し、次の行で同じ `phi4` モデルに `u[0]` を追加設定します。
詳細は [`model ... coeff`](./model_coeff.md) を参照してください。

## 詳細

`model free_energy` は、スカラー場の保存型ダイナミクスで使う化学ポテンシャル $\mu_\alpha$ を定義します。
実際の時間発展では、[`model transport`](./model_transport.md) で指定した移動度 $M_{\alpha}$ と組み合わせて、

```math
M_{\alpha} \nabla^2 \mu_\alpha
```

に対応する項が評価されます。

スカラー場の移流項は `model free_energy` ではなく、[`fix ... nonlinear`](./fix_nonlinear.md) で有効化します。

## 制限・注意

- `order_parameters` は、`model free_energy` より前に指定してください。
- 係数 index `alpha` は `0 <= alpha < order_parameters` の範囲である必要があります。
- `u[alpha]` は `0` 以上である必要があります。
- 現在の組み込みモデルでは、成分間の交差項は未実装です。
- `model free_energy <type> ...` を再度指定すると、そのカテゴリの設定はデフォルト値から作り直されます。既存係数だけを更新したい場合は [`model ... coeff`](./model_coeff.md) を使ってください。

## 関連コマンド

- [`order_parameters`](./order_parameters.md)
- [`model transport`](./model_transport.md)
- [`model ... coeff`](./model_coeff.md)
- [`set order_parameter`](./set_order_parameter.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
