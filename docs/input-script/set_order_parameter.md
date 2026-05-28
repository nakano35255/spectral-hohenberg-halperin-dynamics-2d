# `set order_parameter`

## 目的

スカラー場（オーダーパラメータ）$\psi_\alpha(\boldsymbol{x})$ の初期条件を指定します。

## 形式

```txt
set order_parameter <target> <type> [key value ...]
```

`order_parameter` は `psi` と書くこともできます。

```txt
set psi <target> <type> [key value ...]
```

## 例

```sh
set order_parameter 0 uniform value 0.0
set order_parameter all uniform value 0.0
set order_parameter 0 sine base 0.0 amplitude 0.1 nkx 1 nky 0
set order_parameter 0 gaussian base 0.0 amplitude 0.1 x0 32.0 y0 32.0 sigma 4.0
```

## 引数

- `<target>`
  - 型: integer または string
  - 指定可能な値:
    - `0`, `1`, ...
    - `all`
- `<type>`
  - 型: string
  - 指定可能な値:
    - `uniform`
    - `sine`
    - `gaussian`
- `[key value ...]`
  - 初期条件 type ごとの key-value 引数です。

## デフォルト値

スカラー場の初期条件を指定しない場合、各スカラー場はゼロから開始します。

## 詳細

`<target>` には、初期条件を設定するスカラー場成分を指定します。
`all` を指定すると、全てのスカラー場成分に同じ初期条件を適用します。

スカラー場の成分数は [`order_parameters`](./order_parameters.md) で指定します。

### `uniform`

空間一様なスカラー場を設定します。

```sh
set order_parameter 0 uniform value 0.0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `value` | double | yes | 一様なスカラー場の値 |

### `sine`

正弦波状のスカラー場摂動を設定します。

```sh
set order_parameter 0 sine base 0.0 amplitude 0.1 nkx 1 nky 0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `base` | double | no | 背景値。省略時は `0.0` |
| `amplitude` | double | yes | 正弦波の振幅 |
| `nkx` | integer | yes | $x$ 方向の波数 index |
| `nky` | integer | no | $y$ 方向の波数 index。省略時は `0` |

### `gaussian`

Gaussian 状のスカラー場摂動を設定します。

```sh
set order_parameter 0 gaussian base 0.0 amplitude 0.1 x0 32.0 y0 32.0 sigma 4.0
```

`sigma_x` と `sigma_y` を別々に指定することもできます。

```sh
set order_parameter 0 gaussian base 0.0 amplitude 0.1 x0 32.0 y0 32.0 sigma_x 4.0 sigma_y 8.0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `base` | double | yes | 背景値 |
| `amplitude` | double | yes | Gaussian 摂動の振幅 |
| `x0` | double | yes | 中心位置の $x$ 座標 |
| `y0` | double | yes | 中心位置の $y$ 座標 |
| `sigma` | double | no | $x$, $y$ 方向に共通の幅 |
| `sigma_x` | double | no | $x$ 方向の幅 |
| `sigma_y` | double | no | $y$ 方向の幅 |

## 制限・注意

- `<target>` に成分 index を指定する場合、`0 <= target < order_parameters` である必要があります。
- `all` を使う場合、`order_parameters` は `1` 以上である必要があります。
- `sine` の波数は `0 <= nkx < Nx/2`, `-Ny/2 < nky < Ny/2` の範囲で指定します。
- `sine` では `(nkx, nky) = (0, 0)` は指定できません。
- `gaussian` の `sigma`, `sigma_x`, `sigma_y` は正の値である必要があります。
- 複数の order_parameter 初期条件を同じ成分に指定した場合、後から指定した初期条件がその成分全体を上書きします。
- 異なる成分に指定した初期条件は、互いに独立に適用されます。

## 関連コマンド

- [`order_parameters`](./order_parameters.md)
- [`model free_energy`](./model_free_energy.md)
- [`model transport`](./model_transport.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
