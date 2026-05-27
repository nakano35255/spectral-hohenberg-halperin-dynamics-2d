# `set density`

## 目的

密度場 $\rho(\boldsymbol{x})$ の初期条件を指定します。

## 形式

```txt
set density <type> [key value ...]
```

`density` は `rho` と書くこともできます。

```txt
set rho <type> [key value ...]
```

## 例

```sh
set density uniform value 1.0
set density sine base 1.0 amplitude 0.01 nkx 1 nky 0
set density gaussian base 1.0 amplitude 0.1 x0 32.0 y0 32.0 sigma 4.0
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `uniform`
    - `sine`
    - `gaussian`
- `[key value ...]`
  - 初期条件 type ごとの key-value 引数です。

## デフォルト値

密度場の初期条件を指定しない場合、密度場はゼロから開始します。

ただし、`compressible` / `incompressible` モードでは速度場の定義に密度を使うため、通常は非ゼロの密度初期条件を明示的に指定してください。

## 詳細

### `uniform`

空間一様な密度場を設定します。

```sh
set density uniform value 1.0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `value` | double | yes | 一様密度の値 |

`value` は `0` 以上である必要があります。

### `sine`

正弦波状の密度摂動を設定します。

```sh
set density sine base 1.0 amplitude 0.01 nkx 1 nky 0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `base` | double | yes | 背景密度 |
| `amplitude` | double | yes | 正弦波の振幅 |
| `nkx` | integer | yes | $x$ 方向の波数 index |
| `nky` | integer | no | $y$ 方向の波数 index。省略時は `0` |

密度が負にならないよう、`base >= abs(amplitude)` である必要があります。

### `gaussian`

ガウス型の密度分布を設定します。

```sh
set density gaussian base 1.0 amplitude 0.1 x0 32.0 y0 32.0 sigma 4.0
```

異方的な幅を指定する場合は、`sigma_x` と `sigma_y` を使います。

```sh
set density gaussian base 1.0 amplitude 0.1 x0 32.0 y0 32.0 sigma_x 4.0 sigma_y 8.0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `base` | double | yes | 背景密度 |
| `amplitude` | double | yes | ガウス分布の振幅 |
| `x0` | double | yes | 中心位置の $x$ 座標 |
| `y0` | double | yes | 中心位置の $y$ 座標 |
| `sigma` | double | conditional | 等方的な幅 |
| `sigma_x` | double | conditional | $x$ 方向の幅 |
| `sigma_y` | double | conditional | $y$ 方向の幅 |

`sigma` を指定した場合、`sigma_x = sigma_y = sigma` として扱われます。
`sigma` を指定しない場合は、`sigma_x` と `sigma_y` の両方が必要です。

## 制限・注意

- density には成分 index や `all` target はありません。
- `set density 0 ...` や `set density all ...` という形式は現在の実装では使えません。
- `sine` の波数は `0 <= nkx < Nx/2`, `-Ny/2 < nky < Ny/2` の範囲で指定します。
- `sine` では `(nkx, nky) = (0, 0)` は指定できません。
- 複数の density 初期条件を指定した場合、後から指定した初期条件が密度場全体を上書きします。

## 関連コマンド

- [`grid`](./grid.md)
- [`length`](./length.md)
- [`time_evolution`](./time_evolution.md)
- [`set momentum`](./set_momentum.md)
