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
set density equilibrium/gaussian rho0 1.0 kBT 1.0 cT 10.0 seed 12345
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `uniform`
    - `sine`
    - `gaussian`
    - `equilibrium/gaussian`
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

### `equilibrium/gaussian`

線形圧力に対応する Gaussian 平衡分布から、密度場の初期揺らぎを生成します。

```sh
set density equilibrium/gaussian rho0 1.0 kBT 1.0 cT 10.0 seed 12345
```

ゼロ mode は平均密度 `rho0` に固定され、非ゼロ mode のみがランダムに生成されます。
分散は初期条件コマンドに指定した `cT` から決まる線形圧力係数 $c_T^2$ を使います。
この `cT` は初期条件生成用のパラメータであり、時間発展に使われる [`model thermo`](./model_thermo.md) の `cT` を自動参照しません。

計算格子点数を $N = N_x N_y$、セル面積を $\Delta A = L_x L_y / N$ とします。
ゼロ mode は

```math
\hat{\rho}_{\boldsymbol{0}} = \rho_0 N
```

に固定されます。
各独立な非ゼロ Fourier mode では、$\xi_1,\xi_2$ を互いに独立な標準正規乱数として

```math
\hat{\rho}_{\boldsymbol{k}}
=
\sigma_\rho(\xi_1 + i\xi_2),
\qquad
\sigma_\rho^2
=
\frac{k_B T\,N\,\rho_0}{2c_T^2\Delta A}
```

から生成します。
したがって、非ゼロ mode の分布は

```math
P(\hat{\rho}_{\boldsymbol{k}})
=
\frac{1}{2\pi\sigma_\rho^2}
\exp\left[
-\frac{|\hat{\rho}_{\boldsymbol{k}}|^2}{2\sigma_\rho^2}
\right]
```

です。
実空間の密度場が実数になるよう、共役な mode は Hermitian 対称性で決まります。

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `rho0` | double | yes | 平均密度。ゼロ mode に設定される値 |
| `kBT` | double | yes | 揺らぎの強さ |
| `cT` | double | yes | 線形圧力係数を決める等温音速 |
| `seed` | integer | no | mode ごとの Gaussian 乱数の seed。省略時は `12345` |

`kBT = 0` の場合は、ゼロ mode のみを持つ一様密度になります。

## 制限・注意

- `set density` に指定するパラメータは、密度場の初期条件を生成するためのパラメータです。時間発展で使われる `model` や `fix` のパラメータとは自動的には同期しません。
- density には成分 index や `all` target はありません。
- `set density 0 ...` や `set density all ...` という形式は現在の実装では使えません。
- `sine` の波数は `0 <= nkx < Nx/2`, `-Ny/2 < nky < Ny/2` の範囲で指定します。
- `sine` では `(nkx, nky) = (0, 0)` は指定できません。
- `equilibrium/gaussian` では `rho0 > 0`, `kBT >= 0`, `cT > 0` である必要があります。
- `equilibrium/gaussian` のランダム成分は非ゼロ mode のみに生成されます。ゼロ mode は `rho0` で指定した平均密度に固定されます。
- 複数の density 初期条件を指定した場合、後から指定した初期条件が密度場全体を上書きします。

## 関連コマンド

- [`grid`](./grid.md)
- [`length`](./length.md)
- [`time_evolution`](./time_evolution.md)
- [`set momentum`](./set_momentum.md)
