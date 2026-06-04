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
set order_parameter all equilibrium/gaussian psi0 0.0 kBT 1.0 k0 1.0 seed 12345
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
    - `equilibrium/gaussian`
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

### `equilibrium/gaussian`

線形化された化学ポテンシャルに対応する Gaussian 平衡分布から、スカラー場の初期揺らぎを生成します。

```sh
set order_parameter all equilibrium/gaussian psi0 0.0 kBT 1.0 k0 1.0 seed 12345
```

ゼロ mode は `psi0` に固定され、非ゼロ mode のみがランダムに生成されます。
非ゼロ mode の分散は、初期条件コマンドに指定した

```math
\mu_{\rm lin}(\boldsymbol{k})
=
\left(k0 + k2 |\boldsymbol{k}|^2 + k4 |\boldsymbol{k}|^4\right)\hat{\psi}(\boldsymbol{k})
```

に基づいて決まります。
ここで指定する `k0`, `k2`, `k4` は初期条件生成用のパラメータです。
時間発展に使われる [`model free_energy`](./model_free_energy.md) の係数とは自動的には同期しません。

計算格子点数を $N = N_x N_y$、セル面積を $\Delta A = L_x L_y / N$ とします。
ゼロ mode は

```math
\hat{\psi}_\alpha(\boldsymbol{0}) = \psi_0 N
```

に固定されます。
非ゼロ mode ごとに

```math
A(\boldsymbol{k})
=
k0 + k2|\boldsymbol{k}|^2 + k4|\boldsymbol{k}|^4
```

を定義し、各独立な非ゼロ Fourier mode では

```math
\hat{\psi}_\alpha(\boldsymbol{k})
=
\sigma_\psi(\boldsymbol{k})(\xi_1+i\xi_2),
\qquad
\sigma_\psi^2(\boldsymbol{k})
=
\frac{k_B T\,N}{2A(\boldsymbol{k})\Delta A}
```

から生成します。
ここで $\xi_1,\xi_2$ は互いに独立な標準正規乱数です。
非ゼロ mode の分布は

```math
P(\hat{\psi}_\alpha(\boldsymbol{k}))
=
\frac{1}{2\pi\sigma_\psi^2(\boldsymbol{k})}
\exp\left[
-\frac{|\hat{\psi}_\alpha(\boldsymbol{k})|^2}{2\sigma_\psi^2(\boldsymbol{k})}
\right]
```

です。
実空間のスカラー場が実数になるよう、共役な mode は Hermitian 対称性で決まります。

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `psi0` | double | no | 平均値。ゼロ mode に設定される値。省略時は `0.0` |
| `kBT` | double | yes | 揺らぎの強さ |
| `k0` | double | yes | 線形化学ポテンシャルの $k^0$ 係数 |
| `k2` | double | no | 線形化学ポテンシャルの $k^2$ 係数。省略時は `0.0` |
| `k4` | double | no | 線形化学ポテンシャルの $k^4$ 係数。省略時は `0.0` |
| `seed` | integer | no | mode ごとの Gaussian 乱数の seed。省略時は `12345` |

`all` を指定した場合、同じ `psi0`, `kBT`, `k0`, `k2`, `k4`, `seed` が全成分に適用されます。
乱数列は成分ごとに分けられます。
`kBT = 0` の場合は、ゼロ mode のみを持つ一様なスカラー場になります。

## 制限・注意

- `set order_parameter` に指定するパラメータは、スカラー場の初期条件を生成するためのパラメータです。時間発展で使われる `model` や `fix` のパラメータとは自動的には同期しません。
- `<target>` に成分 index を指定する場合、`0 <= target < order_parameters` である必要があります。
- `all` を使う場合、`order_parameters` は `1` 以上である必要があります。
- `sine` の波数は `0 <= nkx < Nx/2`, `-Ny/2 < nky < Ny/2` の範囲で指定します。
- `sine` では `(nkx, nky) = (0, 0)` は指定できません。
- `gaussian` の `sigma`, `sigma_x`, `sigma_y` は正の値である必要があります。
- `equilibrium/gaussian` では `kBT >= 0` である必要があります。
- `equilibrium/gaussian` では、全ての active な非ゼロ mode で `k0 + k2 k^2 + k4 k^4 > 0` である必要があります。
- `equilibrium/gaussian` のランダム成分は非ゼロ mode のみに生成されます。ゼロ mode は `psi0` で指定した平均値に固定されます。
- 複数の order_parameter 初期条件を同じ成分に指定した場合、後から指定した初期条件がその成分全体を上書きします。
- 異なる成分に指定した初期条件は、互いに独立に適用されます。

## 関連コマンド

- [`order_parameters`](./order_parameters.md)
- [`model free_energy`](./model_free_energy.md)
- [`model transport`](./model_transport.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
