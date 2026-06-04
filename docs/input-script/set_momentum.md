# `set momentum`

## 目的

運動量密度場 $\boldsymbol{j}(\boldsymbol{x}) = (j_x, j_y)$ の初期条件を指定します。

## 形式

```txt
set momentum <target> <type> [key value ...]
```

`momentum` は `j` と書くこともできます。

```txt
set j <target> <type> [key value ...]
```

## 例

```sh
set momentum x uniform value 0.0
set momentum y uniform value 0.0
set momentum all uniform value 0.0
set momentum x sine base 0.0 amplitude 1.0 nkx 1 nky 0
set momentum all random_vorticity n0 8 sigma 2 urms 1.0 seed 12345
set momentum all equilibrium/gaussian/compressible kBT 1.0 seed 12345
set momentum all equilibrium/gaussian/incompressible kBT 1.0 seed 12345
```

## 引数

- `<target>`
  - 型: string
  - 指定可能な値:
    - `x`
    - `y`
    - `0`
    - `1`
    - `all`
- `<type>`
  - 型: string
  - 指定可能な値:
    - `uniform`
    - `sine`
    - `random_vorticity`
    - `equilibrium/gaussian/compressible`
    - `equilibrium/gaussian/incompressible`
- `[key value ...]`
  - 初期条件 type ごとの key-value 引数です。

## デフォルト値

運動量場の初期条件を指定しない場合、運動量場はゼロから開始します。

## 詳細

`<target>` は、初期条件を設定する運動量成分を指定します。

| target | 対象 |
| --- | --- |
| `x`, `0` | $j_x$ |
| `y`, `1` | $j_y$ |
| `all` | $j_x$ と $j_y$ の両方 |

### `uniform`

空間一様な運動量場を設定します。

```sh
set momentum x uniform value 0.0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `value` | double | yes | 一様な運動量密度 |

### `sine`

正弦波状の運動量摂動を設定します。

```sh
set momentum x sine base 0.0 amplitude 1.0 nkx 1 nky 0
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `base` | double | yes | 背景値 |
| `amplitude` | double | yes | 正弦波の振幅 |
| `nkx` | integer | yes | $x$ 方向の波数 index |
| `nky` | integer | no | $y$ 方向の波数 index。省略時は `0` |

### `random_vorticity`

ランダム位相を持つ narrow-band な渦度場から、非圧縮条件を満たす運動量場を設定します。

```sh
set momentum all random_vorticity n0 8 sigma 2 urms 1.0 seed 12345
```

この初期条件では、Fourier mode index

```math
\boldsymbol{n} = (n_x, n_y)
```

を使って、各 mode の物理波数を

```math
\boldsymbol{k}_{\boldsymbol{n}} = \left(\frac{2\pi}{L_x} n_x, \frac{2\pi}{L_y} n_y \right)
```

と定義します。index ${\boldsymbol{n}}$ での渦度の Fourier 成分は

```math
\hat{\omega}_{\boldsymbol{n}} = A \exp\left[ -\frac{(|\boldsymbol{n}|-n_0)^2}{2\sigma^2} \right] e^{i\theta_{\boldsymbol{n}}}
```

として作られます。$\theta_{\boldsymbol{n}}$ は mode ごとのランダム位相です。その後、非ゼロ mode に対して

```math
\hat{u}_x(\boldsymbol{k}) = i\frac{k_y}{|\boldsymbol{k}|^2}\hat{\omega}_{\boldsymbol{n}},
\qquad
\hat{u}_y(\boldsymbol{k}) = -i\frac{k_x}{|\boldsymbol{k}|^2}\hat{\omega}_{\boldsymbol{n}}
```

から速度場を作り、

```math
\hat{\boldsymbol{j}}(\boldsymbol{k}) = \rho_0 \hat{\boldsymbol{u}}(\boldsymbol{k})
```

を設定します。$\rho_0$ は密度場のゼロモードから読み取られます。この構成により、各非ゼロ mode で

```math
\boldsymbol{k}\cdot\hat{\boldsymbol{j}}(\boldsymbol{k}) = 0
```

が満たされます。最後に、実空間での RMS 速度

```math
u_{\rm rms}
=
\sqrt{
\frac{1}{L_xL_y}
\int |\boldsymbol{u}(\boldsymbol{x})|^2\,dA
}
```

が `urms` になるように全体が正規化されます。

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `n0` | double | yes | 初期渦度 band の中心 mode index 半径 |
| `sigma` | double | yes | 初期渦度 band の mode index 幅 |
| `urms` | double | yes | 目標 RMS 速度。`u_rms`, `velocity_rms` も同義 |
| `seed` | integer | no | ランダム位相の seed。省略時は `12345` |

### `equilibrium/gaussian/compressible`

圧縮性流体向けに、運動量密度の各成分を Gaussian 平衡分布から生成します。

```sh
set density uniform value 1.0
set momentum all equilibrium/gaussian/compressible kBT 1.0 seed 12345
```

密度場のゼロ mode から平均密度 $\rho_0$ を読み取り、非ゼロ mode の運動量揺らぎを生成します。
ゼロ mode の運動量は `0` に固定されます。
`target all` では $j_x$ と $j_y$ が独立な Gaussian mode として生成されます。

計算格子点数を $N = N_x N_y$、セル面積を $\Delta A = L_x L_y / N$ とします。
ゼロ mode は

```math
\hat{j}_x(\boldsymbol{0})
=
\hat{j}_y(\boldsymbol{0})
=
0
```

に固定されます。
各独立な非ゼロ Fourier mode では、各成分 $a=x,y$ について

```math
\hat{j}_a(\boldsymbol{k})
=
\sigma_j(\xi_{a,1} + i\xi_{a,2}),
\qquad
\sigma_j^2
=
\frac{k_B T\,N\,\rho_0}{2\Delta A}
```

から生成します。
ここで $\xi_{a,1},\xi_{a,2}$ は互いに独立な標準正規乱数です。
1成分あたりの分布は

```math
P(\hat{j}_a(\boldsymbol{k}))
=
\frac{1}{2\pi\sigma_j^2}
\exp\left[
-\frac{|\hat{j}_a(\boldsymbol{k})|^2}{2\sigma_j^2}
\right]
```

です。
`target all` では $j_x$ と $j_y$ も互いに独立に生成されます。
実空間の運動量場が実数になるよう、共役な mode は Hermitian 対称性で決まります。

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `kBT` | double | yes | 揺らぎの強さ |
| `seed` | integer | no | mode ごとの Gaussian 乱数の seed。省略時は `12345` |

`kBT = 0` の場合は、指定した運動量成分がゼロになります。

### `equilibrium/gaussian/incompressible`

非圧縮性流体向けに、各非ゼロ mode で横波成分だけを持つ Gaussian 平衡分布から運動量密度場を生成します。

```sh
set density uniform value 1.0
set momentum all equilibrium/gaussian/incompressible kBT 1.0 seed 12345
```

密度場のゼロ mode から平均密度 $\rho_0$ を読み取ります。
各非ゼロ mode では、1つの transverse scalar Gaussian mode を生成し、

```math
\hat{j}_T(\boldsymbol{k})
=
\sigma_j(\xi_1+i\xi_2),
\qquad
\sigma_j^2
=
\frac{k_B T\,N\,\rho_0}{2\Delta A}
```

とします。
ここで $N = N_xN_y$、$\Delta A = L_xL_y/N$ であり、$\xi_1,\xi_2$ は互いに独立な標準正規乱数です。
transverse scalar の分布は

```math
P(\hat{j}_T(\boldsymbol{k}))
=
\frac{1}{2\pi\sigma_j^2}
\exp\left[
-\frac{|\hat{j}_T(\boldsymbol{k})|^2}{2\sigma_j^2}
\right]
```

です。
その後、

```math
\hat{j}_x(\boldsymbol{k}) = -\frac{k_y}{|\boldsymbol{k}|}\hat{j}_T(\boldsymbol{k}),
\qquad
\hat{j}_y(\boldsymbol{k}) = \frac{k_x}{|\boldsymbol{k}|}\hat{j}_T(\boldsymbol{k})
```

として $j_x$, $j_y$ に戻します。
この構成により、

```math
\boldsymbol{k}\cdot\hat{\boldsymbol{j}}(\boldsymbol{k}) = 0
```

が満たされます。
ゼロ mode の運動量は `0` に固定されます。
実空間の運動量場が実数になるよう、共役な mode は Hermitian 対称性で決まります。

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `kBT` | double | yes | 揺らぎの強さ |
| `seed` | integer | no | mode ごとの Gaussian 乱数の seed。省略時は `12345` |

非圧縮性モードでは、通常 `target all` で指定してください。
`x` または `y` だけを個別に指定すると、もう一方の成分と組み合わさらないため、初期運動量の横波条件チェックに失敗することがあります。

## 制限・注意

- `set momentum` に指定するパラメータは、運動量密度場の初期条件を生成するためのパラメータです。時間発展で使われる `model` や `fix` のパラメータとは自動的には同期しません。
- `<target>` には `x`, `y`, `0`, `1`, `all` のいずれかを指定します。
- `sine` の波数は `0 <= nkx < Nx/2`, `-Ny/2 < nky < Ny/2` の範囲で指定します。
- `sine` では `(nkx, nky) = (0, 0)` は指定できません。
- `random_vorticity` は、通常 `target all` で指定してください。
- `random_vorticity` は密度のゼロモードから $\rho_0$ を読むため、事前に非ゼロの密度初期条件を与えてください。
- `equilibrium/gaussian/compressible` と `equilibrium/gaussian/incompressible` は、運動量揺らぎの分散を決める補助量として密度のゼロ mode から $\rho_0$ を読みます。事前に正の平均密度を持つ密度初期条件を与えてください。
- `equilibrium/gaussian/compressible` と `equilibrium/gaussian/incompressible` では `kBT >= 0` である必要があります。
- `equilibrium/gaussian/incompressible` は、通常 `target all` で指定してください。
- `equilibrium/gaussian` 系のランダム成分は非ゼロ mode のみに生成されます。ゼロ mode の運動量は `0` に固定されます。
- `target all` を指定すると、同じ type と引数が $j_x$ と $j_y$ の両方に適用されます。
- 複数の momentum 初期条件を同じ成分に指定した場合、後から指定した初期条件がその成分全体を上書きします。
- 異なる成分に指定した初期条件は、互いに独立に適用されます。

## 関連コマンド

- [`set density`](./set_density.md)
- [`grid`](./grid.md)
- [`time_evolution`](./time_evolution.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
