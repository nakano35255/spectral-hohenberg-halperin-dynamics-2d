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

## 制限・注意

- `<target>` には `x`, `y`, `0`, `1`, `all` のいずれかを指定します。
- `sine` の波数は `0 <= nkx < Nx/2`, `-Ny/2 < nky < Ny/2` の範囲で指定します。
- `sine` では `(nkx, nky) = (0, 0)` は指定できません。
- `target all` を指定すると、同じ type と引数が $j_x$ と $j_y$ の両方に適用されます。
- 複数の momentum 初期条件を同じ成分に指定した場合、後から指定した初期条件がその成分全体を上書きします。
- 異なる成分に指定した初期条件は、互いに独立に適用されます。

## 関連コマンド

- [`set density`](./set_density.md)
- [`grid`](./grid.md)
- [`time_evolution`](./time_evolution.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
