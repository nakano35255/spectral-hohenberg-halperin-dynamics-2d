# `fix ... force/sine`

## 目的

運動量方程式またはスカラー場方程式に、正弦波型の外力を加えます。

## 形式

```txt
fix <ID> <target> force/sine <on|off> component <component> axis <axis> nk <integer> amplitude <value>
```

## 例

```sh
fix 3 momentum force/sine on component x axis y nk 1 amplitude 1.0
fix 4 order_parameter force/sine on component 0 axis y nk 1 amplitude 1.0

fix 3 momentum force/sine off
```

## 引数

- `<ID>`
  - 型: string
  - fix の識別子です。
  - 同じ `<ID>` の `force/sine` を再指定すると、古い設定は置き換えられます。
  - `off` を指定すると、同じ `<ID>` の `force/sine` が削除されます。
- `<target>`
  - 型: string
  - 指定可能な値:
    - `momentum`
    - `j`
    - `order_parameter`
    - `psi`
  - `on` では `all` は指定できません。
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `component <component>`
  - 外力を加える成分を指定します。
  - `momentum` または `j` では `x`, `y`, `0`, `1` を指定できます。
  - `order_parameter` または `psi` では、スカラー場成分の整数 index を指定します。
- `axis <axis>`
  - 正弦波が依存する座標軸を指定します。
  - 指定可能な値は `x`, `y`, `0`, `1` です。
- `nk <integer>`
  - 波数番号を指定します。
  - `axis x` なら $k = 2\pi\, nk / L_x$、`axis y` なら $k = 2\pi\, nk / L_y$ です。
- `amplitude <value>`
  - 外力の振幅を指定します。

## デフォルト値

省略時は、正弦波型外力はありません。

`fix ... force/sine on` では、`component`, `axis`, `nk`, `amplitude` はすべて明示的に指定する必要があります。

## 詳細

`fix ... force/sine` は、外力を保存型 flux としてではなく、時間発展方程式の右辺へ直接加えます。

`momentum` または `j` を指定した場合、運動量成分に body force を加えます。

```math
\partial_t j_a = \cdots + A \sin(k r_b)
```

ここで `component` が $a$、`axis` が $b$、`amplitude` が $A$ です。

例えば

```sh
fix 3 momentum force/sine on component x axis y nk 1 amplitude 1.0
```

は、運動量方程式に

```math
\partial_t j_x = \cdots + 1.0 \sin(2\pi y / L_y)
```

を加えます。

`order_parameter` または `psi` を指定した場合、指定したスカラー場成分に source を加えます。

```math
\partial_t \psi_\alpha = \cdots + A \sin(k r_b)
```

この外力は診断用の `FluxBuffer` には記録されません。
そのため、`measure ... flux` で測る flux や stress には、外力そのものは含まれません。

## 制限・注意

- `fix` コマンドは、現在の実装では `run` または `measure` より前に指定する必要があります。
- `target all` は `on` では指定できません。複数の場に外力を加える場合は、複数行の `fix ... force/sine` を使ってください。
- `direction` 引数はありません。現在の `force/sine` は保存型 flux ではなく、RHS に直接足される外力です。
- `nk` は `0` より大きく、指定した `axis` の active grid size の半分より小さい必要があります。
- `quiescent` モードでは、運動量場を時間発展しないため `momentum` target の `force/sine` は使用できません。
- 現在の実装では、外力だけで RHS を作る用途はサポートしていません。対象 field に通常の deterministic RHS が存在する設定と組み合わせて使ってください。

## 関連コマンド

- [`time_evolution`](./time_evolution.md)
- [`grid`](./grid.md)
- [`length`](./length.md)
- [`model transport`](./model_transport.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
- [`fix ... noise`](./fix_noise.md)
- [`fix ... force/gradient`](./fix_force_gradient.md)
- [`measure ... flux`](./measure_flux.md)
