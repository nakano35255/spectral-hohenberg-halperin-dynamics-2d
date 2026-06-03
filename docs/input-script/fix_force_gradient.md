# `fix ... force/gradient`

## 目的

スカラー場方程式に、背景勾配と速度場の結合に対応する外力項を加えます。

## 形式

```txt
fix <ID> order_parameter force/gradient <on|off> component <component> direction <direction> amplitude <value>
```

## 例

```sh
fix 5 order_parameter force/gradient on component 0 direction y amplitude 1.0
fix 6 psi force/gradient on component 0 direction x amplitude 0.5

fix 5 order_parameter force/gradient off
```

## 引数

- `<ID>`
  - 型: string
  - fix の識別子です。
  - 同じ `<ID>` の `force/gradient` を再指定すると、古い設定は置き換えられます。
  - `off` を指定すると、同じ `<ID>` の `force/gradient` が削除されます。
- `<target>`
  - 型: string
  - 指定可能な値:
    - `order_parameter`
    - `psi`
  - `momentum`, `j`, `all` は指定できません。
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `component <component>`
  - 外力を加えるスカラー場成分の整数 index を指定します。
- `direction <direction>`
  - 結合させる速度成分を指定します。
  - 指定可能な値は `x`, `y`, `0`, `1` です。
- `amplitude <value>`
  - 背景勾配の大きさを指定します。

## デフォルト値

省略時は、勾配型外力はありません。

`fix ... force/gradient on` では、`component`, `direction`, `amplitude` はすべて明示的に指定する必要があります。

## 詳細

`fix ... force/gradient` は、保存型 flux としてではなく、スカラー場方程式の右辺へ直接 source として加えられる外力です。

指定したスカラー場成分 $\psi_\alpha$ に対して、以下の項を加えます。

```math
\partial_t \psi_\alpha = \cdots - G v_a
```

ここで `component` が $\alpha$、`direction` が $a$、`amplitude` が $G$ です。

例えば

```sh
fix 5 order_parameter force/gradient on component 0 direction y amplitude 1.0
```

は、スカラー場方程式に

```math
\partial_t \psi_0 = \cdots - 1.0\, v_y
```

を加えます。

この形は、平均場に

```math
\langle \psi \rangle = \psi_0 + G y
```

のような背景勾配があるとき、揺らぎ場の方程式に現れる移流項 $-G v_y$ に対応します。
`direction x` を指定した場合は $-G v_x$、`direction y` を指定した場合は $-G v_y$ が加わります。

`amplitude` は $G$ として解釈され、実際に加わる項は `- amplitude * v_direction` です。
符号を反転したい場合は、負の `amplitude` を指定してください。

非圧縮モードでは、速度場はスペクトル空間の運動量密度から

```math
\boldsymbol{v} = \boldsymbol{j}/\rho_0
```

として評価されます。
圧縮性モードでは、物理空間で

```math
\boldsymbol{v}(\boldsymbol{x}) = \boldsymbol{j}(\boldsymbol{x})/\rho(\boldsymbol{x})
```

を評価し、その後 Fourier 空間へ変換して RHS に加えます。

この外力は診断用の `FluxBuffer` には記録されません。
そのため、`measure ... flux` で測る flux には、外力そのものは含まれません。

## 制限・注意

- `fix` コマンドは、現在の実装では `run` または `measure` より前に指定する必要があります。
- `force/gradient` は `order_parameter` または `psi` target 専用です。
- `order_parameters` は `1` 以上である必要があります。
- `quiescent` モードでは、流体速度を時間発展しないため `force/gradient` は使用できません。
- 現在の `force/gradient` は保存型 flux ではなく、RHS に直接足される非保存型 source です。
- 現在の実装では、外力だけで RHS を作る用途はサポートしていません。対象 field に通常の deterministic RHS が存在する設定と組み合わせて使ってください。

## 関連コマンド

- [`time_evolution`](./time_evolution.md)
- [`order_parameters`](./order_parameters.md)
- [`model transport`](./model_transport.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
- [`fix ... noise`](./fix_noise.md)
- [`fix ... force/sine`](./fix_force_sine.md)
- [`measure ... flux`](./measure_flux.md)
