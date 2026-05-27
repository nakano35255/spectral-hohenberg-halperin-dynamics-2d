# `fix ... nonlinear`

## 目的

運動量方程式およびスカラー場方程式の移流項を有効化または無効化します。

## 形式

```txt
fix <ID> all nonlinear <on|off> <target>...
```

## 例

```sh
fix 1 all nonlinear on momentum
fix 1 all nonlinear on order_parameter
fix 1 all nonlinear on all

fix 1 all nonlinear off momentum
```

## 引数

- `<ID>`
  - 型: string
  - fix の識別子です。
  - 現在の実装では、値は構文上のラベルとして扱われます。
- `all`
  - 対象グループです。
  - 現在指定できる値は `all` のみです。
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `<target>...`
  - 型: string
  - 1つ以上指定します。
  - 指定可能な値:
    - `momentum`
    - `j`
    - `order_parameter`
    - `psi`
    - `all`

## デフォルト値

省略時は、運動量方程式とスカラー場方程式の移流項はいずれも無効です。

```txt
momentum_advection off
order_parameter_advection off
```

## 詳細

`fix ... nonlinear` は、非線形な移流項を使うかどうかを指定します。

`momentum` または `j` を指定すると、運動量方程式の移流項を制御します。

```math
-\partial_b(j_a v_b)
```

`order_parameter` または `psi` を指定すると、スカラー場方程式の移流項を制御します。

```math
-\nabla\cdot(\psi_\alpha \boldsymbol{v})
```

`all` を指定すると、存在する場の移流項をまとめて制御します。
`order_parameters 0` の場合、`all` は運動量の移流項だけを制御します。
`order_parameters` が `1` 以上の場合は、運動量とスカラー場の両方の移流項を制御します。

```sh
fix 1 all nonlinear on all
```

複数の target を並べることもできます。

```sh
fix 1 all nonlinear on momentum order_parameter
```

## 制限・注意

- `fix` コマンドは、現在の実装では `run` または `measure` より前に指定する必要があります。
- `<target>` は1つ以上指定する必要があります。
- `all` 以外の group は未実装です。
- `quiescent` モードでは、流体速度をゼロとして扱うため `fix ... nonlinear` は使用できません。
- `order_parameter` または `psi` を明示的に指定する場合、`order_parameters` は `1` 以上である必要があります。
- `all` は `order_parameters 0` でも使用できます。
- 現在の実装では `<ID>` による個別管理や削除は行いません。

## 関連コマンド

- [`time_evolution`](./time_evolution.md)
- [`order_parameters`](./order_parameters.md)
- [`dealias`](./dealias.md)
- [`fix ... noise`](./fix_noise.md)
