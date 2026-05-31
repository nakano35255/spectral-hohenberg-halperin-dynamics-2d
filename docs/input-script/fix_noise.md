# `fix ... noise`

## 目的

保存型の熱ノイズを有効化または無効化します。

## 形式

```txt
fix <ID> <target> noise <on|off> [seed <integer>] [kBT <value>] [chi <value>]
```

## 例

```sh
fix 2 all noise on seed 12345 kBT 1.0
fix 2 momentum noise on temperature 1.0
fix 2 order_parameter noise on seed 12345 kBT 1.0 chi 2.0
fix 2 all noise off
```

## 引数

- `<ID>`
  - 型: string
  - fix の識別子です。
  - 現在の実装では、値は構文上のラベルとして扱われます。
- `<target>`
  - 型: string
  - 指定可能な値:
    - `momentum`
    - `j`
    - `order_parameter`
    - `psi`
    - `all`
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `seed <integer>`
  - 乱数シードを指定します。
- `kBT <value>`
  - ノイズ強度を決める温度を指定します。
  - `temperature <value>` も同じ意味で使えます。
- `chi <value>`
  - スカラー場ノイズだけに掛ける倍率です。
  - `order_parameter_chi <value>`、`psi_chi <value>` も同じ意味で使えます。

## デフォルト値

省略時はノイズは無効です。

明示的に無効化する場合は、例えば以下のように指定します。

```sh
fix 2 all noise off
```

`fix ... noise on` を指定した場合、明示的に与えなかった値には以下が使われます。

```txt
seed 12345
kBT 1.0
chi 1.0
```

## 詳細

`fix ... noise` は、確率フラックスおよび確率応力に対応する保存型ノイズを有効化します。

`momentum` または `j` を指定すると、運動量方程式の確率応力を制御します。

`order_parameter` または `psi` を指定すると、スカラー場方程式の保存型ノイズを制御します。

`all` を指定すると、存在する場のノイズをまとめて制御します。
`quiescent` モードでは流体場を時間発展しないため、運動量ノイズは使われません。

`chi` はスカラー場ノイズの強度だけを変更します。スカラー場の確率フラックスでは、実効的に `kBT` が `chi * kBT` として使われます。運動量ノイズには影響しません。

ノイズの共分散や方程式上の位置は [基礎方程式](../equations.md) を参照してください。

## 制限・注意

- `fix` コマンドは、現在の実装では `run` または `measure` より前に指定する必要があります。
- `kBT` は `0` 以上である必要があります。
- `chi` は `0` 以上である必要があります。
- `chi` は `momentum` target では指定できません。
- `fix ... noise off` には追加引数を指定できません。
- 現在の実装では `<ID>` による個別管理や削除は行いません。

## 関連コマンド

- [`model transport`](./model_transport.md)
- [`time_evolution`](./time_evolution.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
