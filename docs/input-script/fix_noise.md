# `fix ... noise`

## 目的

保存型の熱ノイズを有効化または無効化します。

## 形式

```txt
fix <ID> all noise <on|off> [seed <integer>] [kBT <value>]
```

## 例

```sh
fix 2 all noise on seed 12345 kBT 1.0
fix 2 all noise on temperature 1.0
fix 2 all noise off
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
- `seed <integer>`
  - 乱数シードを指定します。
- `kBT <value>`
  - ノイズ強度を決める温度を指定します。
  - `temperature <value>` も同じ意味で使えます。

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
```

## 詳細

`fix ... noise` は、確率フラックスおよび確率応力に対応する保存型ノイズを有効化します。

スカラー場が存在する場合、スカラー場方程式には移動度に対応した保存型ノイズが加わります。

流体を時間発展する `compressible` / `incompressible` モードでは、運動量方程式にも粘性係数に対応した確率応力が加わります。
`quiescent` モードでは流体場を時間発展しないため、運動量ノイズは使われません。

ノイズの共分散や方程式上の位置は [基礎方程式](../equations.md) を参照してください。

## 制限・注意

- `fix` コマンドは、現在の実装では `run` または `measure` より前に指定する必要があります。
- `all` 以外の group は未実装です。
- `kBT` は `0` 以上である必要があります。
- `fix ... noise off` には追加引数を指定できません。
- 現在の実装では `<ID>` による個別管理や削除は行いません。

## 関連コマンド

- [`model transport`](./model_transport.md)
- [`time_evolution`](./time_evolution.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
