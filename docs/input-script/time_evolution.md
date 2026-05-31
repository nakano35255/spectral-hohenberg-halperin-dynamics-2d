# `time_evolution`

## 目的

時間発展に使用する積分スキームと、流体場の動作モードを指定します。

## 形式

```txt
time_evolution <type>
```

## 例

```sh
time_evolution srk3/compressible
time_evolution srk3/incompressible
time_evolution srk3/quiescent
time_evolution euler/compressible
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `srk3/compressible`
    - `srk3/incompressible`
    - `srk3/quiescent`
    - `euler/compressible`
    - `euler/incompressible`
    - `euler/quiescent`

## デフォルト値

省略時は以下が使われます。

```txt
time_evolution srk3/compressible
```

ただし、入力スクリプトで明示的に指定することを推奨します。

## 詳細

`time_evolution` は、次の2つをまとめて指定します。

- 時間積分スキーム: `srk3` または `euler`
- 流体モード: `compressible`, `incompressible`, `quiescent`

現在の実装では、積分器だけを指定する独立した `integrator` コマンドはありません。
例えば、SRK3で非圧縮性モードを使う場合は、以下のように指定します。

```sh
time_evolution srk3/incompressible
```

各モードの意味は以下の通りです。

| モード | 内容 |
| --- | --- |
| `compressible` | 密度場、運動量密度場、スカラー場を時間発展します。 |
| `incompressible` | 基準密度を用い、運動量の横波成分を時間発展します。 |
| `quiescent` | 流体速度をゼロとし、スカラー場のみを時間発展します。 |

各モードで解かれる方程式は [基礎方程式](../equations.md) を参照してください。

## 制限・注意

- 指定可能な値は、上記の6種類のみです。
- 未知の `<type>` を指定するとエラーになります。
- `quiescent` モードでは、流体の移流項を有効化できません。
- `fix ... momentum nonlinear on` や `fix ... order_parameter nonlinear on` は、`quiescent` モードでは使用できません。
- 通常の本計算では `srk3/...` を推奨します。
- `euler/...` は実装確認や比較検証に向いています。

## 関連コマンド

- [`timestep`](./timestep.md)
- [`run`](./run.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
- [`fix ... noise`](./fix_noise.md)
