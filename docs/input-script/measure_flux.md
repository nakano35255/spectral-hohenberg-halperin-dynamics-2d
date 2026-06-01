# `measure flux`

## 目的

密度場、スカラー場、運動量方程式で使われたフラックスの空間平均をファイルへ出力します。
`flux` は、時間積分中に評価されたフラックスを蓄積し、そのゼロ波数 mode を取り出して出力します。

## 形式

```txt
measure <ID> flux <on|off> [key value ...]
```

## 例

```sh
measure flx flux on nevery 100 file output/flux.dat
measure rho_flux flux on nevery 100 file output/rho_flux.dat target density
measure psi_flux flux on nevery 100 file output/psi_flux.dat target order_parameter
measure mom_flux flux on nevery 100 file output/momentum_flux.dat target momentum
measure flx flux off
```

## 引数

- `<ID>`
  - 型: string
  - measure の識別子です。
  - 同じ `<ID>` を再指定すると、既存の measure は置き換えられます。
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `nevery <integer>`
  - `on` のとき必須です。
  - 何ステップごとに出力するかを指定します。
- `file <filename>`
  - `on` のとき必須です。
  - 出力ファイル名を指定します。
- `target <value>`
  - 出力するフラックスの種類を指定します。
  - 指定可能な値:
    - `all`
    - `density`
    - `rho`
    - `order_parameter`
    - `psi`
    - `momentum`
    - `j`

## デフォルト値

`measure flux` を指定しない場合、フラックスは出力されません。

`target` を省略した場合は以下が使われます。

```txt
target all
```

## 出力タイミング

`measure flux` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われ、通算ステップ番号が `nevery` で割り切れるときに出力されます。

## 出力内容

出力ファイルは rank 0 が書き込みます。
各列は以下の形式です。

```txt
# step time <flux columns...>
```

`flux` measure は、各フラックスの Fourier 空間ゼロ mode を計算格子点数で割った値を出力します。
これは計算領域上の空間平均に対応します。

## `target density`

密度方程式のフラックスを出力します。

```txt
density_flux_x density_flux_y
```

## `target order_parameter`

スカラー場の各成分について、スカラー場方程式のフラックスを出力します。

```txt
order_parameter0_flux_x order_parameter0_flux_y ...
```

`order_parameters 0` の場合、`target order_parameter` または `target psi` は使用できません。

## `target momentum`

運動量方程式のフラックステンソルを出力します。

```txt
momentum_flux_xx momentum_flux_xy momentum_flux_yy
```

2次元の対称テンソル成分として、`xx`, `xy`, `yy` の3成分を出力します。

## `target all`

密度フラックス、運動量フラックス、および存在する場合はスカラー場フラックスを出力します。
`order_parameters 0` の場合、スカラー場フラックスの列は出力されません。

## 制限・注意

- `nevery` は正の整数である必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- `file` は snapshot のような prefix ではなく、出力ファイル名そのものです。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
- [`fix ... noise`](./fix_noise.md)
- [`order_parameters`](./order_parameters.md)
