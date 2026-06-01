# `measure energetics`

## 目的

系のエネルギー診断量をファイルへ出力します。
`energetics` は、運動エネルギー、自由エネルギー、圧縮性エネルギー、およびそれらの和を計算します。

## 形式

```txt
measure <ID> energetics <on|off> [key value ...]
```

## 例

```sh
measure eng energetics on nevery 100 file output/energetics.dat
measure eng energetics off
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

## デフォルト値

`measure energetics` を指定しない場合、エネルギー診断量は出力されません。

`on` の場合、`nevery` と `file` は必須です。

## 出力タイミング

`measure energetics` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われ、通算ステップ番号が `nevery` で割り切れるときに出力されます。

## 出力内容

出力ファイルは rank 0 が書き込みます。
各列は以下の形式です。

```txt
# step time total_energy kinetic_energy free_energy compressibility_energy
```

- `step`
  - 通算ステップ番号です。
- `time`
  - 通算の物理時刻です。
- `total_energy`
  - `kinetic_energy + free_energy + compressibility_energy` です。
- `kinetic_energy`
  - 運動エネルギーです。
- `free_energy`
  - スカラー場の自由エネルギーです。
- `compressibility_energy`
  - 密度揺らぎに由来する圧縮性エネルギーです。

## 運動エネルギー

運動エネルギーの計算方法は [`time_evolution`](./time_evolution.md) の mode に依存します。

- compressible mode では、実空間で `0.5 * (jx^2 + jy^2) / rho` を積分します。
- incompressible mode では、密度のゼロ波数 mode から参照密度を求め、Fourier 空間の運動量密度から計算します。
- quiescent mode では、運動量場を時間発展しないため `0` になります。

## 自由エネルギー

自由エネルギーは、[`model free_energy`](./model_free_energy.md) で指定したモデルに基づいて計算されます。

Fourier 空間で表せる二次形式の寄与は、各スカラー場成分について `k0`, `k2`, `k4` 係数から計算されます。
例えば `quadratic` モデルでは `k0` の寄与、`phi4` モデルでは `a` と `b` の寄与がここに含まれます。

`phi4` モデルで `u` が非ゼロの場合は、実空間で `0.25 * u * psi^4` の寄与も加算されます。

`order_parameters 0` の場合、自由エネルギーは `0` になります。

## 圧縮性エネルギー

圧縮性エネルギーは compressible mode のときだけ計算されます。
それ以外の mode では `0` になります。

現在の `linear_eos` では、密度の非ゼロ波数 mode から二次の圧縮性エネルギーを計算します。
係数には [`model thermo`](./model_thermo.md) で指定した `cT` から得られる `cT^2` が使われます。

## 制限・注意

- `nevery` は正の整数である必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- `file` は snapshot のような prefix ではなく、出力ファイル名そのものです。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`time_evolution`](./time_evolution.md)
- [`model thermo`](./model_thermo.md)
- [`model free_energy`](./model_free_energy.md)
- [`order_parameters`](./order_parameters.md)
