# `measure ave/profile`

## 目的

指定した物理量の1次元実空間プロファイルを、空間平均および時間平均してファイルへ出力します。

`ave/profile` は、非平衡定常状態で速度場やフラックスの平均プロファイルを測るための measure です。
例えば、正弦波外力で駆動した流れに対して、`vx` と `pi_xy` の1次元プロファイルを測定し、粘性係数の推定に使うことを想定しています。

## 形式

```txt
measure <ID> ave/profile <on|off> axis <x|y> nevery <integer> nblock <integer> file <filename> average <block|running> target <target1> [target2 ...]
```

`target` は必ず最後に指定します。
`target` の後ろに並ぶ token はすべて測定対象として解釈されます。

## 例

```sh
measure prof ave/profile on axis y nevery 100 nblock 1000 file output/profile.dat average block target vx pi_xy
measure prof ave/profile on axis y nevery 100 nblock 1000 file output/profile_running.dat average running target rho psi[0] vx vy Jpsi_y[0] pi_xy

measure prof ave/profile off
```

## 引数

- `<ID>`
  - 型: string
  - measure の識別子です。
  - 同じ `<ID>` を再指定すると、既存の measure は置き換えられます。
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `axis <x|y>`
  - プロファイルの独立変数となる軸を指定します。
  - `axis x` では、`y` 方向に空間平均した $x$ の関数を出力します。
  - `axis y` では、`x` 方向に空間平均した $y$ の関数を出力します。
- `nevery <integer>`
  - 何ステップごとに観測するかを指定します。
- `nblock <integer>`
  - 何ステップごとに平均プロファイルを出力するかを指定します。
  - `nblock` は `nevery` で割り切れる必要があります。
- `file <filename>`
  - 出力ファイル名を指定します。
- `average <block|running>`
  - 時間平均の出力方法を指定します。
  - `block` は、各 block の平均プロファイルをファイル末尾へ追記します。
  - `running` は、これまでに完了した block すべての running average をファイルへ書き直します。
- `target <target1> [target2 ...]`
  - 測定する物理量を指定します。
  - `target` はコマンドの最後に置く必要があります。

## デフォルト値

`measure ave/profile` を指定しない場合、プロファイルは出力されません。

`on` の場合、`axis`, `nevery`, `nblock`, `file`, `average`, `target` はすべて明示的に指定する必要があります。

## 出力タイミング

`measure ave/profile` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われ、通算ステップ番号が `nevery` で割り切れるときに現在のプロファイルを積算します。

`nblock` ステップ分の観測が終わると、平均されたプロファイルを出力します。
例えば、

```txt
nevery 100 nblock 1000
```

では、100 step ごとに観測し、10個の観測値を平均して1つの block profile を作ります。

## 出力内容

出力ファイルは rank 0 が書き込みます。
各行は、1つの profile grid point に対応します。

```txt
# measure ave/profile
# axis <x|y> average <block|running>
# columns block step time index coord <target columns...>
```

出力列は以下の通りです。

- `block`
  - 完了した block の番号です。
- `step`
  - 出力時点の通算ステップ番号です。
- `time`
  - 出力時点の通算物理時刻です。
- `index`
  - profile 軸上の格子 index です。
- `coord`
  - profile 軸上の物理座標です。
- `<target columns...>`
  - `target` で指定した物理量の平均プロファイルです。

`average block` では、block ごとの平均プロファイルが追記されます。
`average running` では、これまでに完了した block の running average が書き直されます。

## 空間平均

`axis x` では、出力は $x$ の関数です。

```math
\overline{q}(x) = \frac{1}{L_y}\int_0^{L_y} q(x,y)\,dy
```

離散的には、同じ `x` index を持つすべての格子点について平均します。

`axis y` では、出力は $y$ の関数です。

```math
\overline{q}(y) = \frac{1}{L_x}\int_0^{L_x} q(x,y)\,dx
```

離散的には、同じ `y` index を持つすべての格子点について平均します。

## Target

以下の target を指定できます。

| target | 意味 |
| --- | --- |
| `rho` | 密度場 $\rho$ |
| `psi[N]` | スカラー場成分 $\psi_N$ |
| `jx` | 運動量密度 $j_x$ |
| `jy` | 運動量密度 $j_y$ |
| `vx` | 速度場 $v_x$ |
| `vy` | 速度場 $v_y$ |
| `Jpsi_x[N]` | スカラー場成分 $\psi_N$ のフラックス $J^\psi_{N,x}$ |
| `Jpsi_y[N]` | スカラー場成分 $\psi_N$ のフラックス $J^\psi_{N,y}$ |
| `pi_xx` | 運動量フラックス $\Pi_{xx}$ |
| `pi_xy` | 運動量フラックス $\Pi_{xy}$ |
| `pi_yy` | 運動量フラックス $\Pi_{yy}$ |

`psi[N]` の例:

```sh
measure prof ave/profile on axis y nevery 100 nblock 1000 file output/profile.dat average block target psi[0] vx pi_xy
```

`pixx`, `pixy`, `piyy` は、それぞれ `pi_xx`, `pi_xy`, `pi_yy` の alias として使用できます。

## Velocity

`vx`, `vy` は、時間発展 mode に応じて以下のように計算されます。

- compressible mode:

```math
\boldsymbol{v}(\boldsymbol{x}) = \boldsymbol{j}(\boldsymbol{x}) / \rho(\boldsymbol{x})
```

- incompressible mode:

```math
\boldsymbol{v}(\boldsymbol{x}) = \boldsymbol{j}(\boldsymbol{x}) / \rho_0
```

ここで $\rho_0$ は密度場の空間平均です。

## Flux

`Jpsi_x[N]`, `Jpsi_y[N]`, `pi_xx`, `pi_xy`, `pi_yy` は、時間積分中に評価されたフラックスを用いて計算されます。
そのため、これらの target を指定すると、内部的に対応する flux の記録が有効になります。

`fix ... force/sine` や `fix ... force/gradient` によって加えられた外力そのものは、運動量フラックスには含まれません。

## 制限・注意

- `target` は必ず最後に指定してください。
- `nevery` は正の整数である必要があります。
- `nblock` は正の整数であり、`nevery` で割り切れる必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- `file` は snapshot のような prefix ではなく、出力ファイル名そのものです。
- `average` は `block` または `running` のいずれかです。
- `order_parameters 0` の場合、`psi[N]`, `Jpsi_x[N]`, `Jpsi_y[N]` target は使用できません。
- quiescent mode では、`jx`, `jy`, `vx`, `vy`, `pi_xx`, `pi_xy`, `pi_yy` は通常の流体場プロファイルとしては使用しないでください。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`measure time_series`](./measure_time_series.md)
- [`measure snapshot`](./measure_snapshot.md)
- [`fix ... force/sine`](./fix_force_sine.md)
- [`fix ... force/gradient`](./fix_force_gradient.md)
- [`time_evolution`](./time_evolution.md)
