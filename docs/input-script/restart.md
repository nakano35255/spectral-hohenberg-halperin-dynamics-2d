# `restart`

## 目的

`restart` は、シミュレーション状態をファイルへ保存し、別の実行でその状態から計算を再開するためのコマンドです。

restart ファイルには、ソルバー内部で使っている Fourier 空間の状態が保存されます。
保存される状態は、密度場、オーダーパラメータ場、運動量密度場を含む spectral state です。

## 形式

```txt
restart off
restart read file <filename>
restart write file <filename>
```

## 例

計算終了時に restart ファイルを書き出します。

```sh
restart write file output/restart.dat

run 100000
```

保存した状態から再開します。

```sh
restart read file output/restart.dat
restart write file output/restart_next.dat

run 100000
```

restart を使わない場合は、明示的に無効化できます。

```sh
restart off
```

## 引数

- `off`
  - restart input と restart output の両方を無効にします。
- `read file <filename>`
  - 指定した restart ファイルを読み込みます。
  - 読み込まれた state、step、time から計算を開始します。
- `write file <filename>`
  - 計算終了時に restart ファイルを書き出します。
  - 出力先ファイルが既に存在する場合は上書きされます。

## 実行タイミング

`restart read` は、シミュレーション開始時に一度だけ実行されます。
`restart read` が有効な場合、`set density`、`set momentum`、`set order_parameter` による初期条件は使いません。

`restart write` は、入力スクリプト中のすべての `run` と `measure` 処理が終わった後に一度だけ実行されます。
複数の `run` を書いた場合でも、各 `run` の終端ではなく、スクリプト全体の最後の state が保存されます。

## 初期条件との関係

`restart read` と `set` コマンドは同時に使えません。

```sh
restart read file output/restart.dat

# これはエラーになります
set density uniform value 1.0
```

restart から再開する input script では、`set density`、`set momentum`、`set order_parameter` を書かないでください。
一方、`model`、`fix`、`measure`、`run` は通常通り指定します。

## チェックされる項目

読み込み時には、restart ファイルの格子サイズと state の field 数が、現在の input script と一致しているかを確認します。

- computational grid size: `Nx`, `Ny`
- R2C spectral grid size: `nkx = Nx / 2 + 1`, `nky = Ny`
- state field count: `num_fields = order_parameters + 3`

ここで `Nx`, `Ny` は [`dealias`](./dealias.md) 適用後の computational grid size です。
active grid size や `dealias` の指定そのものは、別々にはチェックしません。

輸送係数、自由エネルギー係数、ノイズ強度、外力などの model/fix parameter はチェックしません。
restart 後にこれらのパラメータを意図的に変更して計算を続けることができます。

## MPI 並列数

restart ファイルは rank 0 に state を集約して保存されます。
読み込み時には、rank 0 がファイルを読み込み、現在の MPI 分割に合わせて各 rank へ state を配ります。

そのため、書き出し時と読み込み時の MPI 並列数は同じである必要はありません。

## ファイル形式

restart ファイルはテキスト形式です。
ヘッダの後に、global spectral layout のデータが `field kx_index ky_index real imag` の順で出力されます。

```txt
SHHD_RESTART_V1
step <step>
time <time>
nx <Nx>
ny <Ny>
nkx <Nx/2+1>
nky <Ny>
num_order_parameters <N>
num_fields <N+3>
precision text_float64
layout field_major_ky_kx_text
columns field kx_index ky_index real imag
data
<field> <kx_index> <ky_index> <real> <imag>
...
```

field index は以下に対応します。

| field | 内容 |
| --- | --- |
| `0` | density `rho` |
| `1 ... order_parameters` | order parameter `psi` |
| `order_parameters + 1` | momentum `jx` |
| `order_parameters + 2` | momentum `jy` |

`kx_index`, `ky_index` は spectral array 上の index です。
実際の波数値ではありません。

## 制限・注意

- `restart read` と `set` コマンドは同時に使えません。
- `restart write` はスクリプト全体の最後に一度だけ実行されます。
- 出力先ディレクトリは自動作成されません。
- `file` には空白を含めないでください。
- restart ファイルは text format のため、大きな格子ではファイルサイズが大きくなります。
- 書き出し時と読み込み時の MPI 並列数は変更できます。
- restart ファイルは computational grid 上の spectral state を保存します。active range 外の mode もファイルには含まれます。

## 関連コマンド

- [`run`](./run.md)
- [`grid`](./grid.md)
- [`dealias`](./dealias.md)
- [`set density`](./set_density.md)
- [`set momentum`](./set_momentum.md)
- [`set order_parameter`](./set_order_parameter.md)
