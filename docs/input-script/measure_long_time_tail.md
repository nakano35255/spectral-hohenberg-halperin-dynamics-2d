# `measure long_time_tail`

## 目的

`long_time_tail` は、流体変数やオーダーパラメータの空間平均された時間相関関数を測定します。
測定する量は

```math
C_{AB}(\tau)
=
\frac{1}{V}
\int d\boldsymbol{x}\,
\langle \delta A(\boldsymbol{x}, t+\tau)\,
\delta B(\boldsymbol{x}, t) \rangle
```

です。
ここで、ゼロ Fourier mode を除外することで平均成分を常に差し引きます。
そのため、密度場やオーダーパラメータでは connected correlation が出力されます。

## 形式

```txt
measure <ID> long_time_tail <on|off> nevery <integer> nblock <integer> file <filename> [average <block|running>] [cross <on|off>] target <target1> [target2 ...]
```

`target` は必ず最後に指定してください。
`target` 以降の語はすべて測定対象として解釈されます。

## 例

```sh
measure ltt long_time_tail on nevery 10 nblock 10000 file results/long_time_tail.dat average running cross off target jx jy rho
measure ltt_cross long_time_tail on nevery 10 nblock 10000 file results/long_time_tail_cross.dat average running cross on target rho jx jy psi[0]
measure ltt long_time_tail off
```

## 引数

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `nevery` | integer | yes | 何ステップごとに場を保存するか |
| `nblock` | integer | yes | 何ステップを1つの block とするか |
| `file` | string | yes | 出力ファイル名 |
| `average` | string | no | 出力する平均の種類。`block` または `running`。省略時は `running` |
| `cross` | string | no | cross correlation を出力するか。`on` または `off`。省略時は `off` |
| `target` | list | yes | 相関関数を測定する場の一覧 |

`nevery` と `nblock` は正の整数である必要があります。
`nblock` は `nevery` で割り切れる必要があります。
出力先ディレクトリは自動作成されません。

## Block と Lag

`nblock` は simulation step 数で測った block 長です。
1 block 内に保存されるサンプル数、および出力される lag 点数は

```math
n_{\mathrm{lag}} = \frac{\mathrm{nblock}}{\mathrm{nevery}}
```

です。
出力される遅れ時間は

```math
\tau_m = m\,\mathrm{nevery}\,\Delta t,
\qquad
m = 0, 1, \ldots, n_{\mathrm{lag}} - 1
```

です。
したがって最大の遅れ時間は `(nblock - nevery) * dt` になります。

各 lag では、同じ block 内のすべての可能な time origin を平均します。
lag index が `m` のとき、1 block あたり `nlag - m` 個の time origin が使われます。
この数が `average block` で出力される `nsamples` です。
`average running` では、完了した block 数を `nblock_done` として

```math
\mathrm{nsamples}
=
n_{\mathrm{block\_done}}\,
(n_{\mathrm{lag}} - m)
```

が出力されます。

## 出力タイミング

`measure long_time_tail` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われます。
`nevery` と `nblock` は、通算ステップ番号ではなく、その measure が有効化されてからの内部カウンタで判定されます。
同じ measure が有効であり続ける間は、複数の `run` にまたがって同じ block が継続します。
同じ `<ID>` の measure を `off` または再度 `on` にすると、その measure の block は終了または作り直されます。

`nblock` ステップ分の観測が終わると、その block 内に保存された場から時間相関関数を計算して出力します。

## 出力形式

出力ファイルは空白区切りのテキストです。

```txt
# measure long_time_tail
# nevery <nevery> nblock <nblock> average <block|running> cross <on|off>
# columns nsamples tau <correlation columns...>
```

各データ行には、平均に使われたサンプル数 `nsamples`、遅れ時間 `tau`、指定された相関関数の値が出力されます。
`nsamples` は lag に依存します。
`average block` では、その block 内で使われた time origin 数です。
`average running` では、これまでに完了した block すべてを含めた time origin 数です。

`average block` では、block ごとの相関関数がファイル末尾へ追記されます。
各 block の前には

```txt
# block <block>
```

が出力されます。

`average running` では、これまでに完了した block すべての running average がファイルへ書き直されます。

## Target

以下の target を指定できます。

| target | 意味 |
| --- | --- |
| `rho` | 密度場 $\rho$ |
| `jx` | 運動量密度 $j_x$ |
| `jy` | 運動量密度 $j_y$ |
| `psi[N]` | スカラー場成分 $\psi_N$ |

`psi[N]` では、`0 <= N < order_parameters` である必要があります。
`vx`, `vy` は `long_time_tail` の target ではありません。
速度相関ではなく運動量密度相関を測る場合は `jx`, `jy` を使います。

## Correlation Columns

`cross off` では、各 target の auto-correlation だけを出力します。

```sh
measure ltt long_time_tail on nevery 10 nblock 1000 file output/ltt.dat cross off target rho jx jy
```

この場合、列は次のようになります。

```txt
# columns nsamples tau rhorho jxjx jyjy
```

`cross on` では、target の入力順に upper triangular な組を出力します。
target が `rho jx jy psi[0]` の場合、列は次のようになります。

```txt
# columns nsamples tau rhorho rhojx rhojy rhopsi[0] jxjx jxjy jxpsi[0] jyjy jypsi[0] psi[0]psi[0]
```

列名 `AB` は

```math
C_{AB}(\tau)
=
\frac{1}{V}
\int d\boldsymbol{x}\,
\langle \delta A(\boldsymbol{x}, t+\tau)\,
\delta B(\boldsymbol{x}, t) \rangle
```

に対応します。
現在の実装では、平衡での対称性を仮定して `C_AB` と `C_BA` の両方は出力せず、target リストの upper triangular な組だけを出力します。
特定の向きの cross correlation が必要な場合は、later 側にしたい target を先に書いてください。

## Fourier 空間での評価

実装では、各観測時刻に target の Fourier 成分を保存し、block の終わりに相関関数を計算します。
離散的には、各 time origin と lag について

```math
C_{AB}(\tau)
=
\frac{1}{N^2}
\sum_{\boldsymbol{k}\ne 0}
w_{\boldsymbol{k}}\,
\mathrm{Re}
\left[
A_{\boldsymbol{k}}(t+\tau)
B_{\boldsymbol{k}}(t)^*
\right]
```

を評価します。
ここで `N = Nx * Ny`、`w_k` は real-to-complex FFT 表現で省略された共役 mode を補う重みです。
`gx == 0` の mode では `w_k = 1`、`gx > 0` の mode では `w_k = 2` です。

ゼロ mode `(kx, ky) = (0, 0)` は常に除外されます。
これは、各場の空間平均を差し引いた相関関数を測ることに対応します。

## 制限・注意

- `target` は必ず最後に指定してください。
- `nevery` は正の整数である必要があります。
- `nblock` は正の整数であり、`nevery` で割り切れる必要があります。
- `nblock` は保存する lag 点数ではなく、simulation step 数での block 長です。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- `average` は `block` または `running` のいずれかです。
- `cross` は `on` または `off` のいずれかです。
- `order_parameters 0` の場合、`psi[N]` target は使用できません。
- quiescent mode では流体場を時間発展しないため、通常は `psi[N]` の相関関数を測定します。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`measure time_series`](./measure_time_series.md)
- [`measure ave/profile`](./measure_ave_profile.md)
- [`measure yokota_green_kubo`](./measure_yokota_green_kubo.md)
- [`set momentum`](./set_momentum.md)
- [`set order_parameter`](./set_order_parameter.md)
- [`time_evolution`](./time_evolution.md)
