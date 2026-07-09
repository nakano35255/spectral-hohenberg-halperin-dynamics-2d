# `measure correlation/static`

## 目的

`correlation/static` は、密度場、運動量密度、オーダーパラメータの等時刻 Fourier 相関を測定します。
測定する量は、非ゼロ Fourier mode に対する

```math
S_{AB}(\boldsymbol{k})
=
\left\langle
  A(\boldsymbol{k}, t)
  B(\boldsymbol{k}, t)^\ast
\right\rangle
```

です。
出力は内部 FFT の規格化に合わせて、格子点数 $N_g = N_x N_y$ により

```math
\frac{1}{N_g^2}
\left\langle
  A(\boldsymbol{k}, t)
  B(\boldsymbol{k}, t)^\ast
\right\rangle
```

として正規化されます。
ゼロ mode は常に除外されます。
これは、空間平均成分を含まない静的相関を測ることに対応します。

`mode 2d` では非等方性を保ったまま各 Fourier mode の値を出力します。
`mode shell` では同じ $|\boldsymbol{k}|$ shell に入る mode を平均し、等方化された静的相関を出力します。

## 形式

```txt
measure <ID> correlation/static <on|off> nevery <integer> nblock <integer> file <filename> mode <2d|shell> [average <block|running>] [cross <on|off>] target <target1> [target2 ...]
```

`target` は必ず最後に指定してください。
`target` 以降の語はすべて測定対象として解釈されます。

## 例

```sh
measure sc_shell correlation/static on nevery 100 nblock 10000 file results/static_corr_shell.dat mode shell average running cross off target rho jx jy psi[0]
measure sc_2d correlation/static on nevery 100 nblock 10000 file results/static_corr_2d.dat mode 2d average block cross on target rho psi[0]
measure sc_shell correlation/static off
```

定常状態に入ってからだけ測定したい場合は、先に測定なしで `run` し、その後に `measure correlation/static` を有効化します。

```sh
run                 50000

measure             sc_shell correlation/static on nevery 100 nblock 10000 file results/static_corr.dat mode shell average running cross off target rho jx jy psi[0]

run                 50000
```

## 引数

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `nevery` | integer | yes | 何ステップごとに観測するか |
| `nblock` | integer | yes | 何ステップごとに平均相関を出力するか |
| `file` | string | yes | 出力ファイル名 |
| `mode` | string | yes | 出力形式。`2d` または `shell` |
| `average` | string | no | 出力する平均の種類。`block` または `running`。省略時は `running` |
| `cross` | string | no | cross correlation を出力するか。`on` または `off`。省略時は `off` |
| `target` | list | yes | 静的相関を測定する場の一覧 |

`nevery` と `nblock` は正の整数である必要があります。
`nblock` は `nevery` で割り切れる必要があります。
出力先ディレクトリは自動作成されません。

## 出力タイミング

`measure correlation/static` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われます。
`nevery` と `nblock` は、通算ステップ番号ではなく、その measure が有効化されてからの内部カウンタで判定されます。
同じ measure が有効であり続ける間は、複数の `run` にまたがって同じ block が継続します。
同じ `<ID>` の measure を `off` または再度 `on` にすると、その measure の block は終了または作り直されます。

`nblock` ステップ分の観測が終わると、平均された静的相関を出力します。
例えば、

```txt
nevery 20 nblock 200
```

では、20 step ごとに観測し、10個の観測値を平均して1つの block correlation を作ります。

## 平均方法

`average block` では、各 block の平均相関がファイル末尾へ追記されます。
各 block の `samples` は、その block に含まれる観測サンプル数です。

`average running` では、これまでに完了した block すべての running average がファイルへ書き直されます。
この場合、最終ファイルには最後に出力された running average だけが残ります。
各 block の `samples` は、これまでに完了した block すべてに含まれる累積サンプル数です。

## 出力形式

出力ファイルは rank 0 が書き込みます。
最初に共通 header が出力されます。

```txt
# measure correlation/static
# nevery <nevery> nblock <nblock> mode <2d|shell> average <block|running> cross <on|off> normalization fft_grid
```

各 block は、次のコメント行から始まります。

```txt
# block <block> step <step> time <time> samples <samples>
```

`samples` は、その出力に含まれる時間サンプル数です。
`average block` では、その block に含まれるサンプル数です。
`average running` では、これまでに完了した block すべてに含まれる累積サンプル数です。

### `mode 2d`

```txt
# kx ky <correlation columns...>
```

各データ行は、active spectral mask 内の非ゼロ Fourier mode に対応します。
R2C 表現のため、出力される mode は $k_x \ge 0$ 側の表現です。
`mode 2d` では、共役 mode を含めるための重みは掛けません。

`cross off` では、各 target の auto-correlation だけを出力します。
例えば `target rho jx psi[0]` の場合、列は次のようになります。

```txt
# kx ky rho_rho jx_jx psi[0]_psi[0]
```

`cross on` では、target の入力順に upper triangular な組を出力します。
`mode 2d` では cross correlation は一般に複素数なので、実部と虚部を分けて出力します。
例えば `target rho jx psi[0]` の場合、列は次のようになります。

```txt
# kx ky rho_rho_re rho_rho_im rho_jx_re rho_jx_im rho_psi[0]_re rho_psi[0]_im jx_jx_re jx_jx_im jx_psi[0]_re jx_psi[0]_im psi[0]_psi[0]_re psi[0]_psi[0]_im
```

### `mode shell`

```txt
# k count <correlation columns...>
```

各データ行は、1つの shell に対応します。

| column | 意味 |
| --- | --- |
| `k` | shell の代表波数 |
| `count` | R2C 重みを含めた shell 内 mode 数 |
| correlation columns | shell 内 mode 平均された相関 |

shell index は

```math
n = \left\lfloor \frac{|\boldsymbol{k}|}{\Delta k} + \frac{1}{2} \right\rfloor,
\qquad
\Delta k = \min\left(\frac{2\pi}{L_x}, \frac{2\pi}{L_y}\right)
```

で決まり、出力される `k` は $n\Delta k$ です。

R2C 表現で省略された共役 mode を補うため、shell 平均では重み

```math
w(\boldsymbol{k})
=
\begin{cases}
1, & g_x = 0,\\
2, & g_x > 0
\end{cases}
```

を用います。
`count` はこの重みを足し合わせた値です。

`mode shell` の相関列は shell sum ではなく shell average です。
離散的には

```math
S_{AB}(n)
=
\frac{1}{N_g^2}
\frac{
  \sum_{\boldsymbol{k}\in n}
  w(\boldsymbol{k})\,
  \operatorname{Re}
  \left[
    A(\boldsymbol{k})B(\boldsymbol{k})^\ast
  \right]
}{
  \sum_{\boldsymbol{k}\in n} w(\boldsymbol{k})
}
```

を時間平均して出力します。
`cross off` の auto-correlation では、これは $\lvert A(\boldsymbol{k}) \rvert^2$ の shell average です。

`cross on` でも、`mode shell` では実部だけを出力します。
虚部を含めて mode ごとの cross correlation を見たい場合は、`mode 2d` を使ってください。

## Target

以下の target を指定できます。

| target | 意味 |
| --- | --- |
| `rho` | 密度場 $\rho$ |
| `jx` | 運動量密度 $j_x$ |
| `jy` | 運動量密度 $j_y$ |
| `psi[N]` | スカラー場成分 $\psi_N$ |

`psi[N]` では、`0 <= N < order_parameters` である必要があります。
`vx`, `vy` は `correlation/static` の target ではありません。
速度相関ではなく運動量密度相関を測る場合は `jx`, `jy` を使います。
compressible mode では $\boldsymbol{v} = \boldsymbol{j}/\rho$ なので、`jx`, `jy` は速度ではなく運動量密度です。

## 正規化と理論値との比較

出力 header の

```txt
normalization fft_grid
```

は、内部 FFT 係数を $N_g^2$ で割った相関であることを表します。
`A(k)` がこのコードの未規格化 FFT 係数であるとき、`mode 2d` の auto-correlation は

```math
\frac{\langle |A(k)|^2 \rangle}{N_g^2}
```

です。

平衡揺らぎとの比較では、入力された場の初期化やノイズの規格化と同じ Fourier 規格化を用いてください。
例えば、格子間隔を 1 として体積 $V = N_x N_y$ とみなす平衡計算では、線形化された compressible fluctuating hydrodynamics の運動量密度相関は

```math
S_{j_xj_x} = S_{j_yj_y} = \frac{k_B T \rho_0}{V}
```

と比較できます。
格子間隔や物理体積を変える場合は、使用している離散化と初期条件の規格化に合わせて解釈してください。

## 制限・注意

- `target` は必ず最後に指定してください。
- `nevery` は正の整数である必要があります。
- `nblock` は正の整数であり、`nevery` で割り切れる必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- `file` は snapshot のような prefix ではなく、出力ファイル名そのものです。
- `average` は `block` または `running` のいずれかです。
- `mode` は `2d` または `shell` のいずれかです。
- `cross` は `on` または `off` のいずれかです。
- `order_parameters 0` の場合、`psi[N]` target は使用できません。
- ゼロ Fourier mode は常に除外されます。
- `mode shell` の cross correlation は実部だけを出力します。
- compressible nonlinear simulation では、擬スペクトル離散化と $\boldsymbol{v}=\boldsymbol{j}/\rho$ の非線形性により、連続理論の揺動散逸関係が有限解像度で厳密には保たれない場合があります。このずれは、密度揺らぎが小さい準非圧縮領域で弱くなります。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`measure long_time_tail`](./measure_long_time_tail.md)
- [`measure time_series`](./measure_time_series.md)
- [`measure budget/spectrum`](./measure_budget_spectrum.md)
- [`set density`](./set_density.md)
- [`set momentum`](./set_momentum.md)
- [`set order_parameter`](./set_order_parameter.md)
- [`time_evolution`](./time_evolution.md)
