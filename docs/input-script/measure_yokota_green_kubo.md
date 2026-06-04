# `measure yokota_green_kubo`

## 目的

`yokota_green_kubo` は、Yokota の Green-Kubo 型評価に使う運動量フラックスの時間積分量を Fourier mode ごとに出力します。

現在の実装では、運動量フラックスの `pi_xy` 成分を使い、選択された波数 mode ごとに

```math
S(\boldsymbol{k}, \tau)
=
\frac{1}{2 k_B T L_x L_y \tau}
\left|
\int_0^\tau \Pi_{xy}(\boldsymbol{k}, t)\,dt\,\Delta A
\right|^2
```

を block ごとに評価し、その running mean を出力します。

## 形式

```txt
measure <ID> yokota_green_kubo <on|off> nevery <integer> nblock <integer> file <filename> [mode <diagonal|all>]
```

## 例

```sh
fix     1 momentum noise on seed 12345 kBT 1.0
measure gk yokota_green_kubo on nevery 100 nblock 10000 file examples/05_yokota_green_kubo/results/yokota_green_kubo.dat mode diagonal
measure gk yokota_green_kubo off
```

## 引数

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `nevery` | integer | yes | 何ステップごとに積分値を評価するか |
| `nblock` | integer | yes | 何ステップを1つの block とするか |
| `file` | string | yes | 出力ファイル名 |
| `mode` | string | no | 出力する波数 mode の選び方。`diagonal` または `all`。省略時は `diagonal` |

`nevery` と `nblock` は正の整数である必要があります。
`nblock` は `nevery` で割り切れる必要があります。
出力先ディレクトリは自動作成されません。

## 出力タイミング

`measure yokota_green_kubo` は、スクリプト内で有効化された後の `run` に対して作用します。
各時間ステップで `pi_xy` の Fourier 成分を時間積分し、`nevery` ステップごとに各遅れ時間 $\tau$ の block 値を評価します。

`nblock` ステップが終わると `nsample` が1つ増え、これまでの block による running mean を `file` に書き出します。
出力ファイルは block ごとに追記されるのではなく、最新の running mean で上書きされます。

## 出力形式

出力ファイルは空白区切りのテキストです。

```txt
# nsample tau mode_index kx ky S_mean
```

各データ行の列は以下の通りです。

| 列 | 説明 |
| --- | --- |
| `nsample` | 平均に使われた block 数 |
| `tau` | 積分時間 |
| `mode_index` | measure 内部の mode index |
| `kx` | 物理波数の x 成分 |
| `ky` | 物理波数の y 成分 |
| `S_mean` | `S(k,tau)` の block running mean |

## Mode Selection

active grid を `Nx active`, `Ny active` としたとき、

```txt
max_nx = Nx active / 2 - 1
max_ny = Ny active / 2 - 1
```

までの mode が候補になります。
Nyquist 線は含まれません。

### `mode diagonal`

`mode diagonal` では、対角 mode

```math
(n_x, n_y) = (1,1), (2,2), \ldots, (n_{\max}, n_{\max})
```

を出力します。
ここで

```math
n_{\max} = \min(\mathrm{max\_nx}, \mathrm{max\_ny})
```

です。

### `mode all`

`mode all` では、まず y 軸上の

```math
(n_x,n_y) = (0,1), (0,2), \ldots, (0,\mathrm{max\_ny})
```

を出力し、その後

```math
n_x = 1,\ldots,\mathrm{max\_nx},
\qquad
n_y = -\mathrm{max\_ny},\ldots,\mathrm{max\_ny}
```

の全 mode を出力します。

物理波数は

```math
k_x = \frac{2\pi n_x}{L_x},
\qquad
k_y = \frac{2\pi n_y}{L_y}
```

として出力されます。

## 制限・注意

- `fix ... noise` の `kBT` が正である必要があります。
- 現在の実装では、`pi_xy` のみを使います。
- この measure を有効にすると、内部的に momentum flux の記録が有効になります。
- `mode` は `diagonal` または `all` のいずれかです。
- `nblock` は `nevery` で割り切れる必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。
- quiescent mode では流体場を時間発展しないため、通常は使用しません。

## 関連コマンド

- [`run`](./run.md)
- [`measure time_series`](./measure_time_series.md)
- [`fix ... noise`](./fix_noise.md)
- [`set momentum`](./set_momentum.md)
- [`time_evolution`](./time_evolution.md)
