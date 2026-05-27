# `measure snapshot`

## 目的

密度場、スカラー場、運動量密度場のスナップショットをファイルへ出力します。
`snapshot` は、ソルバー内部で使っている計算格子上の場をそのまま確認するための観測量です。

## 形式

```txt
measure <ID> snapshot <on|off> [key value ...]
```

## 例

```sh
measure phys snapshot on nevery 1000 file output/physical space physical
measure spec snapshot on nevery 1000 file output/spectral space spectral
measure phys snapshot off
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
- `file <prefix>`
  - `on` のとき必須です。
  - 出力ファイル名の prefix を指定します。
- `space <value>`
  - 出力する空間を指定します。
  - 指定可能な値:
    - `physical`
    - `real`
    - `spectral`
    - `fourier`
    - `spec`

## デフォルト値

`measure snapshot` を指定しない場合、スナップショットは出力されません。

`space` を省略した場合は以下が使われます。

```txt
space physical
```

## 出力タイミング

`measure snapshot` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われ、通算ステップ番号が `nevery` で割り切れるときに出力されます。

出力ファイル名は、指定した `file` prefix にステップ番号を付けたものになります。

```txt
<prefix>_step<step>.snapshot
```

## 出力格子

`snapshot` は active grid ではなく、padding 後の computational grid に対応するデータを出力します。

ここで **active grid** は [`grid`](./grid.md) で指定する基準格子です。
一方、**computational grid** は [`dealias`](./dealias.md) の指定を反映した FFT 計算用の格子です。
`dealias none` の場合だけ、active grid と computational grid は一致します。

指定値と computational grid の対応は以下の通りです。

| `dealias` | `Nx computational` | `Ny computational` |
| --- | --- | --- |
| `none`, `off` | `Nx active` | `Ny active` |
| `three_halves`, `3/2` | `3 Nx active / 2` | `3 Ny active / 2` |
| `two`, `2` | `2 Nx active` | `2 Ny active` |

physical snapshot の出力格子は `Nx computational x Ny computational` です。
spectral snapshot の出力格子は R2C 表現に対応し、`(floor(Nx computational / 2) + 1) x Ny computational` です。

例えば、

```sh
grid 64 64
dealias two
```

の場合、physical snapshot は `128 x 128` の実空間データを出力します。
spectral snapshot は R2C 表現に対応する `65 x 128` の複素スペクトルを出力します。

padding によって増えた高波数側の mode は、時間発展させる active mode ではありません。
通常、これらの inactive mode は spectral mask によって `0` に保たれます。
そのため snapshot の配列サイズは computational grid に対応しますが、時間発展の自由度は active grid によって決まります。

## `space physical`

`space physical` では、Fourier 空間の状態を computational grid 上へ逆変換した実空間データが出力されます。
出力は rank 0 に集約されます。

ヘッダには computational grid のサイズが出力されます。

```txt
# nx <Nx computational> ny <Ny computational> order_parameters <N>
```

列は以下の順です。

```txt
x y rho psi_com0 ... jx jy
```

`x`, `y` は physical snapshot 上の格子インデックスです。
物理座標ではありません。

## `space spectral`

`space spectral` では、R2C Fourier 空間の複素スペクトルが実部・虚部に分けて出力されます。
出力は rank 0 に集約されます。

ヘッダには computational grid に対応する R2C スペクトルサイズが出力されます。

```txt
# nkx <floor(Nx computational / 2) + 1> nky <Ny computational> order_parameters <N>
```

列は以下の順です。

```txt
kx_index ky_index rho_real rho_imag psi_com0_real psi_com0_imag ... jx_real jx_imag jy_real jy_imag
```

`kx_index`, `ky_index` は spectral snapshot 上の波数インデックスです。
実際の波数値ではありません。

`space spectral` では active range の外側の inactive mode も出力配列に含まれます。
これらは padding 領域に対応する mode であり、通常は `0` です。

## フィールド列

`order_parameters` が `1` 以上の場合、`psi_com0`, `psi_com1`, ... が出力されます。
`order_parameters 0` の場合、`psi_com...` の列は出力されません。

すべての snapshot には、密度 `rho` と運動量密度 `jx`, `jy` が含まれます。

## 制限・注意

- `nevery` は正の整数である必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。例えば `file output/physical` を指定する場合、`output` ディレクトリが存在している必要があります。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。
- `dealias` を有効にしている場合、snapshot の行数は active grid ではなく computational grid に従って増えます。

## 関連コマンド

- [`run`](./run.md)
- [`grid`](./grid.md)
- [`dealias`](./dealias.md)
- [`order_parameters`](./order_parameters.md)
