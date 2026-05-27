# `measure snapshot`

## 目的

密度場、スカラー場、運動量密度場のスナップショットをファイルへ出力します。

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

## 詳細

`measure snapshot` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われ、通算ステップ番号が `nevery` で割り切れるときに出力されます。

出力ファイル名は、指定した `file` prefix にステップ番号を付けたものになります。

```txt
<prefix>_step<step>.snapshot
```

`space physical` では、rank 0 に集約された実空間データが出力されます。
列は以下の順です。

```txt
x y rho psi_com0 ... jx jy
```

`space spectral` では、Fourier 空間の複素スペクトルが実部・虚部に分けて出力されます。
列は以下の順です。

```txt
kx_index ky_index rho_real rho_imag psi_com0_real psi_com0_imag ... jx_real jx_imag jy_real jy_imag
```

`order_parameters 0` の場合、`psi_com...` の列は出力されません。

## 制限・注意

- `nevery` は正の整数である必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。例えば `file output/physical` を指定する場合、`output` ディレクトリが存在している必要があります。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`grid`](./grid.md)
- [`dealias`](./dealias.md)
- [`order_parameters`](./order_parameters.md)
