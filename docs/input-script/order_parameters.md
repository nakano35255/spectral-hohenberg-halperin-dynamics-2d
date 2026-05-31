# `order_parameters`

## 目的

時間発展させるスカラー場（オーダーパラメータ）の成分数を指定します。

## 形式

```txt
order_parameters <num>
```

別名として、単数形の `order_parameter` も使用できます。

```txt
order_parameter <num>
```

## 例

```sh
order_parameters 0
order_parameters 1
order_parameters 2
```

## 引数

- `<num>`
  - 型: integer
  - スカラー場の成分数
  - `0` 以上である必要があります。

## デフォルト値

省略時は以下が使われます。

```txt
order_parameters 0
```

この場合、スカラー場の方程式は解かず、密度場と運動量密度場のみを含む流体計算になります。

## 詳細

`order_parameters` は、状態変数

```math
\psi_\alpha(\boldsymbol{x},t),
\qquad
\alpha = 0,1,\dots,N_\psi-1
```

の成分数 $N_\psi$ を指定します。

`order_parameters 0` とした場合、スカラー場は存在しません。
この設定は、等温の揺らぐNavier-Stokes型ダイナミクスだけを計算したい場合に使います。

`order_parameters 1` 以上を指定した場合は、各成分に対して初期条件、自由エネルギー係数、移動度などを設定できます。

## 制限・注意

- `<num>` は `0` 以上の整数である必要があります。
- `model free_energy` や `model transport` は、指定時点の `order_parameters` の値を用いて係数配列を初期化します。
- そのため、`order_parameters` は `model free_energy` や `model transport` より前に指定してください。
- `fix ... order_parameter nonlinear on` を使う場合は、`order_parameters` が `1` 以上である必要があります。

## 関連コマンド

- [`model free_energy`](./model_free_energy.md)
- [`model transport`](./model_transport.md)
- [`set order_parameter`](./set_order_parameter.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
