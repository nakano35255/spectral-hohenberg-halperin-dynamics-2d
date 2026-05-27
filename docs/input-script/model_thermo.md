# `model thermo`

## 目的

密度場から圧力項を決める熱力学モデルを指定します。

## 形式

```txt
model thermo <type> [key value ...]
```

## 例

```sh
model thermo linear_eos cT 10.0
model thermo linear_eos sound_speed 10.0
```

```sh
model thermo linear_eos
model thermo coeff cT 10.0
```

## 引数

- `<type>`
  - 型: string
  - 指定可能な値:
    - `linear_eos`
- `[key value ...]`
  - 熱力学モデルごとの係数を key-value 形式で指定します。

### `linear_eos`

密度場に対して線形の圧力を与えるモデルです。

```math
p = c_T^2 \rho
```

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `cT` | double | no | 等温音速 |
| `sound_speed` | double | no | `cT` と同じ意味 |

## デフォルト値

`model thermo` を省略した場合、圧力項を持たない熱力学モデルが使われます。

各係数を省略した場合、その係数は `0.0` になります。

## coeff コマンド

`model thermo coeff` を使用すると、すでに指定済みの熱力学モデルに対して係数だけを追加・更新できます。
指定した key の値だけが上書きされ、指定しなかった係数は現在の値を保持します。

```sh
model thermo linear_eos
model thermo coeff cT 10.0
```

この例では、最初の行で `linear_eos` モデルを指定し、次の行で同じ `linear_eos` モデルに `cT` を追加設定します。
詳細は [`model ... coeff`](./model_coeff.md) を参照してください。

## 詳細

`model thermo` は、流体モードの運動量方程式で使う圧力 $p$ を定義します。
実際の時間発展では、圧力勾配として

```math
- \nabla_a p
```

に対応する項が評価されます。
`linear_eos` では、Fourier 空間で

```math
\widehat{p}(\boldsymbol{k})
= c_T^2 \widehat{\rho}(\boldsymbol{k})
```

を用います。

`incompressible` モードでは、運動量方程式の右辺を横波成分へ射影するため、圧力項は丸々更新から取り除かれます。
そのため、圧力項を `model thermo` で明示的に指定する必要はありません。

## 制限・注意

- 現在指定できる熱力学モデルは `linear_eos` のみです。
- `cT` は `0` 以上である必要があります。
- 非線形の密度依存圧力は、現在の実装では未実装です。

## 関連コマンド

- [`time_evolution`](./time_evolution.md)
- [`model ... coeff`](./model_coeff.md)
- [`model transport`](./model_transport.md)
