# `boundary`

## 目的

計算領域の境界条件を指定します。

現在の実装は、2次元周期境界条件のみをサポートします。

## 形式

```txt
boundary <x> <y>
```

## 例

```sh
boundary p p
boundary periodic periodic
```

## 引数

- `<x>`
  - 型: string
  - $x$ 方向の境界条件
  - 指定可能な値:
    - `p`
    - `periodic`
- `<y>`
  - 型: string
  - $y$ 方向の境界条件
  - 指定可能な値:
    - `p`
    - `periodic`

## デフォルト値

省略した場合も、ソルバーは周期境界条件を前提として動作します。
ただし、入力スクリプトの可読性のため、明示的に `boundary p p` を指定することを推奨します。

## 詳細

`boundary` は、計算領域

```math
\Omega = [0,L_x) \times [0,L_y)
```

の境界条件を指定します。
Fourierスペクトル法を用いるため、現在の実装では両方向とも周期境界条件である必要があります。

## 制限・注意

- 非周期境界条件は未実装です。
- `<x>` と `<y>` の両方に `p` または `periodic` を指定してください。
- `boundary p` のように1方向だけを指定する形式は使えません。

## 関連コマンド

- [`dimension`](./dimension.md)
- [`grid`](./grid.md)
- [`length`](./length.md)
