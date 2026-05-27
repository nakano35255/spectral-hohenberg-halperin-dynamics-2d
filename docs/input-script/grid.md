# `grid`

## 目的

2次元計算領域の格子点数を指定します。

## 形式

```txt
grid <Nx> <Ny>
```

## 例

```sh
grid 64 64
grid 128 256
```

## 引数

- `<Nx>`
  - 型: integer
  - $x$ 方向の active grid 数
  - 正の偶数である必要があります。
- `<Ny>`
  - 型: integer
  - $y$ 方向の active grid 数
  - 正の偶数である必要があります。

## デフォルト値

省略時は以下が使われます。

```txt
grid 128 128
```

ただし、入力スクリプトで明示的に指定することを推奨します。

## 詳細

`grid` は、active grid 数を指定します。

ここで **active grid** とは、時間発展させるスペクトル mode の範囲を決める基準格子です。
言い換えると、シミュレーションの解像度としてユーザーが直接指定する格子数です。

一方、`dealias` を有効にした場合、非線形項評価のために実際のFFT計算格子は active grid より大きい padding されたサイズになります。
この padding 後の格子は、aliasing を抑えるための作業格子であり、時間発展させる自由度そのものを増やすものではありません。

例えば、

```sh
grid 64 64
dealias 3/2
```

と指定した場合、active grid は $64 \times 64$ ですが、FFT計算に使われる格子は $96 \times 96$ になります。

スペクトル空間では、active grid に対応する波数だけが時間発展に使われます。
padding により増えた高波数側の mode は、aliasing を抑えるための作業領域として扱われます。

## 制限・注意

- `<Nx>` と `<Ny>` は正の偶数である必要があります。
- FFTの効率のため、2の累乗を指定することを推奨します。
- `dealias` と組み合わせる場合、出力される snapshot の格子数は padding 後の計算格子に対応します。
- MPI の x-slab 分割の都合上、MPIプロセス数は格子サイズに対して大きすぎない必要があります。

## 関連コマンド

- [`length`](./length.md)
- [`dealias`](./dealias.md)
- [`dimension`](./dimension.md)
