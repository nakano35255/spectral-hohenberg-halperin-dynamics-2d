# `dealias`

## 目的

擬スペクトル法で非線形項を評価するときの aliasing を抑えるため、padding ルールを指定します。

## 形式

```txt
dealias <rule>
```

## 例

```sh
dealias none
dealias 3/2
dealias two
```

## 引数

- `<rule>`
  - 型: string
  - 指定可能な値:
    - `none`
    - `off`
    - `three_halves`
    - `3/2`
    - `two`
    - `2`

## デフォルト値

省略時は padding なしです。

```txt
dealias none
```

## 詳細

`dealias` は、非線形項を実空間で評価する際に使用する計算格子の padding ルールを指定します。

このページでいう **active grid** とは、[`grid`](./grid.md) コマンドで指定する基準格子です。
active grid は、時間発展させるスペクトル mode の範囲を決めます。
`dealias` によって増える格子点は、非線形項評価のための padding 領域であり、時間発展させる自由度そのものではありません。

指定値と計算格子サイズの対応は以下の通りです。

| 指定 | 計算格子サイズ |
| --- | --- |
| `none`, `off` | active grid と同じ |
| `three_halves`, `3/2` | active grid の 3/2 倍 |
| `two`, `2` | active grid の 2 倍 |

例えば、

```sh
grid 64 64
dealias 3/2
```

では、active grid は $64 \times 64$、実際のFFT計算格子は $96 \times 96$ になります。

## 制限・注意

- `dealias` は `grid` で指定した active grid から計算格子を決定します。
- `3/2` ルールでは active grid が偶数である必要があります。
- 現在の `grid` コマンド自体も偶数のみを受け付けます。
- `dealias` を有効にすると、snapshot の出力格子も padding 後の計算格子に対応します。
- 移流項などの非線形項を使う計算では、aliasing の影響を抑えるため `3/2` または `two` の利用を検討してください。

## 関連コマンド

- [`grid`](./grid.md)
- [`fix ... nonlinear`](./fix_nonlinear.md)
- [`measure snapshot`](./measure_snapshot.md)
