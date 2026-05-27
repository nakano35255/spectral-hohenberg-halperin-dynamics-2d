# `model ... coeff`

## 目的

すでに指定済みの `model` に対して、係数だけを追加・更新します。

## 形式

```txt
model <category> coeff [key value ...]
```

## 例

```sh
model thermo linear_eos
model thermo coeff cT 10.0

model free_energy phi4
model free_energy coeff k0[0] -3.0 k2[0] 1.0 phi4[0] 1.0

model transport constant eta 1.0
model transport coeff zeta 1.0 M[0,0] 0.5
```

## 引数

- `<category>`
  - 型: string
  - 指定可能な値:
    - `thermo`
    - `free_energy`
    - `transport`
- `[key value ...]`
  - 対象カテゴリで使える係数を key-value 形式で指定します。

## 詳細

`model ... coeff` は、モデルの種類を変更せず、係数だけを更新します。
指定した key の値だけが上書きされ、指定しなかった係数は現在の値を保ちます。

一方で、

```sh
model transport constant eta 1.0
model transport constant zeta 1.0
```

のように `model <category> <type>` を再度指定すると、そのカテゴリの設定はデフォルト値から作り直されます。
この例では、2行目で `eta` はデフォルト値 `0.0` に戻ります。

既存のモデルに係数を足していきたい場合は、次のように書きます。

```sh
model transport constant eta 1.0
model transport coeff zeta 1.0
```

## 制限・注意

- `model <category> coeff ...` を使う前に、対応する `model <category> <type>` を指定しておく必要があります。
- 指定できる key は、選択済みの model type に依存します。
- key-value の組が崩れている場合はエラーになります。

## 関連コマンド

- [`model thermo`](./model_thermo.md)
- [`model free_energy`](./model_free_energy.md)
- [`model transport`](./model_transport.md)
