# `run`

## 目的

指定したステップ数だけ時間発展を実行します。

## 形式

```txt
run <steps>
```

## 例

```sh
run 10000
run 200000
```

## 引数

- `<steps>`
  - 型: integer
  - 実行する時間ステップ数
  - 正の整数である必要があります。

## 詳細

`run` は、読み込まれたシステム設定、物理モデル、初期条件に基づいて、時間発展を実行します。
`measure` コマンドについては、スクリプト内でその `run` より前に有効化されているものが観測されます。

複数回 `run` を書いた場合、ステップ番号と物理時刻はリセットされず、通算で進みます。

```sh
run 10000
run 10000
```

この場合、合計で20,000ステップが実行されます。

`measure` コマンドと組み合わせることで、途中から観測量を有効化または停止できます。

```sh
run 10000

measure phys snapshot on nevery 1000 file output/physical space physical
run 10000

measure phys snapshot off
run 10000
```

この例では、最初の10,000ステップでは観測を行わず、次の10,000ステップで `phys` を出力し、最後の10,000ステップでは `phys` を停止します。

## 制限・注意

- `<steps>` は正の整数である必要があります。
- `0` 以下の値を指定するとエラーになります。
- `run` より前に、少なくとも `model transport` を指定する必要があります。
- `fix` コマンドは、現在の実装では `run` または `measure` より前に指定する必要があります。
- 初期条件は、最初の `run` が始まる前に一度だけ状態へ適用されます。
- 複数の `run` を連続して実行しても、状態は途中で再初期化されません。

## 関連コマンド

- [`timestep`](./timestep.md)
- [`time_evolution`](./time_evolution.md)
- [`measure snapshot`](./measure_snapshot.md)
- [`thermo`](./output.md#thermo)
