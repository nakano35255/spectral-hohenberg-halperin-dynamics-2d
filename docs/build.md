# ビルド方法

このページでは、Spectral Hohenberg-Halperin Dynamics 2D のビルド方法をまとめます。
通常の開発では Docker / Dev Containers を使う方法を推奨します。
ローカル環境で直接ビルドする手順は現在工事中です。


## 必要なライブラリ

本ソルバーのビルドには以下が必要です。

- C++17 対応コンパイラ
- MPI 実装
- FFTW3
- FFTW3 MPI interface
- heFFTe

リポジトリ直下の `Makefile` は、以下の配置を前提にしています。

```text
/usr/local/include/heffte.h
/usr/local/lib/libheffte.so
```

FFTW3 はコンパイラの標準探索パス、または linker が見つけられる場所にあることを想定しています。


> [!NOTE]
> 環境確認用の詳細なコマンドは、次のドキュメントにもあります。
> - [dev_env/README.md](../dev_env/README.md)
> - [dev_env/docker/README.md](../dev_env/docker/README.md)
> - [dev_env/heffte/README.md](../dev_env/heffte/README.md)


## 推奨: Docker 環境でビルドする

Dockerfile には、OpenMPI、FFTW3、FFTW3 MPI interface、FFTW backend 付き heFFTe のインストール手順が入っています。

リポジトリ直下で Docker image を作成します。

```sh
docker build -t spectral-hohenberg-halperin-dynamics-2d .
```

作成した image で shell を起動します。

```sh
docker run --rm -it \
  -v "$PWD":/work \
  -w /work \
  spectral-hohenberg-halperin-dynamics-2d \
  bash
```

container の中で、まず開発環境の簡単な確認を行います。

```sh
make -C dev_env/docker clean
make -C dev_env/docker
mpirun --allow-run-as-root -np 4 dev_env/docker/out.exe

make -C dev_env/heffte clean
make -C dev_env/heffte
mpirun --allow-run-as-root -np 2 dev_env/heffte/out.exe
```

期待される出力には、次のような行が含まれます。

```text
FFTW MPI check OK with 4 ranks
heFFTe FFTW check OK with 2 ranks
```

確認が通ったら、ソルバー本体をビルドします。

```sh
make
```

成功すると、`src/out.exe` が生成されます。


## ローカル環境でビルドする

工事中


## 実行確認

example を実行する場合は、実行ファイルに入力スクリプトを明示的に渡します。

```sh
./src/out.exe examples/01_phi4_phase_separation/input.script
```

MPI 並列で実行する場合も同様です。

```sh
mpirun -np 2 ./src/out.exe examples/02_compressible_navier_stokes/input.script
```

Docker container 内で root user として `mpirun` する場合は、OpenMPI の設定により次のように実行します。

```sh
mpirun --allow-run-as-root -np 2 ./src/out.exe examples/02_compressible_navier_stokes/input.script
```


## クリーンビルド

生成された実行ファイルを削除するには、以下を実行します。

```sh
make clean
```

その後、再度ビルドします。

```sh
make
```
