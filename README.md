# Spectral Hohenberg-Halperin Dynamics 2D

本リポジトリは、2次元周期領域上の Hohenberg-Halperin 型の確率的場のダイナミクスと、等温の揺らぐ Navier-Stokes 型ダイナミクスを、Fourier 擬スペクトル法で時間発展させるための C++/MPI ソルバーです。

## 概要

本ソルバーは、2次元周期境界条件

```math
\Omega = [0, L_x) \times [0, L_y)
```

のもとで、以下の場を時間発展させます。

- 密度場 $\rho(\boldsymbol{x}, t)$
- 運動量密度場 $\boldsymbol{j}(\boldsymbol{x}, t)$ 
- $N$ 成分スカラー場 $\psi_{\alpha}(\boldsymbol{x}, t)$ ( $\alpha = 0,1,\cdots,N-1$ )

具体的には、以下の確率偏微分方程式を数値的に解きます。

```math
\partial_t \rho = -\nabla \cdot \boldsymbol{j}
```
```math
\partial_t \boldsymbol{j} = -\nabla \cdot (\boldsymbol{j} \otimes \boldsymbol{v}) - \nabla p + \eta \nabla^2 \boldsymbol{v} + \zeta \nabla (\nabla \cdot \boldsymbol{v}) + \nabla \cdot \widetilde{\boldsymbol{\Sigma}}
```
```math
\partial_t \psi_{\alpha} = -\nabla \cdot (\psi_{\alpha} \boldsymbol{v}) + M_{\alpha} \nabla^2 \mu_{\alpha} + \nabla \cdot \widetilde{\boldsymbol{F}}_{\alpha}
```

ここで、速度場 $\boldsymbol{v}(\boldsymbol{x},t)$ は

```math
\boldsymbol{v}(\boldsymbol{x},t) = \frac{\boldsymbol{j}(\boldsymbol{x},t)}{\rho(\boldsymbol{x},t)}
```

によって定義されます。

現在の実装では、スカラー場の移動度 $M_{\alpha}$、せん断粘性 $\eta$、体積粘性 $\zeta$ は空間的に一様な定数として扱い、ノイズは加法的な白色ガウスノイズとして導入します。
スカラー場の確率フラックスの統計性は以下で与えられます。

```math
\left\langle \widetilde{F}_{\alpha,a}(\boldsymbol{x},t) \widetilde{F}_{\beta,b}(\boldsymbol{x}',t') \right\rangle
=
2 k_B T M_{\alpha} \delta_{\alpha\beta}\delta_{ab} \delta(\boldsymbol{x}-\boldsymbol{x}') \delta(t-t')
```

確率応力テンソルの統計性は以下で与えられます。

```math
\left\langle \widetilde{\Sigma}_{ab}(\boldsymbol{x},t) \widetilde{\Sigma}_{cd}(\boldsymbol{x}',t') \right\rangle
=
2 k_B T \left[\eta(\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc}) + (\zeta-\eta)\delta_{ab}\delta_{cd} \right]
\delta(\boldsymbol{x}-\boldsymbol{x}') \delta(t-t')
```

ここで、$a,b \in \{x,y\}$ は空間方向を表します。

圧力 $p$ および各成分の化学ポテンシャル $\mu_{\alpha}$ は、ユーザーがコードに独自の定義を容易に追加できる設計となっています。
また、時間積分スキーム、初期条件、および観測量の出力は、入力スクリプトから柔軟に設定・切り替えが可能です。

> [!NOTE]
> **詳細なマニュアル、各コマンドの仕様、チュートリアルについては [Wiki (WIKI.md)](WIKI.md) をご覧ください。**


## 主な機能

- 2次元周期境界条件
- 多成分オーダーパラメータ対応（スカラー場が0成分の場合は揺らぐNavier-Stokes方程式に帰着）
- 圧力 / 化学ポテンシャルのユーザー定義および追加実装が容易なアーキテクチャ
- 圧縮性 / 非圧縮性 / 静止流体モードの切り替え
- heFFTeを用いたMPI並列
- Fourierスペクトル法による空間微分
- 3/2ルール（padding）を用いた擬スペクトル法による非線形項評価
- 時間積分スキーム: Euler法、SRK3法（確率的Runge-Kutta法）
- 保存型加法ノイズ
- 観測量出力: `snapshot`
- 今後の実装予定 (Future works): リスタート機能（状態の保存と読み込み）


## 必要環境

- C++17対応コンパイラ
- MPI
- FFTW3
- heFFTe

必要環境の準備については [dev_env/README.md](dev_env/README.md) も参照してください。
Dockerを用いたローカル開発環境の構築セットも dev_env ディレクトリに同梱されています。


## ビルド

ビルド方法の詳細は [docs/build.md](docs/build.md) に記載しています。

通常はリポジトリ直下で以下を実行します。

```sh
make
```

このコマンドにより、`mpic++` を用いて実行ファイル `src/out.exe` が生成されます。


## 実行

作成した実行ファイルに、入力スクリプトのパスを引数として渡して実行します。

```sh
./src/out.exe path/to/input.script
```

MPIを用いた並列計算を行う場合は、以下のように実行します（例: 2プロセス）。

```sh
mpirun -np 2 ./src/out.exe path/to/input.script
```

## ドキュメント

- [WIKI.md](WIKI.md)
- [方程式と符号規約](docs/equations.md)
- [ビルド方法](docs/build.md)
- [入力スクリプト仕様](docs/input-script/index.md)
- [新しい機能の追加手順](docs/developer/user-add.md)


## ライセンス

このプロジェクトは MIT License のもとで公開されています。詳細は [LICENSE](LICENSE) を参照してください。
