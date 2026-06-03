# 入力スクリプト仕様

本ページでは、Spectral Hohenberg-Halperin Dynamics 2D における入力スクリプトの仕様と記述方法について解説します。
本ソルバーは、入力スクリプトで指定された物理モデル、係数、初期条件、および観測設定に従って、2次元周期領域上の密度場・運動量密度場・スカラー場（オーダーパラメータ）を時間発展させます。

> [!NOTE]
> 初めて入力スクリプトを作成する場合は、`examples` ディレクトリに同梱されている `input.script` をコピーし、必要な箇所を書き換えて使用することをおすすめします。


## 目次

- [入力スクリプトとは](#入力スクリプトとは)
- [基本構文](#基本構文)
- [推奨される記述順序](#推奨される記述順序)
- [コマンド一覧](#コマンド一覧)
  - [System Setup](#system-setup)
  - [Runtime](#runtime)
  - [Model](#model)
  - [Fix](#fix)
  - [Initial Condition](#initial-condition)
  - [Output](#output)
  - [Restart](#restart)
  - [Measure and Run](#measure-and-run)
- [代表例](#代表例)
  - [phi4 phase separation](#phi4-phase-separation)
  - [compressible Navier-Stokes](#compressible-navier-stokes)


## 入力スクリプトとは

入力スクリプトは、シミュレーション対象、数値積分方法、物理モデル、初期条件、観測量、および実行ステップ数を指定するテキストファイルです。

典型的には、以下の内容を順に記述します。

1. 格子数、系の長さ、境界条件の定義
2. スカラー場（オーダーパラメータ）の成分数の指定
3. 時間刻み幅と時間発展モードの選択
4. 圧力、自由エネルギー、輸送係数モデルの設定
5. 非線形移流項とノイズの有効化
6. 各場（密度、運動量、スカラー場）の初期条件の設定
7. ログ出力、出力ファイルの設定
8. `measure` コマンドによる観測量の指定と `run` コマンドによる時間発展の実行

システムやモデルの基本設定を行うコマンド群は、最初の `run` コマンドより前に一度だけ記述します。一方、`measure` コマンドと `run` コマンドは、スクリプトの上から下へ順番に逐次処理されます。


## 基本構文

- **1行1コマンド:** 1行につき1つのコマンドを記述します。
- **区切り文字:** コマンドと引数は空白（スペース）で区切ります。
- **コメント:** 行頭または行中の # 以降の文字列はコメントとして無視されます。
- **空行:** 空行は無視されます。
- **大文字・小文字の区別:** コマンド名、引数名、および on / off などの指定は、大文字と小文字が厳密に区別されます。
- **ファイルパス:** 空白は区切り文字として認識されるため、出力先などのファイルパスには空白を含めないでください。
- **ステップと時刻の通算:** 複数回 run を実行した場合、ステップ番号と物理時刻はリセットされず通算されます。

**記述例:**

```sh
# comment
grid      64 64
length    64.0 64.0
timestep  0.005
```

**measure と run の順序依存性:**
観測設定は、スクリプト内の記述順序に従って有効化・無効化されます。

```sh
measure phys snapshot on nevery 1000 file output/physical space physical
run 10000

measure spec snapshot on nevery 1000 file output/spectral space spectral
run 10000

measure phys snapshot off
run 10000
```

上記の場合、最初の10,000ステップ（1つ目の `run`）では `phys` のみが観測されます。次の10,000ステップでは `phys` と `spec` の両方が観測され、最後の10,000ステップでは `phys` の観測が停止し `spec` のみが観測されます。


## 推奨される記述順序

入力スクリプトは、論理的な依存関係を満たすため、以下の順序で記述することを推奨します。

```sh
# 1. system setup
dimension           2
boundary            p p
grid                64 64
length              64.0 64.0
dealias             3/2
order_parameters    1

# 2. runtime
timestep            0.005
time_evolution      srk3/compressible

# 3. models
model thermo        linear_eos cT 10.0
model free_energy   phi4 a[0] -3.0 b[0] 1.0 u[0] 1.0
model transport     constant eta 1.0 zeta 1.0 M[0,0] 1.0

# 4. optional terms
fix                 1 all nonlinear on
fix                 2 all noise on seed 12345 kBT 1.0

# 5. initial condition
set                 density uniform value 1.0
set                 momentum all uniform value 0.0
set                 order_parameter 0 uniform value 0.0

# 6. outputs and measurements
thermo              observe on progress on nevery 1000
measure             phys snapshot on nevery 1000 file output/physical space physical

# 7. time evolution
run                 100000
```

**重要な依存関係のルール:**

- **`order_parameters` の優先指定:** `model free_energy` や `model transport` は、指定時点のスカラー場成分数に基づいて内部の係数配列を初期化します。必ずこれらのモデル設定より前に `order_parameters` を指定してください。
- **`fix` の指定位置:** `fix` コマンドは、必ず `run` や `measure` より前に指定してください。
- **途中からの観測:** 特定のステップから観測を開始したい場合は、まず観測を無効にした状態で `run`を実行し、その後に `measure ... on` を記述して再度 `run` を実行します。


## コマンド一覧

各コマンドの詳細は、表の「詳細」列のリンク先を参照してください。
(※リンク先の個別ページは順次公開予定です)

### System Setup

系の基本設定（グリッド数や物理サイズなど）を行うコマンド群です。スクリプトの先頭ブロックで指定することが推奨されます。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `dimension` | `dimension 2` | [詳細](./dimension.md) |
| `boundary` | `boundary <x> <y>` | [詳細](./boundary.md) |
| `grid` | `grid <Nx> <Ny>` | [詳細](./grid.md) |
| `length` | `length <Lx> <Ly>` | [詳細](./length.md) |
| `dealias` | `dealias <none\|3/2\|two>` | [詳細](./dealias.md) |
| `order_parameters` | `order_parameters <num>` | [詳細](./order_parameters.md) |


### Runtime

時間刻み幅や積分スキームなど、時間発展に関する基本設定です。System Setup の直後に指定することが推奨されます。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `timestep` | `timestep <dt>` | [詳細](./timestep.md) |
| `time_evolution` | `time_evolution <scheme>/<mode>` | [詳細](./time_evolution.md) |

### Model

数値シミュレーションに用いる物理モデル（圧力、化学ポテンシャル、輸送係数など）を定義します。Runtime の直後に指定することが推奨されます。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `model thermo` | `model thermo linear_eos cT <value>` | [詳細](./model_thermo.md) |
| `model free_energy` | `model free_energy <quadratic\|phi4> [key value ...]` | [詳細](./model_free_energy.md) |
| `model transport` | `model transport constant [key value ...]` | [詳細](./model_transport.md) |
| `model ... coeff` | `model <category> coeff [key value ...]` | [詳細](./model_coeff.md) |

### Fix

時間発展方程式に対する非線形項、確率的ノイズ、外力の有効化・無効化を制御します。Model の直後に指定することが推奨されます。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `fix ... nonlinear` | `fix <ID> <target> nonlinear <on\|off>` | [詳細](./fix_nonlinear.md) |
| `fix ... noise` | `fix <ID> <target> noise <on\|off> [seed <integer>] [kBT <value>] [chi <value>]` | [詳細](./fix_noise.md) |
| `fix ... force/sine` | `fix <ID> <target> force/sine <on\|off> component <component> axis <axis> nk <integer> amplitude <value>` | [詳細](./fix_force_sine.md) |
| `fix ... force/gradient` | `fix <ID> order_parameter force/gradient <on\|off> component <component> direction <direction> amplitude <value>` | [詳細](./fix_force_gradient.md) |

### Initial Condition

各場の初期状態を定義します。Fix の直後に指定することが推奨されます。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `set density` | `set density <uniform\|sine\|gaussian> [key value ...]` | [詳細](./set_density.md) |
| `set momentum` | `set momentum <x\|y\|0\|1\|all> <uniform\|sine> [key value ...]` | [詳細](./set_momentum.md) |
| `set order_parameter` | `set order_parameter <component\|all> <uniform\|sine\|gaussian> [key value ...]` | [詳細](./set_order_parameter.md) |

### Output

シミュレーションの進捗表示（標準出力）の頻度などを設定します。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `thermo` | `thermo observe <on\|off> progress <on\|off> nevery <integer>` | [詳細](./output.md#thermo) |

### Restart

リスタートファイルの設定を行います。（※実際のファイル保存・読み込み処理は今後のアップデートで実装予定です）

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `restart` | `restart off` / `restart on file <filename>` | [詳細](./restart.md#restart) |

### Measure and Run

観測量のファイル出力設定と、シミュレーションの実行（時間発展）を指示します。

| コマンド | 形式 | 詳細 |
| --- | --- | --- |
| `measure snapshot` | `measure <ID> snapshot <on\|off> [key value ...]` | [詳細](./measure_snapshot.md) |
| `measure energetics` | `measure <ID> energetics <on\|off> [key value ...]` | [詳細](./measure_energetics.md) |
| `measure flux` | `measure <ID> flux <on\|off> [key value ...]` | [詳細](./measure_flux.md) |
| `measure ave/profile` | `measure <ID> ave/profile <on\|off> axis <x\|y> nevery <integer> nblock <integer> file <filename> average <block\|running> target <target1> [target2 ...]` | [詳細](./measure_ave_profile.md) |
| `run` | `run <steps>` | [詳細](./run.md) |


## 代表例

### phi4 phase separation

静止流体モードを用いて、1成分の保存型Cahn-Hilliard方程式（相分離ダイナミクス）を計算するスクリプト例です。

```sh
dimension           2
boundary            p p

grid                64 64
length              64.0 64.0
dealias             two

order_parameters    1

timestep            0.005
time_evolution      srk3/quiescent

model free_energy   phi4 a[0] -3.0 b[0] 1.0 u[0] 1.0
model transport     constant M[0,0] 1.0

fix                 1 all noise on seed 12345 kBT 1.0

set                 order_parameter 0 uniform value 0.0

thermo              observe on progress on nevery 1000
measure             phys snapshot on nevery 1000 file examples/02_phi4_phase_separation/phase_separation_physical space physical
measure             spec snapshot on nevery 1000 file examples/02_phi4_phase_separation/phase_separation_spectral space spectral
restart             off

run                 200000
```

### compressible Navier-Stokes

スカラー場を除外し、等温の揺らぐ圧縮性Navier-Stokes方程式を計算するスクリプト例です。

```sh
dimension           2
boundary            p p

grid                64 64
length              64.0 64.0
dealias             3/2

order_parameters    0

timestep            0.005
time_evolution      srk3/compressible

model thermo        linear_eos cT 10.0
model transport     constant eta 1.0 zeta 1.0

fix                 1 momentum nonlinear on
fix                 2 all noise on kBT 1.0 seed 12345

set                 density uniform value 1.0
set                 momentum x uniform value 0.0
set                 momentum y uniform value 0.0

thermo              observe on progress on nevery 1000
measure             phys snapshot on nevery 1000 file examples/09_compressible_navier_stokes/compressible_ns_physical space physical
measure             spec snapshot on nevery 1000 file examples/09_compressible_navier_stokes/compressible_ns_spectral space spectral
restart             off

run                 50000
```
