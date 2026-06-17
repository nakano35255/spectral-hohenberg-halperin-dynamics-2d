# `measure budget/spectrum`

## 目的

`budget/spectrum` は、指定したオーダーパラメータ成分 $\psi_a$ の Fourier mode ごとに、決定論的な時間発展項が $\lvert \psi_a(\boldsymbol{k}) \rvert^2$ の変化へどのように寄与するかをスペクトルとして出力します。

オーダーパラメータ方程式の右辺を

```math
\partial_t \psi_a(\boldsymbol{k})
= T_a(\boldsymbol{k}) + D_a(\boldsymbol{k}) + P_a(\boldsymbol{k})
```

のように分け、それぞれに対して

```math
B_X(\boldsymbol{k})
= \frac{w(\boldsymbol{k})}{N_g^2}
  \operatorname{Re}
  \left[
    \psi_a(\boldsymbol{k})^\ast X_a(\boldsymbol{k})
  \right]
```

を時間平均して出力します。
ここで $X$ は `transfer`, `dissipation`, `production` のいずれか、$N_g = N_x N_y$ は計算格子点数、$w(\boldsymbol{k})$ は R2C 表現の共役 mode を含めるための重みです。
ゼロ mode は出力しません。

各項の意味は以下です。

- `transfer`
  - 速度場によるオーダーパラメータ移流項です。
  - `fix ... order_parameter nonlinear on` が有効なとき、
    $T_a(\boldsymbol{k}) = -i\boldsymbol{k}\cdot\mathcal{F}\{\psi_a\boldsymbol{v}\}$ を評価します。
- `dissipation`
  - chemical potential による散逸項です。
  - $D_a(\boldsymbol{k}) = -M_a k^2 \mu_a(\boldsymbol{k})$ を評価します。
  - linear chemical potential と、自由エネルギーモデルが持つ実空間 chemical potential の両方に対応します。
- `production`
  - 外力による生成項です。
  - 現在は、対象成分に作用する `fix ... order_parameter force/sine` と `fix ... order_parameter force/gradient` を含めます。

`mode 2d` では非等方性を保ったまま $\boldsymbol{k}$ ごとの値を出力します。
`mode shell` では同じ $|\boldsymbol{k}|$ shell に入る mode を足し合わせ、角度積分した spectrum として出力します。

> [!NOTE]
> この measure は `PASSIVE_SCALAR` package に含まれます。
> 使用するには、ビルド前に package を有効化してください。
>
> ```sh
> make yes-PASSIVE-SCALAR
> make clean
> make
> ```
>
> 一時的に有効化してビルドする場合は、次のように指定することもできます。
>
> ```sh
> make PKG_PASSIVE_SCALAR=1
> ```

実行例として、[`examples/03_2`](../../examples/03_2/README.md) に小さい steady-state check 用の入力スクリプトと解析スクリプトがあります。

## 形式

```txt
measure <ID> budget/spectrum <on|off> component <integer> nevery <integer> nblock <integer> file <filename> mode <2d|shell> average <block|running>
```

## 例

```sh
measure bs_shell budget/spectrum on component 0 nevery 20 nblock 200 file output/budget_shell.dat mode shell average block
measure bs_2d budget/spectrum on component 0 nevery 50 nblock 250 file output/budget_2d.dat mode 2d average running

measure bs_shell budget/spectrum off
```

定常状態に入ってからだけ測定したい場合は、先に測定なしで `run` し、その後に `measure budget/spectrum` を有効化します。

```sh
run                 1000

measure             bs_shell budget/spectrum on component 0 nevery 20 nblock 200 file output/budget_shell.dat mode shell average block

run                 1000
```

## 引数

- `<ID>`
  - 型: string
  - measure の識別子です。
  - 同じ `<ID>` を再指定すると、既存の measure は置き換えられます。
- `<on|off>`
  - 型: string
  - `on` で有効化、`off` で無効化します。
- `component <integer>`
  - 対象にするオーダーパラメータ成分を指定します。
  - `0 <= component < order_parameters` である必要があります。
- `nevery <integer>`
  - 何ステップごとに観測するかを指定します。
- `nblock <integer>`
  - 何ステップごとに平均スペクトルを出力するかを指定します。
  - `nblock` は `nevery` で割り切れる必要があります。
- `file <filename>`
  - 出力ファイル名を指定します。
- `mode <2d|shell>`
  - 出力するスペクトルの形式を指定します。
  - `2d` は active spectral mask 内の各 Fourier mode を出力します。
  - `shell` は同じ shell に入る Fourier mode の寄与を足し合わせて出力します。
- `average <block|running>`
  - 時間平均の出力方法を指定します。
  - `block` は、各 block の平均スペクトルをファイル末尾へ追記します。
  - `running` は、これまでに完了した block すべての running average をファイルへ書き直します。

## デフォルト値

`measure budget/spectrum` を指定しない場合、budget spectrum は出力されません。

`on` の場合、`component`, `nevery`, `nblock`, `file`, `mode`, `average` はすべて明示的に指定する必要があります。

## 出力タイミング

`measure budget/spectrum` は、スクリプト内で有効化された後の `run` に対して作用します。
観測は各ステップの時間発展後に行われます。
`nevery` と `nblock` は、通算ステップ番号ではなく、その measure が有効化されてからの内部カウンタで判定されます。
同じ measure が有効であり続ける間は、複数の `run` にまたがって同じ block が継続します。
同じ `<ID>` の measure を `off` または再度 `on` にすると、その measure の block は終了または作り直されます。

`nblock` ステップ分の観測が終わると、平均された budget spectrum を出力します。
例えば、

```txt
nevery 20 nblock 200
```

では、20 step ごとに観測し、10個の観測値を平均して1つの block spectrum を作ります。

## 出力内容

出力ファイルは rank 0 が書き込みます。
最初に共通 header が出力されます。

```txt
# measure budget/spectrum
# component <component> nevery <integer> nblock <integer> mode <2d|shell> average <block|running> normalization per_volume
```

各 block は、次のコメント行から始まります。

```txt
# block <block> step <step> time <time> samples <samples>
```

`samples` は、その出力に含まれる時間サンプル数です。
`average block` では、その block に含まれるサンプル数です。
`average running` では、これまでに完了した block すべてに含まれる累積サンプル数です。

### `mode 2d`

```txt
# kx ky transfer dissipation production total
```

各データ行は、active spectral mask 内の非ゼロ Fourier mode に対応します。

| column | 意味 |
| --- | --- |
| `kx` | x 方向波数 |
| `ky` | y 方向波数 |
| `transfer` | 移流項の寄与 |
| `dissipation` | chemical potential 散逸項の寄与 |
| `production` | 外力生成項の寄与 |
| `total` | `transfer + dissipation + production` |

R2C 表現のため、出力される mode は $k_x \ge 0$ 側の表現です。
`kx > 0` の行には、共役 mode の寄与を含める重みがすでに掛かっています。

### `mode shell`

```txt
# k count transfer dissipation production total
```

各データ行は、1つの shell に対応します。

| column | 意味 |
| --- | --- |
| `k` | shell の代表波数 |
| `count` | R2C 重みを含めた shell 内 mode 数 |
| `transfer` | shell 内で足し合わせた移流項の寄与 |
| `dissipation` | shell 内で足し合わせた chemical potential 散逸項の寄与 |
| `production` | shell 内で足し合わせた外力生成項の寄与 |
| `total` | `transfer + dissipation + production` |

shell index は

```math
n = \left\lfloor \frac{|\boldsymbol{k}|}{\Delta k} + \frac{1}{2} \right\rfloor,
\qquad
\Delta k = \min\left(\frac{2\pi}{L_x}, \frac{2\pi}{L_y}\right)
```

で決まり、出力される `k` は $n\Delta k$ です。

`mode shell` の `transfer`, `dissipation`, `production`, `total` は shell average ではなく shell sum です。
shell 内 mode あたりの平均が必要な場合は、後処理で `count` で割ってください。

## 計算される項

### Transfer

`transfer` は、`fix ... order_parameter nonlinear on` が指定されている場合に計算されます。
実装では保存形の移流項

```math
T_a(\boldsymbol{k})
= -i\boldsymbol{k}\cdot\mathcal{F}
  \left[
    \psi_a(\boldsymbol{x})\boldsymbol{v}(\boldsymbol{x})
  \right]
```

を用います。
`fix ... order_parameter nonlinear off` の場合、`transfer` は 0 です。

速度場は時間発展 mode に応じて以下のように評価されます。

- incompressible mode:

```math
\boldsymbol{v}(\boldsymbol{x})
= \boldsymbol{j}(\boldsymbol{x}) / \rho_0
```

- compressible mode:

```math
\boldsymbol{v}(\boldsymbol{x})
= \boldsymbol{j}(\boldsymbol{x}) / \rho(\boldsymbol{x})
```

### Dissipation

`dissipation` は、オーダーパラメータ mobility と chemical potential から計算されます。

```math
D_a(\boldsymbol{k})
= -M_a k^2 \mu_a(\boldsymbol{k})
```

linear chemical potential については、

```math
\mu_a(\boldsymbol{k})
=
\left(
  k_0 + k_2 k^2 + k_4 k^4
\right)
\psi_a(\boldsymbol{k})
```

を用います。
自由エネルギーモデルが実空間 chemical potential を持つ場合は、その Fourier transform も加えます。

### Production

`production` は、対象成分に作用する外力項から計算されます。

- `fix ... order_parameter force/sine`
  - sine force の Fourier 成分を直接加えます。
- `fix ... order_parameter force/gradient`
  - 指定された方向の速度 Fourier 成分に比例する項を加えます。

対象成分に外力が指定されていない場合、`production` は 0 です。

## 制限・注意

- `budget/spectrum` を使うには、`PASSIVE_SCALAR` package を有効にしてビルドする必要があります。
- package が無効な実行ファイルでは、`budget/spectrum` は未登録の measure type になります。
- `nevery` は正の整数である必要があります。
- `nblock` は正の整数であり、`nevery` で割り切れる必要があります。
- `file` には空白を含めないでください。
- 出力先ディレクトリは自動作成されません。
- `file` は snapshot のような prefix ではなく、出力ファイル名そのものです。
- `average` は `block` または `running` のいずれかです。
- `mode` は `2d` または `shell` のいずれかです。
- `order_parameters 0` の場合、この measure は使用できません。
- `component` は存在するオーダーパラメータ成分を指定してください。
- 現在の実装では、オーダーパラメータ noise の寄与は budget に含めていません。そのため、`fix ... order_parameter noise on` と併用するとエラーになります。
- `transfer` または `force/gradient` による `production` は速度場を必要とするため、quiescent mode では使用できません。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。

## 関連コマンド

- [`run`](./run.md)
- [`measure time_series`](./measure_time_series.md)
- [`measure ave/profile`](./measure_ave_profile.md)
- [`fix order_parameter nonlinear`](./fix_nonlinear.md)
- [`fix ... force/sine`](./fix_force_sine.md)
- [`fix ... force/gradient`](./fix_force_gradient.md)
- [`time_evolution`](./time_evolution.md)
