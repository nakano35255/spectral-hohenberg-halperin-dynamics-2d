# `measure time_series`

## 目的

`time_series` は、時間発展中のスカラー診断量を1つの時系列ファイルへ出力します。
エネルギー、スカラー場フラックス、運動量フラックス、および Fourier 空間の二乗和を、`target` で個別に選択できます。

## 形式

```txt
measure <ID> time_series <on|off> nevery <integer> file <filename> target <target1> [target2 ...]
```

`target` は必ず最後に指定してください。`target` 以降の語はすべて出力対象として解釈されます。

## 例

```sh
measure ts time_series on nevery 100 file output/time_series.dat target E_T E_K E_psi E_C
measure spec time_series on nevery 100 file output/spectral.dat target |rho|^2 |d_rho|^2 |j|^2 |d_j|^2 |psi[0]|^2 |d_psi[0]|^2
measure flx time_series on nevery 100 file output/flux.dat target Jpsi[0]_x Jpsi[0]_y pi_xx pi_xy pi_yy
measure ts time_series off
```

## 引数

| key | 型 | 必須 | 説明 |
| --- | --- | --- | --- |
| `nevery` | integer | yes | 何ステップごとに出力するか |
| `file` | string | yes | 出力ファイル名 |
| `target` | list | yes | 出力する診断量の一覧 |

`nevery` は正の整数である必要があります。
出力先ディレクトリは自動作成されません。

## 出力形式

出力ファイルは空白区切りのテキストです。
1行目は列名を表す header です。

```txt
# step time <target1> <target2> ...
```

各データ行には、ステップ番号、時刻、指定した target の値が `target` で指定した順番のまま出力されます。

## Target

### Energy

| target | 説明 |
| --- | --- |
| `E_T` | 全エネルギー `E_K + E_psi + E_C` |
| `E_K` | 運動エネルギー |
| `E_psi` | スカラー場の自由エネルギー寄与 |
| `E_C` | 圧縮性エネルギー |

`E_psi` は、自由エネルギーモデルから計算されるスカラー場 sector のエネルギーです。
`quiescent` モードでは流体場を時間発展しないため、`E_K` と `E_C` は通常 `0` になります。

### Flux and Stress

| target | 説明 |
| --- | --- |
| `Jpsi[N]_x` | スカラー場成分 `N` の x 方向フラックスの空間平均 |
| `Jpsi[N]_y` | スカラー場成分 `N` の y 方向フラックスの空間平均 |
| `pi_xx` | 運動量フラックス/stress の xx 成分の空間平均 |
| `pi_xy` | 運動量フラックス/stress の xy 成分の空間平均 |
| `pi_yy` | 運動量フラックス/stress の yy 成分の空間平均 |

これらは、時間積分中に評価された flux buffer の Fourier 空間ゼロ mode を、計算格子点数 `Nx * Ny` で割った値です。
`Jpsi[N]_x` と `Jpsi[N]_y` では、`0 <= N < order_parameters` である必要があります。

`fix ... force/sine` や `fix ... force/gradient` によって加えられた外力そのものは、保存型 flux としては記録されません。
そのため、上記の flux/stress target には外力そのものは含まれません。

### Spectral Diagnostics

| target | 説明 |
| --- | --- |
| `|rho|^2` | 密度場の Fourier 空間二乗和 |
| `|d_rho|^2` | 密度場の勾配二乗和 |
| `|j|^2` | 運動量場の Fourier 空間二乗和 |
| `|d_j|^2` | 運動量場の勾配二乗和 |
| `|psi[N]|^2` | スカラー場成分 `N` の Fourier 空間二乗和 |
| `|d_psi[N]|^2` | スカラー場成分 `N` の勾配二乗和 |

spectral diagnostics ではゼロ mode `(kx, ky) = (0, 0)` を除外します。
`d_` が付く target では、各 Fourier mode の二乗和に `k^2` を掛けた量を出力します。

## 制限・注意

- `target` は必ず最後に指定してください。
- `file` には空白を含めないでください。
- 同じ `<ID>` で `measure ... off` を指定すると、その measure は停止します。
- 同じ `<ID>` で別の `measure ... on` を指定すると、既存の measure は新しい設定に置き換えられます。
- `Jpsi[N]` と `psi[N]` 系の target は、`order_parameters 0` では使用できません。
- `quiescent` モードでは運動量場を時間発展しないため、`pi_xx`, `pi_xy`, `pi_yy`, `|j|^2`, `|d_j|^2` は通常の流体場診断としては使用しないでください。

## 関連コマンド

- [`run`](./run.md)
- [`measure snapshot`](./measure_snapshot.md)
- [`measure ave/profile`](./measure_ave_profile.md)
- [`fix ... force/sine`](./fix_force_sine.md)
- [`fix ... force/gradient`](./fix_force_gradient.md)
- [`time_evolution`](./time_evolution.md)
