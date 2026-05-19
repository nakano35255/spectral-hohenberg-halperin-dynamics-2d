# 初期条件を追加する手順

このページは、新しい初期条件を追加するときに「どのファイルを作り、どこを書き換えるか」を示すテンプレートです。

例として、density 初期条件 `hoge` を追加します。

input script では次のように使う想定です。

```text
set density 0 hoge value 1.0
set density all hoge value 1.0
```

## 作るファイル

density 初期条件なら、この 3 ファイルを作ります。

```text
src/initial_condition_hoge_density_style.h
src/initial_condition_hoge_density.h
src/initial_condition_hoge_density.cc
```

momentum 初期条件なら、`density` を `momentum` に置き換えます。

```text
src/initial_condition_hoge_momentum_style.h
src/initial_condition_hoge_momentum.h
src/initial_condition_hoge_momentum.cc
```

## 1. style.h を作る

`src/initial_condition_hoge_density_style.h`

このファイルには、parser 周りを書きます。

主に編集する場所は次です。

- `TYPE_NAME`: input script で使う名前
- `HogeDensityInitialConditionCommand`: parse 後に `Params` へ保存する値
- `parse_command()`: `value` などの引数を読む
- `create_initial_condition()`: 実装クラスを作る

```cpp
#ifndef SFI_INITIAL_CONDITION_HOGE_DENSITY_STYLE_H
#define SFI_INITIAL_CONDITION_HOGE_DENSITY_STYLE_H

#include "initial_condition_registry.h"
#include "initial_condition_hoge_density.h"

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>

// EDIT: input script で使う type 名。
// CommandBase::type と Style::type_name() は必ず同じ文字列にする。
namespace hoge_density_initial_condition {
    inline const std::string TYPE_NAME = "hoge";
}

// EDIT: parser が読んだ情報を保存する command。
// Project 側では、この command が Params.initial.density_commands に残る。
struct HogeDensityInitialConditionCommand : DensityInitialConditionCommandBase {
    double value = 0.0;

    HogeDensityInitialConditionCommand(int component_in, double value_in)
        : value(value_in) {
        type = hoge_density_initial_condition::TYPE_NAME;
        component = component_in;
    }

    void print(std::ostream& os) const override {
        os << "  Density: Component " << component
           << "              : " << type
           << " value " << value << '\n';
    }
};

// EDIT: parser と registry から呼ばれる style。
class HogeDensityInitialConditionStyle : public DensityInitialConditionStyle {
public:
    const std::string& type_name() const override {
        return hoge_density_initial_condition::TYPE_NAME;
    }

    std::shared_ptr<DensityInitialConditionCommandBase> parse_command(
        int component,
        const InitialConditionArgs& args,
        const Params& params
    ) const override {
        (void)params;

        // EDIT: input script の key を読む。
        // 例: set density 0 hoge value 1.0
        const double value = std::stod(args.get_required("value"));

        // EDIT: 必要なら入力値をチェックする。
        if (value < 0.0) {
            throw std::runtime_error("set density hoge requires nonnegative value");
        }

        return std::make_shared<HogeDensityInitialConditionCommand>(
            component,
            value
        );
    }

    std::unique_ptr<DensityInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    ) const override {
        return std::make_unique<HogeDensityInitialCondition>(params, command);
    }
};

#endif
```

`TYPE_NAME` は 2 か所で使います。

- `HogeDensityInitialConditionCommand::type`: parse 後に保存された command が、自分の type を覚えるため
- `HogeDensityInitialConditionStyle::type_name()`: registry に登録するときの名前

この 2 つがずれると、parse 後に command から style を取り直すときに壊れます。

## 2. initial_condition.h を作る

`src/initial_condition_hoge_density.h`

このファイルには、実際に `State` へ書き込むクラスの宣言を書きます。

主に編集する場所は次です。

- member 変数: command から受け取った値を保持する
- `apply()`: spectral `State` に初期条件を書き込む

```cpp
#ifndef SFI_INITIAL_CONDITION_HOGE_DENSITY_H
#define SFI_INITIAL_CONDITION_HOGE_DENSITY_H

#include "initial_condition.h"

#include <memory>

class HogeDensityInitialCondition : public DensityInitialCondition {
private:
    // EDIT: 初期条件に必要な値を持つ。
    double value_ = 0.0;

public:
    HogeDensityInitialCondition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    );

    void apply(State& state, const Domain2D& domain) const override;
};

#endif
```

## 3. initial_condition.cc を作る

`src/initial_condition_hoge_density.cc`

このファイルには、command から値を取り出す処理と、`State` への書き込みを書きます。

主に編集する場所は次です。

- constructor: command を具体型へ cast して値を読む
- `apply()`: spectral coefficient を書く
- MPI 分割: その rank が対象 mode を持っているか確認する

```cpp
#include "initial_condition_hoge_density.h"
#include "initial_condition_hoge_density_style.h"

#include <stdexcept>

HogeDensityInitialCondition::HogeDensityInitialCondition(
    const Params& params,
    std::shared_ptr<const DensityInitialConditionCommandBase> command
) : DensityInitialCondition(params, command) {
    (void)params;

    const auto cfg =
        std::dynamic_pointer_cast<const HogeDensityInitialConditionCommand>(command);
    if (!cfg) {
        throw std::runtime_error("HogeDensityInitialCondition: invalid command type.");
    }

    value_ = cfg->value;
}

void HogeDensityInitialCondition::apply(State& state, const Domain2D& domain) const {
    const Box2D& box = domain.spectral_box();

    // EDIT: 書き込みたい spectral mode を決める。
    const int kx = 0;
    const int ky = 0;

    // EDIT: MPI 分割されているので、この rank が mode を持つときだけ書く。
    if (box.low[0] <= kx && kx <= box.high[0] &&
        box.low[1] <= ky && ky <= box.high[1]) {
        // EDIT: density なら rho_hat に書く。
        state.rho_hat(component_, kx, ky) = Complex(value_, 0.0);
    }
}
```

空間一様な density の場合は、zero mode に `value * Nx * Ny` を書きます。

```cpp
const double amplitude =
    value_
    * static_cast<double>(domain.nx_global())
    * static_cast<double>(domain.ny_global());

state.rho_hat(component_, 0, 0) = Complex(amplitude, 0.0);
```

`State::clear()` を呼んでから初期条件を適用する前提なら、書かない mode は 0 のままです。

## 4. registry に登録する

`src/initial_condition_registry_builtin.cc`

作った style を include して、registry に登録します。

```cpp
#include "initial_condition_registry_builtin.h"
#include "initial_condition_uniform_density_style.h"
#include "initial_condition_hoge_density_style.h"

InitialConditionRegistry build_initial_condition_registry() {
    InitialConditionRegistry registry;

    registry.register_density_style(std::make_unique<UniformDensityInitialConditionStyle>());
    registry.register_density_style(std::make_unique<HogeDensityInitialConditionStyle>());

    return registry;
}
```

momentum 初期条件なら、こちらを使います。

```cpp
registry.register_momentum_style(std::make_unique<HogeMomentumInitialConditionStyle>());
```

## 5. Makefile に .cc を追加する

`Makefile`

`SRCS` に実装ファイルを追加します。

```make
SRCS := \
    src/main.cc \
    src/initial_condition_registry_builtin.cc \
    src/initial_condition_uniform_density.cc \
    src/initial_condition_hoge_density.cc \
    ...
```

`style.h` は header-only なので、通常は `style.cc` を作りません。

## Momentum の場合

density から momentum に変えるときは、次を置き換えます。

```text
DensityInitialConditionCommandBase -> MomentumInitialConditionCommandBase
DensityInitialConditionStyle       -> MomentumInitialConditionStyle
DensityInitialCondition            -> MomentumInitialCondition
register_density_style             -> register_momentum_style
params.initial.density_commands    -> params.initial.momentum_commands
```

`apply()` では density なら `rho_hat`、momentum なら momentum 用の spectral component に書きます。

## 最低限の確認リスト

- `TYPE_NAME` が input script の type 名と同じ
- command constructor で `type` と `component` をセットしている
- `style.type_name()` が同じ `TYPE_NAME` を返している
- `parse_command()` で必要な key を読んでいる
- constructor で command を正しい具体型へ cast している
- `apply()` が physical ではなく spectral `State` に直接書いている
- 対象 mode が local `spectral_box()` に含まれるときだけ書いている
- registry に style を登録している
- `Makefile` に `.cc` を追加している
