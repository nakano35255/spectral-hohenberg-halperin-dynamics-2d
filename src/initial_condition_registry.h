#ifndef SFI_INITIAL_CONDITION_REGISTRY_H
#define SFI_INITIAL_CONDITION_REGISTRY_H

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Params;
class DensityInitialCondition;
class MomentumInitialCondition;

struct InitialConditionCommandBase {
    std::string type;
    virtual ~InitialConditionCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;
};

struct DensityInitialConditionCommandBase : InitialConditionCommandBase {
    int component = -1;
};

struct MomentumInitialConditionCommandBase : InitialConditionCommandBase {
    int component = -1;
};


struct InitialConditionArgs {
    std::vector<std::pair<std::string, std::string>> entries;

    std::string get_required(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        throw std::runtime_error("Missing required initial condition argument: " + key);
    }

    std::string get_or(const std::string& key, const std::string& default_value) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        return default_value;
    }

    bool has(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return true;
            }
        }
        return false;
    }
};

class DensityInitialConditionStyle {
public:
    virtual ~DensityInitialConditionStyle() = default;
    virtual const std::string& type_name() const = 0;

    virtual std::shared_ptr<DensityInitialConditionCommandBase> parse_command(
        int component,
        const InitialConditionArgs& args,
        const Params& params
    ) const = 0;

    virtual std::unique_ptr<DensityInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const DensityInitialConditionCommandBase> command
    ) const = 0;
};

class MomentumInitialConditionStyle {
public:
    virtual ~MomentumInitialConditionStyle() = default;
    virtual const std::string& type_name() const = 0;

    virtual std::shared_ptr<MomentumInitialConditionCommandBase> parse_command(
        int components,
        const InitialConditionArgs& args,
        const Params& params
    ) const = 0;

    virtual std::unique_ptr<MomentumInitialCondition> create_initial_condition(
        const Params& params,
        std::shared_ptr<const MomentumInitialConditionCommandBase> command
    ) const = 0;
};


class InitialConditionRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<DensityInitialConditionStyle>> density_styles;
    std::unordered_map<std::string, std::unique_ptr<MomentumInitialConditionStyle>> momentum_styles;

public:
    void register_density_style(std::unique_ptr<DensityInitialConditionStyle> style) {
        const std::string name = style->type_name();
        auto result = density_styles.emplace(name, std::move(style));
        if (!result.second) {
            throw std::runtime_error("Duplicate density initial condition style registration: " + name);
        }
    }

    void register_momentum_style(std::unique_ptr<MomentumInitialConditionStyle> style) {
        const std::string name = style->type_name();
        auto result = momentum_styles.emplace(name, std::move(style));
        if (!result.second) {
            throw std::runtime_error("Duplicate momentum initial condition style registration: " + name);
        }
    }

    const DensityInitialConditionStyle& get_density(const std::string& type) const {
        auto it = density_styles.find(type);
        if (it == density_styles.end()) {
            throw std::runtime_error("Unknown density initial condition type: " + type);
        }
        return *(it->second);
    }

    const MomentumInitialConditionStyle& get_momentum(const std::string& type) const {
        auto it = momentum_styles.find(type);
        if (it == momentum_styles.end()) {
            throw std::runtime_error("Unknown momentum initial condition type: " + type);
        }
        return *(it->second);
    }
};

#endif
