#ifndef SFI_MODEL_FREE_ENERGY_REGISTRY_H
#define SFI_MODEL_FREE_ENERGY_REGISTRY_H

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Params;
class FreeEnergy;

struct FreeEnergyCommandBase {
    std::string type;

    virtual ~FreeEnergyCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;
};

struct FreeEnergyArgs {
    std::vector<std::pair<std::string, std::string>> entries;

    bool has(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return true;
            }
        }
        return false;
    }

    std::string get_required(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        throw std::runtime_error("Missing required free energy model argument: " + key);
    }

    std::string get_or(const std::string& key, const std::string& default_value) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        return default_value;
    }
};

class FreeEnergyStyle {
public:
    virtual ~FreeEnergyStyle() = default;
    virtual const std::string& type_name() const = 0;

    virtual std::shared_ptr<FreeEnergyCommandBase> create_default_command(const Params& params) const = 0;

    virtual void update_command(FreeEnergyCommandBase& command, const FreeEnergyArgs& args, const Params& params) const = 0;

    virtual std::unique_ptr<FreeEnergy> create(const Params& params, std::shared_ptr<const FreeEnergyCommandBase> command) const = 0;
};

class FreeEnergyRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<FreeEnergyStyle>> free_energy_styles;

public:
    void register_free_energy_style(std::unique_ptr<FreeEnergyStyle> style) {
        const std::string name = style->type_name();
        auto result = free_energy_styles.emplace(name, std::move(style));
        if (!result.second) {
            throw std::runtime_error("Duplicate free energy style registration: " + name);
        }
    }

    const FreeEnergyStyle& get_free_energy(const std::string& type) const {
        auto it = free_energy_styles.find(type);
        if (it == free_energy_styles.end()) {
            throw std::runtime_error("Unknown free energy type: " + type);
        }
        return *(it->second);
    }
};

#endif
