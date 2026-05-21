#ifndef SFI_MODEL_THERMODYNAMICS_REGISTRY_H
#define SFI_MODEL_THERMODYNAMICS_REGISTRY_H

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Params;
class ThermodynamicsModel;

struct ThermodynamicsModelCommandBase {
    std::string type;
    virtual ~ThermodynamicsModelCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;
};

struct ThermodynamicsModelArgs {
    std::vector<std::pair<std::string, std::string>> entries;

    std::string get_required(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        throw std::runtime_error("Missing required thermodynamics model argument: " + key);
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

class ThermodynamicsModelStyle {
public:
    virtual ~ThermodynamicsModelStyle() = default;
    virtual const std::string& type_name() const = 0;

    virtual std::shared_ptr<ThermodynamicsModelCommandBase> create_default_command(
        const Params& params
    ) const = 0;

    virtual void update_command(
        ThermodynamicsModelCommandBase& command,
        const ThermodynamicsModelArgs& args,
        const Params& params
    ) const = 0;

    virtual std::unique_ptr<ThermodynamicsModel> create_model(
        const Params& params,
        std::shared_ptr<const ThermodynamicsModelCommandBase> command
    ) const = 0;
};


class ThermodynamicsModelRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<ThermodynamicsModelStyle>> thermo_styles;

public:
    void register_thermo_style(std::unique_ptr<ThermodynamicsModelStyle> style) {
        const std::string name = style->type_name();
        auto result = thermo_styles.emplace(name, std::move(style));
        if (!result.second) {
            throw std::runtime_error("Duplicate thermodynamics style registration: " + name);
        }
    };
    const ThermodynamicsModelStyle& get_thermo(const std::string& type) const {
        auto it = thermo_styles.find(type);
        if (it == thermo_styles.end()) {
            throw std::runtime_error("Unknown thermodynamics type: " + type);
        }
        return *(it->second);
    }
};

#endif
