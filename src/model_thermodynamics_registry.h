#ifndef SFI_MODEL_THERMODYNAMICS_REGISTRY_H
#define SFI_MODEL_THERMODYNAMICS_REGISTRY_H

#include <iomanip>
#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Params;
class Thermodynamics;

struct ThermodynamicsCommandBase {
    std::string type;
    double sound_speed = 0.0;

    virtual ~ThermodynamicsCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;

    void print_common(std::ostream& os) const {
        os << "  "
           << std::left << std::setw(25)
           << "Thermodynamics" << ": " << type
           << " cT " << sound_speed
           << '\n';
    }
};

struct ThermodynamicsArgs {
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

class ThermodynamicsStyle {
protected:
    bool update_common_command(ThermodynamicsCommandBase& command, const std::string& key, const std::string& value) const {
        if (key == "cT" || key == "sound_speed") {
            const double sound_speed = std::stod(value);
            if (sound_speed < 0.0) {
                throw std::runtime_error("model thermo " + command.type + " requires nonnegative " + key + ".");
            }
            command.sound_speed = sound_speed;
            return true;
        }
        return false;
    }

public:
    virtual ~ThermodynamicsStyle() = default;
    virtual const std::string& type_name() const = 0;

    virtual std::shared_ptr<ThermodynamicsCommandBase> create_default_command(const Params& params) const = 0;

    virtual void update_command(ThermodynamicsCommandBase& command, const ThermodynamicsArgs& args, const Params& params) const = 0;

    virtual std::unique_ptr<Thermodynamics> create(const Params& params, std::shared_ptr<const ThermodynamicsCommandBase> command) const = 0;
};


class ThermodynamicsRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<ThermodynamicsStyle>> thermo_styles;

public:
    void register_thermo_style(std::unique_ptr<ThermodynamicsStyle> style) {
        const std::string name = style->type_name();
        auto result = thermo_styles.emplace(name, std::move(style));
        if (!result.second) {
            throw std::runtime_error("Duplicate thermodynamics style registration: " + name);
        }
    };
    const ThermodynamicsStyle& get_thermo(const std::string& type) const {
        auto it = thermo_styles.find(type);
        if (it == thermo_styles.end()) {
            throw std::runtime_error("Unknown thermodynamics type: " + type);
        }
        return *(it->second);
    }
};

#endif
