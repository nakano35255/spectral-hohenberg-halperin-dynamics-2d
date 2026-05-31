#ifndef SHHD_MEASURE_REGISTRY_H
#define SHHD_MEASURE_REGISTRY_H

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Params;
class Domain2D;
class Measure;

class Thermodynamics;
class FreeEnergy;
class TransportCoefficient;

struct MeasureCommandBase {
    std::string id;
    std::string type;
    bool enabled = false;

    virtual ~MeasureCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;
};

struct MeasureArgs {
    std::vector<std::pair<std::string, std::string>> entries;

    std::string get_required(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        throw std::runtime_error("Missing required measure argument: " + key);
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

class MeasureStyle {
public:
    virtual ~MeasureStyle() = default;

    virtual const std::string& type_name() const = 0;

    virtual std::shared_ptr<MeasureCommandBase> parse_command(
        const std::string& id,
        bool enabled,
        const MeasureArgs& args,
        const Params& params
    ) const = 0;

    virtual std::unique_ptr<Measure> create_measure(
        const Params& params,
        const Domain2D& domain,
        const Thermodynamics& thermodynamics,
        const FreeEnergy& free_energy,
        const TransportCoefficient& transport_coefficient,
        std::shared_ptr<const MeasureCommandBase> command
    ) const = 0;
};

class MeasureRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<MeasureStyle>> styles;

public:
    void register_style(std::unique_ptr<MeasureStyle> style) {
        const std::string name = style->type_name();
        auto [it, inserted] = styles.emplace(name, std::move(style));
        if (!inserted) {
            throw std::runtime_error("Duplicate measure style registration: " + name);
        }
    }

    const MeasureStyle& get(const std::string& type) const {
        auto it = styles.find(type);
        if (it == styles.end()) {
            throw std::runtime_error("Unknown measure type: " + type);
        }
        return *(it->second);
    }
};

#endif
