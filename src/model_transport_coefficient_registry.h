#ifndef SFI_MODEL_TRANSPORT_COEFFICIENT_REGISTRY_H
#define SFI_MODEL_TRANSPORT_COEFFICIENT_REGISTRY_H

#include "model_requirement.h"

#include <memory>
#include <ostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Params;
class TransportCoefficient;

struct TransportCoefficientCommandBase {
    std::string type;
    virtual ~TransportCoefficientCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;
};

struct TransportCoefficientArgs {
    std::vector<std::pair<std::string, std::string>> entries;

    std::string get_required(const std::string& key) const {
        for (const auto& kv : entries) {
            if (kv.first == key) {
                return kv.second;
            }
        }
        throw std::runtime_error("Missing required transport coefficient model argument: " + key);
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

class TransportCoefficientStyle {
public:
    virtual ~TransportCoefficientStyle() = default;
    virtual const std::string& type_name() const = 0;

    virtual void validate_command(const TransportCoefficientCommandBase& command, const TransportCoefficientRequirement& requirement, const Params& params) const = 0;

    virtual std::shared_ptr<TransportCoefficientCommandBase> create_default_command(const Params& params) const = 0;

    virtual void update_command(TransportCoefficientCommandBase& command, const TransportCoefficientArgs& args, const Params& params) const = 0;

    virtual std::unique_ptr<TransportCoefficient> create(const Params& params, std::shared_ptr<const TransportCoefficientCommandBase> command) const = 0;
};


class TransportCoefficientRegistry {
private:
    std::unordered_map<std::string, std::unique_ptr<TransportCoefficientStyle>> transport_styles;

public:
    void register_transport_style(std::unique_ptr<TransportCoefficientStyle> style) {
        const std::string name = style->type_name();
        auto result = transport_styles.emplace(name, std::move(style));
        if (!result.second) {
            throw std::runtime_error("Duplicate Transport Coefficients style registration: " + name);
        }
    };
    const TransportCoefficientStyle& get_transport(const std::string& type) const {
        auto it = transport_styles.find(type);
        if (it == transport_styles.end()) {
            throw std::runtime_error("Unknown Transport Coefficients type: " + type);
        }
        return *(it->second);
    }
};

#endif
