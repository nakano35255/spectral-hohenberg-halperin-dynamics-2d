#ifndef SHHD_MODEL_TRANSPORT_COEFFICIENT_REGISTRY_H
#define SHHD_MODEL_TRANSPORT_COEFFICIENT_REGISTRY_H

#include <cctype>
#include <cmath>
#include <cstddef>
#include <iomanip>
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
    int num_order_parameters = 0;
    double shear_viscosity = 0.0;
    double bulk_viscosity = 0.0;
    std::vector<double> order_parameter_mobility;

    virtual ~TransportCoefficientCommandBase() = default;
    virtual void print(std::ostream& os) const = 0;

    void print_common(std::ostream& os) const {
        os << "  "
           << std::left << std::setw(25)
           << "Transport Coefficients" << ": " << type
           << " eta " << shear_viscosity
           << " zeta " << bulk_viscosity;

        for (int q = 0; q < num_order_parameters; ++q) {
            os << " M[" << q << ',' << q << "] "
               << order_parameter_mobility[static_cast<std::size_t>(q)];
        }

        os << '\n';
    }
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
protected:
    static bool parse_mobility_key(const std::string& key, int& row, int& col) {
        if (key.size() < 6 || key[0] != 'M' || key[1] != '[' || key.back() != ']') {
            return false;
        }

        std::size_t comma = key.find(',', 2);
        if (comma == std::string::npos || comma + 1 >= key.size() - 1) {
            return false;
        }

        for (std::size_t i = 2; i < comma; ++i) {
            if (!std::isdigit(static_cast<unsigned char>(key[i]))) {
                return false;
            }
        }
        for (std::size_t i = comma + 1; i < key.size() - 1; ++i) {
            if (!std::isdigit(static_cast<unsigned char>(key[i]))) {
                return false;
            }
        }

        row = std::stoi(key.substr(2, comma - 2));
        col = std::stoi(key.substr(comma + 1, key.size() - comma - 2));
        return true;
    }

    bool update_common_command(TransportCoefficientCommandBase& command, const std::string& key, const std::string& raw_value) const {
        if (key == "eta" || key == "shear_viscosity") {
            const double value = std::stod(raw_value);
            if (value < 0.0) {
                throw std::runtime_error("model transport " + command.type + " requires nonnegative " + key + ".");
            }
            command.shear_viscosity = value;
            return true;
        }

        if (key == "zeta" || key == "bulk_viscosity") {
            const double value = std::stod(raw_value);
            if (value < 0.0) {
                throw std::runtime_error("model transport " + command.type + " requires nonnegative " + key + ".");
            }
            command.bulk_viscosity = value;
            return true;
        }

        int row = -1;
        int col = -1;
        if (parse_mobility_key(key, row, col)) {
            const double value = std::stod(raw_value);
            if (row < 0 || row >= command.num_order_parameters || col < 0 || col >= command.num_order_parameters) {
                throw std::runtime_error("model transport " + command.type + " mobility index out of range: " + key);
            }

            if (row != col) {
                if (std::abs(value) > 0.0) {
                    throw std::runtime_error("model transport " + command.type + " only supports diagonal mobility; " + key + " must be zero.");
                }
                return true;
            }

            if (value < 0.0) {
                throw std::runtime_error("model transport " + command.type + " requires nonnegative " + key + ".");
            }
            command.order_parameter_mobility[static_cast<std::size_t>(row)] = value;
            return true;
        }

        return false;
    }

public:
    virtual ~TransportCoefficientStyle() = default;
    virtual const std::string& type_name() const = 0;

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
