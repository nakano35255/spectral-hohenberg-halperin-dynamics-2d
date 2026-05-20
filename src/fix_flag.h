#ifndef SFI_FIX_FLAG_H
#define SFI_FIX_FLAG_H

#include <cstdint>
#include <stdexcept>
#include <string>

enum class FixFlag : std::uint32_t {
     Noise          = 1u << 0,
     Shear          = 1u << 1,
     Nonlinear      = 1u << 2,
     Barodiffusion  = 1u << 3,
     Quiescent      = 1u << 4,
     Incompressible = 1u << 5
};

enum class FixArgKind {
     None,
     Seed,
     Rate
};

struct FixSpec {
     const char* name;
     FixFlag flag;
     FixArgKind arg_kind;
};

inline const FixSpec FIX_SPECS[] = {
     {"noise", FixFlag::Noise, FixArgKind::Seed},
     {"shear", FixFlag::Shear, FixArgKind::Rate},
     {"nonlinear", FixFlag::Nonlinear, FixArgKind::None},
     {"barodiffusion", FixFlag::Barodiffusion, FixArgKind::None},
     {"quiescent", FixFlag::Quiescent, FixArgKind::None},
     {"incompressible", FixFlag::Incompressible, FixArgKind::None},
};

inline const FixSpec& find_fix_spec(const std::string& name) {
     for (const auto& spec : FIX_SPECS) {
          if (name == spec.name) {
               return spec;
          }
     }
     throw std::runtime_error("unknown fix style: " + name);
}

inline std::uint32_t fix_bit(FixFlag flag) {
     return static_cast<std::uint32_t>(flag);
}

#endif
