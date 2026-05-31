#ifndef SHHD_FCALCULATOR_FIELD_REQUESTS_AND_LAYOUT_H
#define SHHD_FCALCULATOR_FIELD_REQUESTS_AND_LAYOUT_H

#include "fcalculator_dynamics_mode.h"
#include "model_free_energy.h"
#include "model_thermodynamics.h"
#include "model_transport_coefficient.h"
#include "simulationinfo.h"

#include <stdexcept>

struct PhysicalFieldRequest {
     bool need_psi = false;
     bool need_rho = false;
     bool need_j = false;
     bool need_velocity = false;

     bool any() const {
          return need_psi || need_rho || need_j || need_velocity;
     }

     void include(const PhysicalFieldRequest& other) {
          need_psi = need_psi || other.need_psi;
          need_rho = need_rho || other.need_rho;
          need_j = need_j || other.need_j;
          need_velocity = need_velocity || other.need_velocity;
     }
};

struct PhysicalRHSFieldRequests {
     PhysicalFieldRequest psi_det;
     PhysicalFieldRequest j_det;

     static PhysicalRHSFieldRequests build(
          const Params& params,
          DynamicsMode dynamics_mode,
          const FreeEnergy& free_energy,
          const Thermodynamics& thermodynamics,
          const TransportCoefficient& transport_coefficient
     ) {
          PhysicalRHSFieldRequests requests;
          requests.psi_det = build_psi_det_request(params, free_energy);
          requests.j_det = build_j_det_request(params, dynamics_mode, thermodynamics, transport_coefficient);
          return requests;
     }

     PhysicalFieldRequest total() const {
          PhysicalFieldRequest request;
          request.include(psi_det);
          request.include(j_det);
          return request;
     }

private:
     static PhysicalFieldRequest build_psi_det_request(const Params& params, const FreeEnergy& free_energy) {
          PhysicalFieldRequest request;

          if (params.physics.num_order_parameters <= 0) {
               return request;
          }

          if (free_energy.has_physical_chemical_potential()) {
               request.need_psi = true;
          }

          if (params.fix.order_parameter_advection) {
               request.need_psi = true;
               request.need_j = true;
               request.need_velocity = true;
          }

          return request;
     }

     static PhysicalFieldRequest build_j_det_request(
          const Params& params,
          DynamicsMode dynamics_mode,
          const Thermodynamics& thermodynamics,
          const TransportCoefficient& transport_coefficient
     ) {
          PhysicalFieldRequest request;

          if (is_quiescent_mode(dynamics_mode)) {
               return request;
          }

          const bool needs_physical_viscosity =
               transport_coefficient.has_physical_shear_viscosity() || transport_coefficient.has_physical_bulk_viscosity();

          if (needs_physical_viscosity) {
               throw std::runtime_error("physical viscosity is not implemented yet.");
          }

          const bool has_constant_viscosity =
               transport_coefficient.shear_viscosity() != 0.0 || transport_coefficient.bulk_viscosity() != 0.0;

          if (is_compressible_mode(dynamics_mode) && has_constant_viscosity) {
               request.need_j = true;
               request.need_velocity = true;
          }

          if (params.fix.momentum_advection) {
               request.need_j = true;
               request.need_velocity = true;
          }

          if (thermodynamics.has_physical_pressure()) {
               if (!is_compressible_mode(dynamics_mode)) {
                    throw std::runtime_error("physical pressure requires compressible time evolution.");
               }
               request.need_rho = true;
          }

          return request;
     }
};

struct PhysicalFieldOffsets {
     int psi = -1;
     int rho = -1;
     int jx = -1;
     int jy = -1;
     int vx = -1;
     int vy = -1;
};

struct PhysicalFieldLayout {
     bool transform_psi = false;
     bool transform_rho = false;
     bool transform_j = false;
     bool make_velocity = false;

     PhysicalFieldOffsets offsets;

     int ntransform = 0;
     int nphysical = 0;

     static PhysicalFieldLayout from_requests(DynamicsMode dynamics_mode, const PhysicalRHSFieldRequests& requests, int num_order_parameters) {
          return from_request(dynamics_mode, requests.total(), num_order_parameters);
     }

     static PhysicalFieldLayout from_request(DynamicsMode dynamics_mode, const PhysicalFieldRequest& request, int num_order_parameters) {
          validate_request(dynamics_mode, request, num_order_parameters);

          PhysicalFieldLayout layout;
          layout.transform_psi = request.need_psi;
          layout.transform_j = request.need_j;
          layout.make_velocity = request.need_velocity;
          layout.transform_rho =
               request.need_rho || (is_compressible_mode(dynamics_mode) && request.need_velocity);

          layout.assign_offsets(num_order_parameters);
          return layout;
     }

private:
     static void validate_request(
          DynamicsMode dynamics_mode,
          const PhysicalFieldRequest& request,
          int num_order_parameters
     ) {
          if (num_order_parameters < 0) {
               throw std::runtime_error("PhysicalFieldLayout requires nonnegative num_order_parameters.");
          }

          if (request.need_psi && num_order_parameters <= 0) {
               throw std::runtime_error("physical psi requested with no order parameters.");
          }

          if (is_quiescent_mode(dynamics_mode) &&
               (request.need_rho || request.need_j || request.need_velocity)) {
               throw std::runtime_error("quiescent workspace cannot request hydrodynamic physical fields.");
          }

          if (request.need_velocity && !request.need_j) {
               throw std::runtime_error("velocity request requires physical momentum.");
          }

          if (!is_compressible_mode(dynamics_mode) && request.need_rho) {
               throw std::runtime_error("density IFFT is only supported for compressible time evolution.");
          }
     }

     void assign_offsets(int num_order_parameters) {
          if (transform_psi) {
               offsets.psi = ntransform;
               ntransform += num_order_parameters;
          }

          if (transform_rho) {
               offsets.rho = ntransform;
               ntransform += 1;
          }

          if (transform_j) {
               offsets.jx = ntransform;
               offsets.jy = ntransform + 1;
               ntransform += 2;
          }

          nphysical = ntransform;

          if (make_velocity) {
               offsets.vx = nphysical;
               offsets.vy = nphysical + 1;
               nphysical += 2;
          }
     }
};

#endif