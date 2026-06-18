#include "restart_io.h"

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include <mpi.h>

namespace {
constexpr const char* RESTART_MAGIC = "SHHD_RESTART_V1";
constexpr int RESTART_STATE_TAG = 7301;

struct RestartHeader {
    int step = 0;
    double time = 0.0;
    int nx = -1;
    int ny = -1;
    int nkx = -1;
    int nky = -1;
    int num_order_parameters = -1;
    int num_fields = -1;
};

int expected_num_fields(const Params& params) {
    const int nfields = params.physics.num_order_parameters + 3;
    if (nfields < 3) {
        throw std::runtime_error("RestartIO: invalid number of fields.");
    }
    return nfields;
}

std::size_t double_count_from_complex_count(std::size_t complex_count) {
    if (complex_count > std::numeric_limits<std::size_t>::max() / 2) {
        throw std::runtime_error("RestartIO: payload size overflow.");
    }
    return 2 * complex_count;
}

int checked_int(std::size_t value, const std::string& name) {
    if (value > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
        throw std::runtime_error("RestartIO: " + name + " exceeds MPI int range.");
    }
    return static_cast<int>(value);
}

std::size_t total_count_from_counts_and_displs(
    const std::vector<int>& counts,
    const std::vector<int>& displs
) {
    if (counts.empty()) {
        return 0;
    }
    return static_cast<std::size_t>(displs.back())
         + static_cast<std::size_t>(counts.back());
}

void build_spectral_counts_and_displs(
    const Domain2D& domain,
    int nfields,
    std::vector<int>& counts,
    std::vector<int>& displs
) {
    counts.assign(static_cast<std::size_t>(domain.size()), 0);
    displs.assign(static_cast<std::size_t>(domain.size()), 0);

    std::size_t offset = 0;
    for (int rank = 0; rank < domain.size(); ++rank) {
        const Box2D box = domain.spectral_box_for_rank(rank);
        const std::size_t local_field_size =
            static_cast<std::size_t>(box.size[0])
          * static_cast<std::size_t>(box.size[1]);
        const std::size_t local_complex_count =
            static_cast<std::size_t>(nfields) * local_field_size;
        const std::size_t local_double_count =
            double_count_from_complex_count(local_complex_count);

        counts[static_cast<std::size_t>(rank)] =
            checked_int(local_double_count, "MPI count");
        displs[static_cast<std::size_t>(rank)] =
            checked_int(offset, "MPI displacement");
        offset += local_double_count;
    }
}

std::vector<double> local_state_to_double_pairs(
    const State& state,
    const Domain2D& domain,
    int nfields
) {
    const std::size_t local_complex_count =
        static_cast<std::size_t>(nfields) * domain.spectral_size();
    std::vector<double> local_data(
        double_count_from_complex_count(local_complex_count)
    );

    const Complex* spectral = state.data();
    for (std::size_t i = 0; i < local_complex_count; ++i) {
        local_data[2 * i] = spectral[i].real();
        local_data[2 * i + 1] = spectral[i].imag();
    }
    return local_data;
}

std::vector<double> gathered_to_global_layout(
    const std::vector<double>& gathered,
    const std::vector<int>& displs,
    const Domain2D& domain,
    int nfields
) {
    const int nkx = domain.nx_global() / 2 + 1;
    const int nky = domain.ny_global();
    const std::size_t global_field_size =
        static_cast<std::size_t>(nkx) * static_cast<std::size_t>(nky);

    std::vector<double> global_data(
        double_count_from_complex_count(
            static_cast<std::size_t>(nfields) * global_field_size
        )
    );

    for (int rank = 0; rank < domain.size(); ++rank) {
        const Box2D box = domain.spectral_box_for_rank(rank);
        const int local_nkx = box.size[0];
        const int local_nky = box.size[1];
        const std::size_t local_field_size =
            static_cast<std::size_t>(local_nkx)
          * static_cast<std::size_t>(local_nky);
        const std::size_t rank_offset =
            static_cast<std::size_t>(displs[static_cast<std::size_t>(rank)]);

        for (int field = 0; field < nfields; ++field) {
            for (int ly = 0; ly < local_nky; ++ly) {
                const int ky = box.low[1] + ly;
                for (int lx = 0; lx < local_nkx; ++lx) {
                    const int kx = box.low[0] + lx;
                    const std::size_t local_index =
                        static_cast<std::size_t>(field) * local_field_size
                      + static_cast<std::size_t>(ly) * static_cast<std::size_t>(local_nkx)
                      + static_cast<std::size_t>(lx);
                    const std::size_t global_index =
                        static_cast<std::size_t>(field) * global_field_size
                      + static_cast<std::size_t>(ky) * static_cast<std::size_t>(nkx)
                      + static_cast<std::size_t>(kx);

                    global_data[2 * global_index] =
                        gathered[rank_offset + 2 * local_index];
                    global_data[2 * global_index + 1] =
                        gathered[rank_offset + 2 * local_index + 1];
                }
            }
        }
    }

    return global_data;
}

std::vector<double> global_to_scatter_layout(
    const std::vector<double>& global_data,
    const Domain2D& domain,
    int nfields,
    const std::vector<int>& counts,
    const std::vector<int>& displs
) {
    const int nkx = domain.nx_global() / 2 + 1;
    const std::size_t global_field_size =
        static_cast<std::size_t>(nkx)
      * static_cast<std::size_t>(domain.ny_global());

    std::vector<double> scatter_data(
        total_count_from_counts_and_displs(counts, displs)
    );

    for (int rank = 0; rank < domain.size(); ++rank) {
        const Box2D box = domain.spectral_box_for_rank(rank);
        const int local_nkx = box.size[0];
        const int local_nky = box.size[1];
        const std::size_t local_field_size =
            static_cast<std::size_t>(local_nkx)
          * static_cast<std::size_t>(local_nky);
        const std::size_t rank_offset =
            static_cast<std::size_t>(displs[static_cast<std::size_t>(rank)]);

        for (int field = 0; field < nfields; ++field) {
            for (int ly = 0; ly < local_nky; ++ly) {
                const int ky = box.low[1] + ly;
                for (int lx = 0; lx < local_nkx; ++lx) {
                    const int kx = box.low[0] + lx;
                    const std::size_t local_index =
                        static_cast<std::size_t>(field) * local_field_size
                      + static_cast<std::size_t>(ly) * static_cast<std::size_t>(local_nkx)
                      + static_cast<std::size_t>(lx);
                    const std::size_t global_index =
                        static_cast<std::size_t>(field) * global_field_size
                      + static_cast<std::size_t>(ky) * static_cast<std::size_t>(nkx)
                      + static_cast<std::size_t>(kx);

                    scatter_data[rank_offset + 2 * local_index] =
                        global_data[2 * global_index];
                    scatter_data[rank_offset + 2 * local_index + 1] =
                        global_data[2 * global_index + 1];
                }
            }
        }
    }

    return scatter_data;
}

void write_header(
    std::ofstream& out,
    const Params& params,
    const Domain2D& domain,
    int nfields,
    int step,
    double time
) {
    out << RESTART_MAGIC << '\n';
    out << "step " << step << '\n';
    out << "time " << std::scientific << std::setprecision(17) << time << '\n';
    out << "nx " << domain.nx_global() << '\n';
    out << "ny " << domain.ny_global() << '\n';
    out << "nkx " << domain.nx_global() / 2 + 1 << '\n';
    out << "nky " << domain.ny_global() << '\n';
    out << "num_order_parameters " << params.physics.num_order_parameters << '\n';
    out << "num_fields " << nfields << '\n';
    out << "precision text_float64" << '\n';
    out << "layout field_major_ky_kx_text" << '\n';
    out << "columns field kx_index ky_index real imag" << '\n';
    out << "data" << '\n';
}

RestartHeader read_header(std::ifstream& in) {
    RestartHeader header;

    std::string token;
    if (!(in >> token) || token != RESTART_MAGIC) {
        throw std::runtime_error("RestartIO: invalid restart magic.");
    }

    bool found_data = false;
    while (in >> token) {
        if (token == "data") {
            found_data = true;
            break;
        }
        if (token == "step") {
            in >> header.step;
        } else if (token == "time") {
            in >> header.time;
        } else if (token == "nx") {
            in >> header.nx;
        } else if (token == "ny") {
            in >> header.ny;
        } else if (token == "nkx") {
            in >> header.nkx;
        } else if (token == "nky") {
            in >> header.nky;
        } else if (token == "num_order_parameters") {
            in >> header.num_order_parameters;
        } else if (token == "num_fields") {
            in >> header.num_fields;
        } else if (token == "precision") {
            std::string value;
            in >> value;
            if (value != "text_float64") {
                throw std::runtime_error(
                    "RestartIO: unsupported restart precision: " + value
                );
            }
        } else if (token == "layout") {
            std::string value;
            in >> value;
            if (value != "field_major_ky_kx_text") {
                throw std::runtime_error(
                    "RestartIO: unsupported restart layout: " + value
                );
            }
        } else if (token == "columns") {
            std::string field;
            std::string kx;
            std::string ky;
            std::string real;
            std::string imag;
            in >> field >> kx >> ky >> real >> imag;
        } else {
            throw std::runtime_error("RestartIO: unknown restart header key: " + token);
        }

        if (!in) {
            throw std::runtime_error("RestartIO: malformed restart header.");
        }
    }

    if (!found_data) {
        throw std::runtime_error("RestartIO: restart file has no data marker.");
    }

    return header;
}

void validate_header(
    const RestartHeader& header,
    const Params& params,
    const Domain2D& domain
) {
    if (header.nx != domain.nx_global() || header.ny != domain.ny_global()) {
        throw std::runtime_error("RestartIO: grid size mismatch.");
    }
    if (header.nkx != domain.nx_global() / 2 + 1
        || header.nky != domain.ny_global()) {
        throw std::runtime_error("RestartIO: spectral grid size mismatch.");
    }

    const int nfields = expected_num_fields(params);
    if (header.num_fields != nfields) {
        throw std::runtime_error("RestartIO: state field count mismatch.");
    }
}

void write_global_text_payload(
    std::ofstream& out,
    const std::vector<double>& global_data,
    const Domain2D& domain,
    int nfields
) {
    const int nkx = domain.nx_global() / 2 + 1;
    const int nky = domain.ny_global();
    const std::size_t global_field_size =
        static_cast<std::size_t>(nkx) * static_cast<std::size_t>(nky);

    out << std::scientific << std::setprecision(17);
    for (int field = 0; field < nfields; ++field) {
        for (int ky = 0; ky < nky; ++ky) {
            for (int kx = 0; kx < nkx; ++kx) {
                const std::size_t global_index =
                    static_cast<std::size_t>(field) * global_field_size
                  + static_cast<std::size_t>(ky) * static_cast<std::size_t>(nkx)
                  + static_cast<std::size_t>(kx);
                out << field << ' '
                    << kx << ' '
                    << ky << ' '
                    << global_data[2 * global_index] << ' '
                    << global_data[2 * global_index + 1] << '\n';
            }
        }
    }
}

std::vector<double> read_global_text_payload(
    std::ifstream& in,
    const Domain2D& domain,
    int nfields
) {
    const int nkx = domain.nx_global() / 2 + 1;
    const int nky = domain.ny_global();
    const std::size_t global_field_size =
        static_cast<std::size_t>(nkx) * static_cast<std::size_t>(nky);

    std::vector<double> global_data(
        double_count_from_complex_count(
            static_cast<std::size_t>(nfields) * global_field_size
        )
    );

    for (int expected_field = 0; expected_field < nfields; ++expected_field) {
        for (int expected_ky = 0; expected_ky < nky; ++expected_ky) {
            for (int expected_kx = 0; expected_kx < nkx; ++expected_kx) {
                int field = -1;
                int kx = -1;
                int ky = -1;
                double real = 0.0;
                double imag = 0.0;

                if (!(in >> field >> kx >> ky >> real >> imag)) {
                    throw std::runtime_error("RestartIO: restart payload is truncated.");
                }
                if (field != expected_field || kx != expected_kx || ky != expected_ky) {
                    throw std::runtime_error("RestartIO: restart payload order mismatch.");
                }

                const std::size_t global_index =
                    static_cast<std::size_t>(field) * global_field_size
                  + static_cast<std::size_t>(ky) * static_cast<std::size_t>(nkx)
                  + static_cast<std::size_t>(kx);
                global_data[2 * global_index] = real;
                global_data[2 * global_index + 1] = imag;
            }
        }
    }

    return global_data;
}

void broadcast_status_or_throw(
    const Domain2D& domain,
    int ok,
    std::string message,
    const std::string& fallback
) {
    MPI_Bcast(&ok, 1, MPI_INT, 0, domain.comm());

    int length = 0;
    if (domain.rank() == 0) {
        length = checked_int(message.size(), "error message size");
    }
    MPI_Bcast(&length, 1, MPI_INT, 0, domain.comm());

    if (domain.rank() != 0) {
        message.resize(static_cast<std::size_t>(length));
    }
    if (length > 0) {
        MPI_Bcast(&message[0], length, MPI_CHAR, 0, domain.comm());
    }

    if (ok == 0) {
        throw std::runtime_error(message.empty() ? fallback : message);
    }
}
} // namespace

void write_restart_file(
    const std::string& file,
    const Params& params,
    const Domain2D& domain,
    const State& state,
    int step,
    double time
) {
    const int nfields = expected_num_fields(params);
    const std::vector<double> local_data =
        local_state_to_double_pairs(state, domain, nfields);
    const int local_count = checked_int(local_data.size(), "local restart payload size");

    std::vector<int> counts;
    std::vector<int> displs;
    std::vector<double> gathered;

    if (domain.rank() == 0) {
        build_spectral_counts_and_displs(domain, nfields, counts, displs);
        gathered.resize(total_count_from_counts_and_displs(counts, displs));
    }

    if (domain.rank() == 0) {
        std::copy(
            local_data.begin(),
            local_data.end(),
            gathered.begin() + displs[0]
        );
        for (int rank = 1; rank < domain.size(); ++rank) {
            MPI_Recv(
                gathered.data() + displs[static_cast<std::size_t>(rank)],
                counts[static_cast<std::size_t>(rank)],
                MPI_DOUBLE,
                rank,
                RESTART_STATE_TAG,
                domain.comm(),
                MPI_STATUS_IGNORE
            );
        }
    } else {
        MPI_Send(
            local_data.data(),
            local_count,
            MPI_DOUBLE,
            0,
            RESTART_STATE_TAG,
            domain.comm()
        );
    }

    int write_ok = 1;
    std::string error_message;

    if (domain.rank() == 0) {
        try {
            const std::vector<double> global_data =
                gathered_to_global_layout(gathered, displs, domain, nfields);

            std::ofstream out(file, std::ios::out | std::ios::trunc);
            if (!out) {
                throw std::runtime_error("RestartIO: cannot open restart file: " + file);
            }

            write_header(out, params, domain, nfields, step, time);
            write_global_text_payload(out, global_data, domain, nfields);

            if (!out) {
                throw std::runtime_error("RestartIO: failed to write restart file: " + file);
            }
        } catch (const std::exception& e) {
            write_ok = 0;
            error_message = e.what();
        }
    }

    broadcast_status_or_throw(
        domain,
        write_ok,
        error_message,
        "RestartIO: failed to write restart file."
    );
}

void read_restart_file(
    const std::string& file,
    const Params& params,
    const Domain2D& domain,
    State& state,
    int& step,
    double& time
) {
    const int nfields = expected_num_fields(params);

    RestartHeader header;
    std::vector<double> global_data;

    int read_ok = 1;
    std::string error_message;

    if (domain.rank() == 0) {
        try {
            std::ifstream in(file);
            if (!in) {
                throw std::runtime_error("RestartIO: cannot open restart file: " + file);
            }

            header = read_header(in);
            validate_header(header, params, domain);
            global_data = read_global_text_payload(in, domain, nfields);
        } catch (const std::exception& e) {
            read_ok = 0;
            error_message = e.what();
        }
    }

    broadcast_status_or_throw(
        domain,
        read_ok,
        error_message,
        "RestartIO: failed to read restart file."
    );

    if (domain.rank() == 0) {
        step = header.step;
        time = header.time;
    }
    MPI_Bcast(&step, 1, MPI_INT, 0, domain.comm());
    MPI_Bcast(&time, 1, MPI_DOUBLE, 0, domain.comm());

    std::vector<int> counts;
    std::vector<int> displs;
    std::vector<double> scatter_data;

    if (domain.rank() == 0) {
        build_spectral_counts_and_displs(domain, nfields, counts, displs);
        scatter_data =
            global_to_scatter_layout(global_data, domain, nfields, counts, displs);
    }

    const std::size_t local_complex_count =
        static_cast<std::size_t>(nfields) * domain.spectral_size();
    std::vector<double> local_data(
        double_count_from_complex_count(local_complex_count)
    );
    const int local_count = checked_int(local_data.size(), "local restart payload size");

    if (domain.rank() == 0) {
        std::copy(
            scatter_data.begin() + displs[0],
            scatter_data.begin() + displs[0] + counts[0],
            local_data.begin()
        );
        for (int rank = 1; rank < domain.size(); ++rank) {
            MPI_Send(
                scatter_data.data() + displs[static_cast<std::size_t>(rank)],
                counts[static_cast<std::size_t>(rank)],
                MPI_DOUBLE,
                rank,
                RESTART_STATE_TAG,
                domain.comm()
            );
        }
    } else {
        MPI_Recv(
            local_data.data(),
            local_count,
            MPI_DOUBLE,
            0,
            RESTART_STATE_TAG,
            domain.comm(),
            MPI_STATUS_IGNORE
        );
    }

    Complex* spectral = state.data();
    for (std::size_t i = 0; i < local_complex_count; ++i) {
        spectral[i] = Complex(local_data[2 * i], local_data[2 * i + 1]);
    }
}
