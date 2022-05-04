//
// Created by mho on 5/4/22.
//

#pragma once

#include <ctiprd/config.h>

namespace ctiprd::binding {
namespace py = pybind11;

template<typename System>
class Trajectory {
    using dtype = typename System::dtype;

public:
    template<typename dtype>
    using np_array = py::array_t<dtype, py::array::c_style>;

    const auto &time() const {
        return time_;
    }

    const auto &positions() const {
        return positions_;
    }

    const auto &types() const {
        return particleTypes_;
    }

    template<typename Integrator, typename Pool>
    void record(std::size_t step, const Integrator &integrator, config::PoolPtr<Pool> pool) {
        py::gil_scoped_acquire acquire;

        const auto nParticles = integrator.particles()->size();
        std::array shapeTraj {nParticles, static_cast<std::size_t>(2)};
        std::array shapeTypes {nParticles};
        positions_.emplace_back(shapeTraj);
        particleTypes_.emplace_back(shapeTypes);
        time_.emplace_back(step);

        auto &traj = positions_.back();
        auto &types = particleTypes_.back();
        auto* trajBegin = traj.mutable_data(0);
        auto* trajEnd = trajBegin + traj.size();
        std::fill(trajBegin, trajEnd, std::numeric_limits<dtype>::infinity());
        auto* typesBegin = types.mutable_data(0);
        auto* typesEnd = typesBegin + types.size();
        std::fill(typesBegin, typesEnd, System::types.size()); // out of bounds type (invalid)

        auto futures = integrator.particles()->forEachParticle([&trajBegin, &typesBegin](auto particleId, const auto &pos, const auto& type, const auto &) {
            typesBegin[particleId] = type;
            trajBegin[2*particleId] = pos[0];
            trajBegin[2*particleId+1] = pos[1];
        }, pool);
        for(auto &future : futures) {
            future.wait();
        }
    }

private:
    std::vector<std::size_t> time_;
    std::vector<np_array<dtype>> positions_;
    std::vector<np_array<std::size_t>> particleTypes_;
};

}
