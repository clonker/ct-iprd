//
// Created by mho on 3/15/22.
//

#include <ctiprd/config.h>
#include <ctiprd/systems/double_well.h>
#include <ctiprd/cpu/integrators/EulerMaruyama.h>

using System = ctiprd::systems::DoubleWell<float>;

int main() {
    auto pool = ctiprd::config::make_pool(5);
    System system {};
    auto integrator = ctiprd::cpu::integrator::EulerMaruyama<System>{system, pool};
    integrator.particles()->addParticle({{0., 0.}}, "A");
    integrator.step(1e-3);
    return 0;
}
