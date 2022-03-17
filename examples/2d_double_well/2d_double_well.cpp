//
// Created by mho on 3/15/22.
//

#include <ctiprd/systems/double_well.h>

using System = ctiprd::systems::DoubleWell<float>;

int main() {
    auto pool = ctiprd::config::make_pool(5);
    auto integrator = System::Integrator{pool};
    integrator.particles()->addParticle({{0., 0.}});
    integrator.step(1e-3);
    return 0;
}
