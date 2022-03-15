//
// Created by mho on 3/15/22.
//

#include "ctiprd/ParticleCollection.h"

using Collection = ctiprd::ParticleCollection<2, float>;

int main() {
    auto pool = ctiprd::config::make_pool(5);
    return 0;
}
