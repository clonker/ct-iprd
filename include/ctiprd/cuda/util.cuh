#pragma once

#include <cuda.h>
#include <utility>

namespace ctiprd::cuda {

struct LaunchConfig {

    LaunchConfig() = default;

    int gridSize {0};
    int blockSize {0};
    unsigned int sharedSize {0};
    CUstream stream {nullptr};
};

template<typename F, typename... Args>
void launch(F kernel, LaunchConfig cfg, Args&&... args) {
    kernel<<<cfg.gridSize, cfg.blockSize, cfg.sharedSize, cfg.stream>>>(std::forward<Args>(args)...);
}

}
// CUDAHOSTCXX=/home/mho/miniconda3/envs/qmlw/bin/x86_64-conda-linux-gnu-g++ ??
