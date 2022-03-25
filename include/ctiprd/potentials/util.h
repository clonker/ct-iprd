//
// Created by mho on 3/25/22.
//

namespace ctiprd::potentials {

template<typename dtype, typename... ExternalPotentials>
struct Cutoff;

template<typename dtype>
struct Cutoff<dtype> {
    constexpr static dtype value = 0;
};

}
