//
// Created by mho on 4/6/22.
//

#pragma once

#include <cstddef>
#include <vector>
#include <string_view>

#include <nanobind/nanobind.h>
#include <nanobind/tensor.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/tuple.h>

namespace ctiprd::binding {

namespace nb = nanobind;

template<typename dtype, std::size_t... shape>
using np_array = nb::tensor<dtype, nb::shape<shape...>, nb::numpy, nb::c_contig, nb::device::cpu>;

template<typename dtype>
void exportBaseTypes(nb::module_ &module) {
    nb::class_<ParticleType<dtype>>(module, "ParticleType")
        .template def_readonly("name", &ParticleType<dtype>::name)
        .template def_readonly("diffusion_constant", &ParticleType<dtype>::diffusionConstant);

    nb::class_<reactions::doi::Catalysis<dtype>>(module, "Catalysis")
            .template def_readwrite("catalyst", &reactions::doi::Catalysis<dtype>::catalyst)
            .template def_readwrite("eductType", &reactions::doi::Catalysis<dtype>::eductType)
            .template def_readwrite("productType", &reactions::doi::Catalysis<dtype>::productType)
            .template def_readwrite("reactionRadius", &reactions::doi::Catalysis<dtype>::reactionRadius)
            .template def_readwrite("rate", &reactions::doi::Catalysis<dtype>::rate);

    nb::class_<reactions::doi::Decay<dtype>>(module, "Decay")
            .template def_readwrite("eductType", &reactions::doi::Decay<dtype>::eductType)
            .template def_readwrite("rate", &reactions::doi::Decay<dtype>::rate);

    nb::class_<reactions::doi::Fusion<dtype>>(module, "Fusion")
            .template def_readwrite("eductType1", &reactions::doi::Fusion<dtype>::eductType1)
            .template def_readwrite("eductType2", &reactions::doi::Fusion<dtype>::eductType2)
            .template def_readwrite("productType", &reactions::doi::Fusion<dtype>::productType)
            .template def_readwrite("reactionRadius", &reactions::doi::Fusion<dtype>::reactionRadius)
            .template def_readwrite("rate", &reactions::doi::Fusion<dtype>::rate)
            .template def_readwrite("w1", &reactions::doi::Fusion<dtype>::w1)
            .template def_readwrite("w2", &reactions::doi::Fusion<dtype>::w2);

    nb::class_<reactions::doi::Fission<dtype>>(module, "Fission")
            .template def_readwrite("educt_type", &reactions::doi::Fission<dtype>::eductType)
            .template def_readwrite("product_type_1", &reactions::doi::Fission<dtype>::productType1)
            .template def_readwrite("product_type_2", &reactions::doi::Fission<dtype>::productType2)
            .template def_readwrite("product_distance", &reactions::doi::Fission<dtype>::productDistance)
            .template def_readwrite("rate", &reactions::doi::Fission<dtype>::rate);

}

template<typename System>
void exportSystem(nb::module_ &module, std::string_view name) {
    std::string strName {name};
    auto clazz = nb::class_<System>(module, strName.c_str()).def(nb::init<>());
    clazz.template def_readonly_static("periodic", &System::periodic);
    clazz.template def_readonly_static("dim", &System::DIM);
    clazz.template def_property_readonly("box_size", [](const System &self) {
        std::vector<typename System::dtype> out {begin(self.boxSize), end(self.boxSize)};
        return out;
    });
    clazz.template def_property_readonly("particle_types", [](const System &self) {
        return std::vector {begin(self.types), end(self.types)};
    });
    clazz.template def_property_readonly("reactions_o1", [](const System &self) {
        return self.reactionsO1;
    });
    clazz.template def_property_readonly("reactions_o2", [](const System &self) {
        return self.reactionsO2;
    });
}

}

