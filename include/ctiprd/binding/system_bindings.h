//
// Created by mho on 4/6/22.
//

#pragma once

#include <cstddef>
#include <vector>
#include <string_view>

#include <nanobind/nanobind.h>
#include <nanobind/tensor.h>

namespace ctiprd::binding {

namespace nb = nanobind;

template<typename dtype, std::size_t... shape>
using np_array = nb::tensor<dtype, nb::shape<shape...>, nb::c_contig, nb::device::cpu>;

template<typename dtype>
void exportBaseTypes(nb::module_ &module) {
    nb::class_<ParticleType<dtype>>(module, "ParticleType")
        .template def_readonly("name", &ParticleType<dtype>::name)
        .template def_readonly("diffusion_constant", &ParticleType<dtype>::diffusionConstant);
}

template<typename System>
void exportSystem(nb::module_ &module, std::string_view name) {
    std::string strName {name};
    auto clazz = nb::class_<System>(module, strName.c_str());
    clazz.template def_readonly_static("periodic", &System::periodic);
    clazz.template def_readonly_static("dim", &System::DIM);
    clazz.template def_property_readonly("box_size", [](const System &self) {
        np_array<typename System::dtype, System::DIM> out {};
        std::copy(begin(self.boxSize), end(self.boxSize), reinterpret_cast<typename System::dtype*>(out.data()));
        return out;
    });
    clazz.template def_readonly_static("particle_types", &System::types);

}

}

