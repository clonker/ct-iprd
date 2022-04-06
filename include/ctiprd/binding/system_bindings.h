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

template<typename dtype>
using np_array = nb::tensor<dtype, nb::c_contig>;

template<typename System, typename Module>
void exportSystem(Module &module, std::string_view name) {

    auto clazz = nb::class_<System>(module, std::string(name));
    clazz.template def_readonly_static("periodic", &System::periodic);
    clazz.template def_readonly_static("dim", &System::DIM);
    clazz.template def_property_readonly("box_size", [](const System &self) {
        np_array<typename System::dtype> out (std::vector<std::size_t>{System::DIM});
        std::copy(begin(self.DIM), end(self.DIM), out.template mutable_data(0));
        return out;
    });
}

}

