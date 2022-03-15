from conans import ConanFile

class CTiPRD(ConanFile):
    options = {}
    name = "CTiPRD"
    version = "0.1"
    requires = (
        "pybind11/2.9.1",
        "spdlog/1.9.2",

        "catch2/2.13.7"
    )
    generators = "cmake", "gcc", "txt", "cmake_find_package"

