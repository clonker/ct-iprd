from conans import ConanFile


class CTiPRDTests(ConanFile):
    options = {}
    name = "CTiPRDTests"
    version = "0.1"
    requires = (
        "catch2/2.13.7"
    )
    generators = "cmake", "gcc", "txt", "cmake_find_package"
