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

    def configure(self):
        self.options["spdlog"].header_only = True
        self.options["fmt"].header_only = True
        super().configure()
