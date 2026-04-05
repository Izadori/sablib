from conan import ConanFile
from conan.tools.cmake import CMake, cmake_layout
import os
import shutil


class SablibConan(ConanFile):
    name = "sablib"
    version = "0.3.2"
    package_type = "library"

    license = "MIT"
    url = "https://github.com/Izadori/sablib"
    description = "Signal smoothing and baseline estimation library"
    topics = ("signal-processing", "smoothing", "baseline")
    package_type = "static-library"
    settings = "os", "compiler", "build_type", "arch"

    exports_sources = (
        "CMakeLists.txt",
        "baseline/*",
        "smoothing/*",
        "misc/*",
        "sablib.h",
        "cmake/*"
    )

    generators = "CMakeToolchain", "CMakeDeps"

    def requirements(self):
        self.requires("eigen/3.4.0")

    def layout(self):
        cmake_layout(self, src_folder="../..")  # リポジトリルートを直接参照

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["sablib"]
