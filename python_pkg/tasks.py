import invoke
import os
import platform
import shlex
import subprocess
import sys
import sysconfig
from pathlib import Path


def print_banner(msg):
    print("==================================================")
    print("= {} ".format(msg))


def _python_executable():
    return os.environ.get("PYTHON", sys.executable)


def _compiler():
    return os.environ.get("CXX", os.environ.get("CC", "g++"))


def _boost_path():
    return os.environ.get(
        "BOOST_PATH",
        str((Path(__file__).resolve().parents[2] / "BOOST" / "boost_1_81_0").resolve()),
    )


def _python_include():
    return os.environ.get("PYTHON_INCLUDE", sysconfig.get_paths()["include"])


def _extension_suffix():
    return sysconfig.get_config_var("EXT_SUFFIX") or ".so"


def _pybind11_includes():
    py = _python_executable()
    return subprocess.check_output(
        [py, "-m", "pybind11", "--includes"],
        text=True,
    ).strip()


def _shared_link_flags():
    if platform.system() == "Darwin":
        return "-L. -lqreach -Wl,-rpath,@loader_path -undefined dynamic_lookup"
    return "-L. -lqreach -Wl,-rpath,."


@invoke.task()
def build_qreach(c):
    """Build the shared library for the C++ code."""
    print_banner("Building C++ Library")
    c.run(
        "cd .. && "
        f"BOOST_PATH={shlex.quote(_boost_path())} "
        f"CXX={shlex.quote(_compiler())} "
        "make all && "
        "cd python_pkg/ && cp ../libqreach.so ."
    )
    print("* Complete")


def compile_python_module(cpp_name, extension_name):
    cmd = " ".join([
        shlex.quote(_compiler()),
        "-g -O3 -std=c++2a -w -shared -Wall -Wextra -DHAVE_CONFIG_H",
        "-Werror -Wunused-but-set-variable -fPIC",
        f"-I{shlex.quote(_python_include())}",
        _pybind11_includes(),
        f"-I{shlex.quote(_boost_path())}",
        "-I../",
        shlex.quote(cpp_name),
        f"-o {shlex.quote(extension_name + _extension_suffix())}",
        _shared_link_flags(),
    ])
    invoke.run(cmd)


@invoke.task()
def clean_qreach(c):
    print_banner("Clean qreach")
    c.run(
        "cd .. && "
        f"BOOST_PATH={shlex.quote(_boost_path())} "
        f"CXX={shlex.quote(_compiler())} "
        "make clean && cd python_pkg/"
    )
    print("* Complete")


@invoke.task()
def build_pybind11(c):
    """Build the pybind11 wrapper library."""
    print_banner("Building PyBind11 Module")
    compile_python_module("qreach_python_wrapper.cpp", "pyqreach")
    print("* Complete")


@invoke.task()
def test_pybind11(c):
    """Run the script to test PyBind11."""
    print_banner("Testing PyBind11 Module")
    c.run(f"{shlex.quote(_python_executable())} pybind11_test.py", pty=True)
