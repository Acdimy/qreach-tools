conda activate qreach && cd python_pkg &&\
export PYTHON_INCLUDE=`python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])"` \
export BOOST_PATH="~/boost_1_81_0"

invoke build-qreach && invoke build-pybind11
