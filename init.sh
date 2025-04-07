conda activate quasimodo && cd python_pkg &&\
alias python=python3\
export PYTHON_INCLUDE=`python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])"` \
export BOOST_PATH="~/stab_dd/boost_1_81_0"

invoke build-qreach && invoke build-pybind11
