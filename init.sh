conda activate qmc && cd python_pkg &&\
alias python=python3\
export PYTHON_INCLUDE=`python -c "from sysconfig import get_paths as gp; print(gp()[\"include\"])"` \
export BOOST_PATH="~/stab_dd/boost_1_81_0"
export NUSMV_LIBRARY_PATH="~/stab_dd/NuSMV-2.7.0/share/nusmv"
export PATH="~/stab_dd/NuSMV-2.7.0-linux64/bin:$PATH"

invoke build-qreach && invoke build-pybind11


g++ -g -O3 -std=c++20 -w -I. -I$BOOST_PATH -I.cflobdd/CFLOBDD -I.cflobdd/CFLOBDD/Solver/uwr/bit_vector/ -I.cflobdd/CFLOBDD/Solver/uwr/assert/ -I.cflobdd/CFLOBDD/Solver/uwr/matrix/ -I.cflobdd/CFLOBDD/Solver/uwr/parsing/ -lm\
 *.cpp cflobdd/CFLOBDD/Solver/uwr/bit_vector/*.cpp cflobdd/CFLOBDD/Solver/uwr/parsing/*.cpp cflobdd/CFLOBDD/*.cpp -o test
