# source this file: dependency prefix for LiPolGen (HepMC3 3.3.0, LHAPDF 6.5.5, PYTHIA 8.317)
export LIPOLGEN_DEPS=$(cd "$(dirname "${BASH_SOURCE[0]}")/../deps/install" && pwd)
export PATH=$LIPOLGEN_DEPS/bin:$PATH
export LD_LIBRARY_PATH=$LIPOLGEN_DEPS/lib:$LD_LIBRARY_PATH
export PYTHIA8DATA=$LIPOLGEN_DEPS/share/Pythia8/xmldoc
export LHAPDF_DATA_PATH=$LIPOLGEN_DEPS/share/LHAPDF:${LHAPDF_DATA_PATH}
export PYTHONPATH=$LIPOLGEN_DEPS/lib:$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/build/python:$PYTHONPATH
