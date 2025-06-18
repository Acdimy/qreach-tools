OPENQASM 2.0;
include "qelib1.inc";
qreg q[8];
creg c[8];
x q[0];
cx q[0], q[1];
measure q[0] -> c[1];
if (c[1] == 0) measure q[2] -> c[5];