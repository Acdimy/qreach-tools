#! /bin/bash

cd python_pkg
mkdir output
echo "---Eval info:  grover   #Qubit:  7   Input dim:  1   Channel:  Unitary  ---" > output/result.txt
python test_new.py grover 4 1 >> output/result.txt
echo "---Eval info:  grover   #Qubit:  15   Input dim:  1   Channel:  Unitary  ---" >> output/result.txt
python test_new.py grover 8 1 >> output/result.txt
echo "---Eval info:  grover   #Qubit:  31   Input dim:  1   Channel:  Unitary  ---" >> output/result.txt
python test_new.py grover 16 1 >> output/result.txt
echo "---Eval info:  grover   #Qubit:  63   Input dim:  1   Channel:  Unitary  ---" >> output/result.txt
python test_new.py grover 32 1 >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  3   Input dim:  1   Channel:  Unitary  ---" >> output/result.txt
python test_new.py qrw 2 1 >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  5   Input dim:  1   Channel:  Unitary  ---" >> output/result.txt
python test_new.py qrw 3 1 >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  7   Input dim:  1   Channel:  Unitary  ---" >> output/result.txt
python test_new.py qrw 5 1 >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  9   Input dim:  2   Channel:  Unitary  ---" >> output/result.txt
python test_new.py qrw 7 2 >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  7   Input dim:  2   Channel:  Noise  ---" >> output/result.txt
python test_new.py qrw 5 2 ad >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  9   Input dim:  2   Channel:  Noise  ---" >> output/result.txt
python test_new.py qrw 7 2 ad >> output/result.txt
echo "---Eval info:  QRW   #Qubit:  10   Input dim:  2   Channel:  Noise  ---" >> output/result.txt
python test_new.py qrw 8 2 ad >> output/result.txt
echo "---Eval info:  RUS   #Qubit:  3   Input dim:  1   Channel:  Measure  ---" >> output/result.txt
python test_new.py rus 1 1 >> output/result.txt
echo "---Eval info:  RUS   #Qubit:  2   Input dim:  1   Channel:  Measure  ---" >> output/result.txt
python test_new.py rus 2 1 >> output/result.txt
echo "---Eval info:  RUS   #Qubit:  2   Input dim:  1   Channel:  Measure  ---" >> output/result.txt
python test_new.py rus 3 1 >> output/result.txt
sleep infinity