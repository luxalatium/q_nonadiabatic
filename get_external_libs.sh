# Download required external libraries

mkdir -p libraries

# EIGEN (https://eigen.tuxfamily.org/)
wget https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz
tar -xf eigen-3.4.0.tar.gz
mv eigen-3.4.0/Eigen libraries 
rm -r eigen-3.4.0 eigen-3.4.0.tar.gz

# PCG (https://www.pcg-random.org/)
mkdir libraries/PCG
wget https://www.pcg-random.org/downloads/pcg-cpp-0.98.zip
unzip pcg-cpp-0.98.zip
mv pcg-cpp-0.98/include libraries/PCG/pcg-cpp-0.98
rm -r pcg-cpp-0.98 pcg-cpp-0.98.zip
