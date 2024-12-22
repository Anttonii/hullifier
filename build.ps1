md -Force lib
cd raylib/src
make
cp ./libraylib.a ../../lib
cd ../../
make