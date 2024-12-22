main:
	gcc -o hullifier.exe main.cpp -O2 -L./lib/ -I./raylib/src/ -I./raygui/src/ -lraylib -lgdi32 -lwinmm -lstdc++
