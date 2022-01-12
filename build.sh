#/bin/sh

CFLAGS="`pkg-config --cflags sdl2`"
LIBS="`pkg-config --libs sdl2`"

cc $CFLAGS -o particles main.c $LIBS
echo "DONE"
