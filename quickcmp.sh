#!/bin/sh

# Get SDL
if [ ! -d sdljunk ]
then
	mkdir -p sdljunk
	
	# Download dev libs
	wget -c http://www.libsdl.org/release/SDL-devel-1.2.14-mingw32.tar.gz
	wget -c http://www.libsdl.org/projects/SDL_mixer/release/SDL_mixer-devel-1.2.12-VC.zip
	wget -c http://www.libsdl.org/projects/SDL_net/release/SDL_net-devel-1.2.8-VC.zip
	
	# Handle SDL
	cd sdljunk
	tar --strip-components=1 -xvvf ../SDL-devel-1.2.14-mingw32.tar.gz
	
	# Download Mixer
	unzip ../SDL_mixer-devel-1.2.12-VC.zip
	cp -r -f SDL_mixer-1.2.12/* .
	rm -rf SDL_mixer-1.2.12
	
	# Download Net
	unzip ../SDL_net-devel-1.2.8-VC.zip
	cp -r -f SDL_net-1.2.8/* .
	rm -rf SDL_net-1.2.8
	
	# Get out of here
	cd ..
fi

# Find the compiler game!
	# Debian Mingw32
if which i586-mingw32msvc-gcc
then
	HOSTIS="i586-mingw32msvc"
		
	# Mingw-w64 PC target (compat w/ 98 still; hopefully)
elif which i686-pc-mingw32-gcc
then
	HOSTIS="i686-pc-mingw32"
	
	# Mingw-w64 extended target (incompat with 98)
elif which i686-w64-mingw32-gcc
then
	HOSTIS="i686-w64-mingw32"
else
	HOSTIS=""
fi

# Recompile
#  --build=`uname -m`-linux  is not needed
./autogen.sh --with-sdl-exec-prefix=sdljunk --host=$HOSTIS --disable-sdltest --with-sdl-prefix=sdljunk CFLAGS="-D__WIN32__ -D_WIN32 -lSDL -Isdljunk/include -Isdljunk/include/SDL -I../sdljunk/include -I../sdljunk/include/SDL -I../../sdljunk/include -I../../sdljunk/include/SDL" LDFLAGS="-Lsdljunk/lib -Lsdljunk/lib/x86 -L../sdljunk/lib -L../sdljunk/lib/x86  -L../../sdljunk/lib -L../../sdljunk/lib/x86 -lSDL"
make clean

sed 's/-mwindows/-mconsole/g' < Makefile > $$
sed 's/-lSDLmain//g' < $$ > Makefile
sed 's/-lSDL/zzzzxcvbn/g' < Makefile > $$
sed 's/zzzzxcvbn/-lmingw32 -lSDLmain -lSDL/g' < $$ > Makefile

sed 's/-mwindows/-mconsole/g' < src/Makefile > $$
sed 's/-lSDLmain//g' < $$ > src/Makefile
sed 's/-lSDL/zzzzxcvbn/g' < src/Makefile > $$
sed 's/zzzzxcvbn/-lmingw32 -lSDLmain -lSDL/g' < $$ > src/Makefile

sed 's/-mwindows/-mconsole/g' < setup/Makefile > $$
sed 's/-lSDLmain//g' < $$ > setup/Makefile
sed 's/-lSDL/zzzzxcvbn/g' < setup/Makefile > $$
sed 's/zzzzxcvbn/-lmingw32 -lSDLmain -lSDL/g' < $$ > setup/Makefile

make -j2

# Make bin dir
rm -rf chocorenderlimits-win32
mkdir -p chocorenderlimits-win32

# Copy stuff there
cp src/chocolate-doom.exe chocorenderlimits-win32
cp src/chocolate-server.exe chocorenderlimits-win32
cp setup/chocolate-setup.exe chocorenderlimits-win32
cp sdljunk/bin/*.dll chocorenderlimits-win32
cp sdljunk/lib/x86/*.dll chocorenderlimits-win32
cp RENDERLIMITS chocorenderlimits-win32
cp README chocorenderlimits-win32
cp INSTALL chocorenderlimits-win32
cp COPYING chocorenderlimits-win32

# ZIP it up
zip chocorenderlimits-win32.zip chocorenderlimits-win32/*

#./configure --host=i586-mingw32msvc --build=`uname -m`-linux

