set -e

cd /cygdrive/c/codes/hybird/submarine/HYBIRD/cmake-build-release
/cygdrive/c/Users/ci1apa/AppData/Local/JetBrains/CLion2024.2/cygwin_cmake/bin/cmake.exe --regenerate-during-build -S$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
