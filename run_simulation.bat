@echo off
echo FDTD Laser Simulation Runner
echo =============================

echo.
echo 1. Compiling C++ simulation...
g++ -O3 -std=c++11 -o Laser.exe Laser.cpp
if %errorlevel% neq 0 (
    echo Error: Compilation failed!
    pause
    exit /b 1
)
echo Compilation successful!

echo.
echo 2. Running simulation...
echo This will take a few minutes...
Laser.exe
if %errorlevel% neq 0 (
    echo Error: Simulation failed!
    pause
    exit /b 1
)

echo.
echo 3. Simulation complete! Files generated:
dir /B *.dat | findstr "field"
echo.
echo 4. Running Python visualization...
python quick_visualize.py

echo.
echo 5. All done! Check the generated PNG files for results.
pause
