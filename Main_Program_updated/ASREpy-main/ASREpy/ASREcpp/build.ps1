Set-Location .\src\build\

cmake .. -G "Visual Studio 17 2022" -A x64

cmake --build . --config Release -j 4

cmake --install .

Set-Location ..\..