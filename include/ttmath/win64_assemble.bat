rem make sure this is a proper path to the 64 bit assembler
"C:\Program Files (x86)\Microsoft Visual Studio 14.0\VC\bin\amd64\ml64.exe" /c ttmathuint_x86_64_msvc.asm
rem ml64.exe will produce ttmathuint_x86_64_msvc.obj which should be added (linked) to your project

rem or you can assemble with debug info
rem ml64.exe /c /Zd /Zi ttmathuint_x86_64_msvc.asm

rem be nice, most Windows users just click on the file
pause
