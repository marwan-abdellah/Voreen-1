@echo off
rem Uncomment, if qt is not in your path
rem set path=c:\qt\4.3.4\bin

rem Select your Visual Studio version
rem set QMAKESPEC=win32-msvc.net
rem set QMAKESPEC=win32-msvc2005
set set QMAKESPEC=win32-msvc2008

cd src\core
qmake.exe -tp vc "CONFIG+=flat voreenlib_core.pro
cd ..\qt
qmake.exe -tp vc "CONFIG+=flat" voreenlib_qt.pro
cd ..\..\apps\rptgui
qmake.exe -tp vc "CONFIG+=flat" rptgui.pro

rem Optional
rem cd ..\..\apps\voreendev
rem qmake.exe -tp vc "CONFIG+=flat" voreendev.pro
rem cd ..\..\apps\voreenpub
rem qmake.exe -tp vc "CONFIG+=flat" voreenpub.pro
rem cd ..\..\apps\voltool
rem qmake.exe -tp vc "CONFIG+=flat" voltool.pro
rem cd ..\..\apps\voreenus
rem qmake.exe -tp vc "CONFIG+=flat" voreenus.pro
rem cd ..\..\apps\glyphapp
rem qmake.exe -tp vc "CONFIG+=flat" glyphapp.pro

cd ..\..