@echo off
set DEPENDENCIES_PREFIX=%cd%\dependencies
set BOOST_PREFIX=%DEPENDENCIES_PREFIX%\boost
set INSTALL_PREFIX=%DEPENDENCIES_PREFIX%\install
call "C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\vcvarsall.bat" amd64


if exist "%INSTALL_PREFIX%" (goto done)


pushd "%BOOST_PREFIX%"
call bootstrap.bat
b2 headers
set BOOST_LIBRARIES=--with-filesystem --with-date_time
b2 --toolset=msvc-12.0 --prefix="%INSTALL_PREFIX%\debug" --layout=tagged %BOOST_LIBRARIES% address-model=64 link=static runtime-link=shared variant=debug install
b2 --toolset=msvc-12.0 --prefix="%INSTALL_PREFIX%\release" --layout=tagged %BOOST_LIBRARIES% address-model=64 link=static runtime-link=shared variant=release install
popd


:done
