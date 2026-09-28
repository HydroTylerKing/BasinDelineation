@echo off
setlocal
rem Basin Delineator - Windows launcher. Double-click to start.
cd /d "%~dp0"
set "VPY=.venv\Scripts\python.exe"
if exist "%VPY%" (
  "%VPY%" -c "import sys" >nul 2>&1 && goto run
  echo Existing environment is from another location - rebuilding ...
  rmdir /s /q .venv
  if exist "tools\python" rmdir /s /q "tools\python"
)

echo First run: setting up (one time only) ...
echo.

rem --- 1) Use an installed Python 3.10+ if there is one --------------------------------
set "PYEXE="
py -3 -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>&1 && set "PYEXE=py -3"
if not defined PYEXE python -c "import sys; sys.exit(0 if sys.version_info >= (3, 10) else 1)" >nul 2>&1 && set "PYEXE=python"
if not defined PYEXE goto private

echo Using the Python installed on this computer.
%PYEXE% -m venv .venv
if not exist "%VPY%" goto private
"%VPY%" -m pip install --upgrade pip
"%VPY%" -m pip install -r requirements.txt || goto fail
goto run

rem --- 2) No Python: get a private copy inside this folder (nothing installed system-wide)
:private
echo No Python found - downloading a private copy for this app (about 60 MB) ...
if not exist "tools\uv.exe" (
  powershell -NoProfile -ExecutionPolicy Bypass -Command "$ErrorActionPreference='Stop'; $ProgressPreference='SilentlyContinue'; [Net.ServicePointManager]::SecurityProtocol='Tls12'; Invoke-WebRequest 'https://github.com/astral-sh/uv/releases/latest/download/uv-x86_64-pc-windows-msvc.zip' -OutFile 'uv.zip'; Expand-Archive 'uv.zip' -DestinationPath 'tools' -Force; Remove-Item 'uv.zip'"
)
if not exist "tools\uv.exe" goto fail
set "UV_PYTHON_INSTALL_DIR=%~dp0tools\python"
set "UV_CACHE_DIR=%~dp0tools\cache"
if exist .venv rmdir /s /q .venv
"tools\uv.exe" venv .venv --python 3.12 || goto fail
echo Installing libraries ...
"tools\uv.exe" pip install --python "%VPY%" -r requirements.txt || goto fail

:run
"%VPY%" basin_app.py %*
pause
exit /b

:fail
echo.
echo Setup did not finish. Check your internet connection and double-click the launcher again.
echo If it keeps failing, install Python 3.12 from https://www.python.org/downloads/
echo (tick "Add python.exe to PATH" in the installer), then run the launcher again.
if exist .venv rmdir /s /q .venv
pause
exit /b 1
