@echo off
if not "%minimized%"=="" goto :minimized
set minimized=true
start /min cmd /C "%~dpnx0"
exit

:minimized
chcp 65001
echo Running Python script...
python "%~dp0%DNA script.py"
if %ERRORLEVEL% NEQ 0 (
    echo Error detected, script did not run successfully.
    pause
) else (
    echo Script executed successfully.
    exit
)
