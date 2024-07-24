@echo off
setlocal

set BATCH_DIR=%~dp0
chcp 65001
echo Checking for Python installation...
where python
if %ERRORLEVEL% NEQ 0 (
    echo Python is not installed or not added to the system PATH.
    pause
    exit
)

echo Checking for pip...
where pip
if %ERRORLEVEL% NEQ 0 (
    echo pip is not installed. Attempting to install pip...
    python -m ensurepip
    if %ERRORLEVEL% NEQ 0 (
        echo Failed to install pip.
        pause
        exit
    )
    echo pip installed successfully.
)

echo Updating pip...
python -m pip install --upgrade pip
if %ERRORLEVEL% NEQ 0 (
    echo Failed to update pip.
    pause
    exit
)

echo Installing required packages...
pip install -r "%BATCH_DIR%requirements.txt"
if %ERRORLEVEL% NEQ 0 (
    echo Failed to install required packages.
    pause
    exit
)

echo Updating or creating shortcut with new icon...
%SystemRoot%\system32\WindowsPowerShell\v1.0\powershell.exe -Command "$ws = New-Object -ComObject WScript.Shell; $shortcutPath = '%BATCH_DIR%DNA analyzer application.lnk'; $s = $ws.CreateShortcut($shortcutPath); $s.TargetPath = '%BATCH_DIR%DNA analyzer.bat'; $s.IconLocation = '%BATCH_DIR%DNA logo.ico'; $s.Save();"

echo Python installation and package requirements met. Running script...
python "%BATCH_DIR%DNA script.py"
if %ERRORLEVEL% NEQ 0 (
    echo Error detected, script did not run successfully.
    pause
) else (
    echo Script executed successfully.
    exit
)
