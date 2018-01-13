@ECHO OFF
ECHO Searching for install path of latest version of R in registry...
SETLOCAL ENABLEEXTENSIONS REM This line will reset path when return
SET RKEY=
SET RPATH=
FOR /F "tokens=* skip=2" %%L IN ('reg.exe QUERY HKLM\Software\R-core\R /f * /k ^| sort') DO (
    IF NOT "%%~L"=="" SET "RKEY=%%~L"
)
IF NOT DEFINED RKEY (
    ECHO Unable to query registry key HKLM\Software\R-core\R
    EXIT /B 1
)
FOR /F "tokens=2* skip=2" %%A IN ('REG QUERY %RKEY% /v "installPath"') DO (
    IF NOT "%%~B"=="" SET "RPATH=%%~B"
)
IF NOT DEFINED RPATH (
    ECHO Unable to query registry value %RKEY%\installPath
    EXIT /B 2
)
IF NOT EXIST "%RPATH%" (
    ECHO Found path for R (%RPATH%^) does not exist
    EXIT /B 3
)
SET OLDPATH=%PATH%
IF "%PROCESSOR_ARCHITECTURE%"=="AMD64" (
    SET PATH=%RPATH%\bin\x64;%OLDPATH%
    ECHO Found %RPATH%\bin\x64
) ELSE (
    SET PATH=%RPATH%\bin\i386;%OLDPATH%
    ECHO Found %RPATH%\bin\i386
)

echo %cd%
pushd %~dp0
echo %cd%


	git add --all
	git commit -m "Simulations added."
	git push 
	