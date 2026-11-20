@echo off
REM ============================================================
REM  WDM_and_HQI : run scripts 01 to 11 on Windows (FULL=0 by default)
REM  Put this file in the  code\  folder and run it. Nothing to edit:
REM  folders are found from the location of this file.
REM  Each script writes results\log_XX.txt; a one-line verdict per script
REM  (exit code, expected output file, elapsed time) goes to results\run_summary.txt.
REM  Scripts are never skipped: a failure is recorded and the next one runs.
REM  99_verify_all.R is not part of this run (see README).
REM ============================================================
setlocal EnableDelayedExpansion

REM --- folders, derived from where this .bat sits ---
set "CODE_DIR=%~dp0"
set "CODE_DIR=%CODE_DIR:~0,-1%"
for %%I in ("%CODE_DIR%\..") do set "ROOT=%%~fI"
set "LOCAL_DATA=%ROOT%\data"
set "LOCAL_OUT=%ROOT%\results"
set "FULL=0"
if not exist "%LOCAL_OUT%" mkdir "%LOCAL_OUT%"

REM --- find Rscript.exe ---
set "RS="
for /f "delims=" %%P in ('where Rscript.exe 2^>nul') do if not defined RS set "RS=%%P"
if not defined RS for /d %%D in ("%ProgramFiles%\R\R-*") do if exist "%%D\bin\Rscript.exe" set "RS=%%D\bin\Rscript.exe"
if not defined RS (
  echo Rscript.exe not found. Install R from https://cran.r-project.org/ and run again.
  pause
  exit /b 1
)

REM R wants forward slashes
set "CODE_DIR=%CODE_DIR:\=/%"
set "LOCAL_DATA=%LOCAL_DATA:\=/%"
set "LOCAL_OUT=%LOCAL_OUT:\=/%"
set "OUT_DIR=%LOCAL_OUT%"
set "SUMMARY=%LOCAL_OUT%/run_summary.txt"
set /a NFAIL=0
set "FAILED_LIST="

echo Rscript  : %RS%
echo code     : %CODE_DIR%
echo data     : %LOCAL_DATA%
echo results  : %LOCAL_OUT%
echo.
cd /d "%~dp0"
> "%SUMMARY%" echo run started %date% %time%  (FULL=%FULL%)
>> "%SUMMARY%" echo script                                  exit  output  elapsed
>> "%SUMMARY%" echo --------------------------------------  ----  ------  -------

REM             script                              expected output file (one per script)
call :run 01_exact_simul_figures.R          ExtDataTable1_exact_moments.csv
call :run 02_closed_forms.R                 closed_form_check.csv
call :run 03_supp_table1.R                  SuppTable1_simulation_vs_exact.csv
call :run 04_realdata.R                     realdata_supercond_rank.csv
call :run 05_table2_extdata_table2.R        Fig4.png
call :run 06_supp_table2_forest.R           SuppTable2_forest_variance.csv
call :run 07_supp_table3_copula.R           SuppTable3_copula_entry.csv
call :run 08_supp_tables4_5_split_rules.R   SuppTable5_occupancy.csv
call :run 09_supp_table6_zeta_limits.R      SuppTable6_zeta_limits_n4096.csv
call :run 10_supp_tables7_9_maximal_state.R SuppTable9_maximal_state.csv
call :run 11_supp_note11_third_moment.R     third_moment_check.csv

echo.
if %NFAIL%==0 (
  >> "%SUMMARY%" echo run finished %date% %time%: all 11 scripts ok
  echo All 11 scripts finished without error. Results are in %LOCAL_OUT%
) else (
  >> "%SUMMARY%" echo run finished %date% %time%: %NFAIL% problem(s) in%FAILED_LIST%
  echo Finished with %NFAIL% problem(s) in%FAILED_LIST%
  echo See %SUMMARY% and the log files.
)
echo (summary: %SUMMARY%)
pause
exit /b %NFAIL%

:run
set "S=%~1"
set "EXPECT=%~2"
set "L=%LOCAL_OUT%/log_%S:~0,-2%.txt"
set "T0=%time%"
echo [%time%] %S% ...
"%RS%" "%S%" > "%L%" 2>&1
set "RC=%errorlevel%"
>> "%L%" echo.
>> "%L%" echo === exit code %RC% ===
if exist "%LOCAL_OUT%/%EXPECT%" (set "OUTOK=ok") else (set "OUTOK=MISSING")
call :elapsed "%T0%"
set "LINE=%S%"
set "LINE=!LINE!                                        "
set "LINE=!LINE:~0,38!"
if "%RC%"=="0" (set "RCS=0   ") else (set "RCS=%RC%   ")
set "RCS=!RCS:~0,4!"
>> "%SUMMARY%" echo !LINE!  !RCS!  %OUTOK%     %ELAPSED%
if "%RC%"=="0" if "%OUTOK%"=="ok" (
  echo    ok  ^(%ELAPSED%^)
  exit /b 0
)
set /a NFAIL+=1
set "FAILED_LIST=%FAILED_LIST% %S:~0,2%"
if not "%RC%"=="0" (echo    FAILED - exit code %RC%, see %L%) else (echo    WARNING - exit code 0 but %EXPECT% not written, see %L%)
exit /b 0

:elapsed
set "A=%~1"
set "A=%A: =0%"
set "B=%time%"
set "B=%B: =0%"
for /f "tokens=1-4 delims=:.," %%a in ("%A%") do set /a S0=(1%%a-100)*3600+(1%%b-100)*60+(1%%c-100)
for /f "tokens=1-4 delims=:.," %%a in ("%B%") do set /a S1=(1%%a-100)*3600+(1%%b-100)*60+(1%%c-100)
set /a DT=S1-S0
if %DT% lss 0 set /a DT+=86400
set /a MM=DT/60, SS=DT%%60
if %SS% lss 10 (set "ELAPSED=%MM%m0%SS%s") else (set "ELAPSED=%MM%m%SS%s")
exit /b 0
