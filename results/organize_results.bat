@echo off
setlocal
cd /d "%~dp0"
for %%D in (figures tables source_data checks intermediate logs) do if not exist "%%D" mkdir "%%D"
move /y Fig*.pdf figures\ >nul 2>&1
move /y Fig*.png figures\ >nul 2>&1
move /y ExtData_Fig*.pdf figures\ >nul 2>&1
move /y ExtData_Fig*.png figures\ >nul 2>&1
move /y Table*.csv tables\ >nul 2>&1
move /y ExtDataTable*.csv tables\ >nul 2>&1
move /y SuppTable*.csv tables\ >nul 2>&1
move /y closed_form_check.csv checks\ >nul 2>&1
move /y verify_closed_forms.csv checks\ >nul 2>&1
move /y third_moment_check.csv checks\ >nul 2>&1
move /y realdata_*.csv intermediate\ >nul 2>&1
move /y copula_*.csv intermediate\ >nul 2>&1
move /y copula_*.rds intermediate\ >nul 2>&1
move /y sim_ch6.rds intermediate\ >nul 2>&1
move /y log_*.txt logs\ >nul 2>&1
move /y run_summary.txt logs\ >nul 2>&1
echo done.
dir /b
pause