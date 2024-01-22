:: for %%a in (shell_scripts/*) do echo %%a

cd shell_scripts
::dir /b *.sh

for %%a in (*.sh) do (
	wsl bash %%a
	)

cmd /k