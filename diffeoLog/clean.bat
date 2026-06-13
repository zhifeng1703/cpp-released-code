@echo off
echo Deleting all .mexw64 files in the current directory and subdirectories...

for /r %%f in (*.mexw64 ) do (
    echo Deleting %%f
    del "%%f"
)
echo Done.
pause