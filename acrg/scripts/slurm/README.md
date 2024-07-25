To run this script, you must be in an environment with pandas available.

To see what options are available, call
``` sh
python make_slurm_array.py -h
```

Using the following options
``` sh
python make_slurm_array.py -J "slurm_script_test" --ini-path "/user/work/bm13805/ini_files" --ini-file "test.ini" --job-root "/user/work/bm13805/" --python-venv "/user/work/bm13805/.venv_inversions_standard" --year 2018 --nmonths 48 --kwargs '{"min_model_error": 20.0}' --ncpu 4 --mem "40gb" --time "12:00:00"
```
creates the files `config.txt` and `inversion_wrapper.sh`, which are present in this repo.

They are created in the folder `/user/work/bm13805/slurm_script_test`, which is the `--job-root` argument followed by the `-J` (job name) argument.

This script logs the git branch and commit used, as well all packages (with versions) from the virtual environment.
