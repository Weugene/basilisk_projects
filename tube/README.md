# How to use Taylor flow project

## Python + Paraview
If you want to use Paraview pvpython, then:
1. install the same `python` version as `pvpython`.
2. install virtual environment: `/path/to/installed/python -m venv .venv`
3. before activation of `.venv` add this to the activate script `.venv/bin/activate`:
`export PYTHONPATH=$(python -c "import site; print(site.getsitepackages()[0])")`
4. activate `source .venv/bin/activate`
5. install packages to the current venv: `pip install -r requirements.txt`

