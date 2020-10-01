# first install all dependencies
pip install -r requirements.txt

# build all the cython code
cd src
python3 setup.py build_ext --inplace