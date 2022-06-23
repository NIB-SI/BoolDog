import os, tempfile, shutil

from pathlib import Path
from distutils.core import setup, Extension
import numpy

class temp_file():
    def __init__(self, file_name, method, folder):
        path = Path(folder.tmpdir, file_name)
        self.file_obj = open(path, method)
    def __enter__(self):
        return self.file_obj
    def __exit__(self, type, value, traceback):
        self.file_obj.close()

class temp_dir():

    # Ensure the file is read/write by the creator only
    saved_umask = os.umask(0077)

    def __init__():
        self.tmpdir = tempfile.mkdtemp()

    def create_file(file_name="test.txt"):
        path = Path(self.tmpdir, predictable_filename)
        print(path)
        try:
        with open(path, "w") as tmp:
            tmp.write("secrets!")
    except IOError as e:
        print 'IOError'
    else:
        os.remove(path)
    finally:
        os.umask(saved_umask)
        os.rmdir(tmpdir)




try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()


cIAA_func_module = Extension('_cIAA_func',
                             ['cIAA_func.c', 'cIAA_func_wrap.c'],
                             include_dirs=[numpy_include])

setup(name='cIAA_func',
      ext_modules = [cIAA_func_module],
      py_modules=["cIAA_func"],
      script_name = 'setup.py',
      script_args=["build_ext", "--inplace"]
)