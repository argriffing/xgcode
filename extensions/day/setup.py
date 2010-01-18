from distutils.core import setup, Extension

if __name__ == '__main__':
    setup(name="day", version="1.0", ext_modules=[Extension("day", ["daymodule.c"])])

