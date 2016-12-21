#!/usr/bin/env python

from setuptools import setup, find_packages, Extension

extensions = [
    Extension('cyhll.hyperloglog', ['cyhll/hyperloglog.pyx']),
    Extension('cyhll.murmer3', ['cyhll/murmer3.pyx']),
]

try:
    from Cython.Build import cythonize
    extensions = cythonize(extensions,
                           compiler_directives={
                               'language_level': 3
                           })
except ImportError:
    pass


if __name__ == '__main__':
    setup(
        name='cyhll',
        version='0.1.4',
        url='https://github.com/bndl/cyhll',
        description='Hyperloglog in Cython',
        long_description=open('README.rst').read(),
        author='Frens Jan Rumph',
        author_email='mail@frensjan.nl',

        packages=(
            find_packages()
        ),

        include_package_data=True,
        zip_safe=False,

        install_requires=[],
        extras_require=dict(
            dev=[
                'cython<0.25',
                'pytest',
            ],
        ),

        ext_modules=extensions,

        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Developers',
            'Operating System :: OS Independent',
            'Programming Language :: Python',
            'Programming Language :: Python :: 3.3',
            'Programming Language :: Python :: 3.4',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
        ],
    )
