# Copyright (C) 2018 OpenEye Scientific Software, Inc.
#
# THIS CODE IS PROPRIETARY TO OPENEYE SCIENTIFIC SOFTWARE INC AND IS SUBJECT
# TO THE FULL PROTECTION OF COPYRIGHT AND TRADESECRET LAW.  IT MAY BE USED
# ONLY PURSUANT TO A VALID AND CURRENT LICENSE FROM OPENEYE AND SUBJECT TO
# THE TERMS AND CONDITIONS THEREIN.  ALL OTHER USE IS STRICTLY PROHIBITED.
# PLEASE CONTACT OPENEYE AT LEGAL@EYESOPEN.COM IF YOU HAVE ANY QUESTIONS
# ABOUT THIS WARNING.

# !/usr/bin/env python
import re
import ast

from setuptools import setup, find_packages

try:  # for pip >= 10
    from pip._internal.req import parse_requirements
except ImportError:  # for pip <= 9.0.3
    from pip.req import parse_requirements


def get_reqs(reqs):
    return [str(ir.req) for ir in reqs]


try:
    install_reqs = get_reqs(parse_requirements("requirements.txt"))
except TypeError:
    try:  # For pip >= 10
        from pip._internal.download import PipSession
    except ImportError:  # For pip <= 9.0.3
        from pip.download import PipSession
    install_reqs = get_reqs(parse_requirements("requirements.txt", session=PipSession()))


def get_version():
    _version_re = re.compile(r'__version__\s+=\s+(.*)')
    with open('MDOrion/__init__.py', 'rb') as f:
        version = str(ast.literal_eval(_version_re.search(f.read().decode('utf-8')).group(1)))
        return version


setup(
    name="MDOrion",
    version='0.7.10',
    packages=find_packages(exclude=['tests*']),
    include_package_data=True,
    author="Christopher Bayly, Gaetano Calabro, Nathan M. Lim, John Chodera, ",
    author_email="bayly@eyesopen.com",
    description='Orion cubes to perform MD and MD analysis',
    install_requires=install_reqs,
    license='Other/Proprietary License',
    classifiers=[
        "Development Status :: 1 - Planning",
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Chemistry',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.5',
    ],
    zip_safe=False
)
