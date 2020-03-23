from distutils.core import setup
from os import path
this_directory = path.abspath(path.dirname(__file__))

with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

setup(
    name="cctk",
    packages=["cctk", "cctk.data", "cctk.groups"],
#    include_package_data=True,
    package_data={"cctk.data": ["*"], "cctk.groups": ["*"],},
    version="v0.1.3",
    license="Apache 2.O",
    description="computational chemistry toolkit",
    author="Corin Wagen and Eugene Kwan",
    author_email="corin.wagen@gmail.com",
    url="https://github.com/ekwan/cctk",
    download_url="https://github.com/ekwan/cctk/archive/v0.1.3.tar.gz",
    install_requires=["numpy", "networkx", "importlib_resources"],
    long_description=long_description,
    long_description_content_type='text/markdown',
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
