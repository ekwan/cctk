from distutils.core import setup

setup(
    name="cctk",
    packages=["cctk"],
    version="v0.1.0",
    license="Apache 2.O",
    description="computational chemistry toolkit",
    author="Corin Wagen and Eugene Kwan",
    author_email="corin.wagen@gmail.com",
    url="https://github.com/ekwan/cctk",
    download_url="https://github.com/ekwan/cctk/archive/v0.1.0.tar.gz",
    install_requires=["numpy", "networkx",],
    classifiers=[
        "Development Status :: 4 - Beta",
        "License :: OSI Approved :: Apache Software License",
        "Intended Audience :: Science/Research",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
)
