import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="cctk",
    version="0.0.1",
    author="Eugene Kwan and Corin Wagen",
    author_email="author@example.com",
    description="Computational Chemistry Tool Kit",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ekwan/cctk",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
